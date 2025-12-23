module Storage

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP, DataFrames, CSV

# Use internal modules
using ..Utils

# Export variables and functions
export StorageUnit, load_data, stochastic_capex_model!, toCSV_stochastic_capex

"""
Storage unit represents an energy storage system in the power system.
# Fields:
- id: ID of the storage system
- storage: name of storage unit. It could be any string.
- technology: technology type of the storage system, for example, "lithium-ion", "flow". It could be any string.
- bus: name of the bus where the storage system is connected to.
- invest_cost: investment cost per MW of power capacity (USD/MW)
- exist_power_cap: pre-existing power capacity of the storage system (MW)
- exist_energy_cap: pre-existing energy capacity of the storage system (MWh)
- var_om_cost: variable operation and maintenance cost (USD/MW)
- efficiency: round-trip efficiency of the storage system (between 0 and 1)
- duration: duration of the storage system at full power (hours)
"""
struct StorageUnit
    id:: Int64
    storage:: String
    technology:: String
    bus:: String
    cap_existing_energy_MWh:: Float64
    cap_existing_power_MW:: Float64
    cap_max_power_MW:: Float64
    cost_fixed_energy_USDperkWh:: Float64
    cost_fixed_power_USDperkW:: Float64
    cost_variable_USDperMWh:: Float64
    duration_hr:: Float64
    efficiency_charge:: Float64
    efficiency_discharge:: Float64
    c0_USD:: Float64
    c1_USDperMWh:: Float64
    c2_USDperMWh2:: Float64
end

function load_data(inputs_dir::String)

    # Load energy storage units using CSVs
    filename = "storage.csv"
    print(" > $filename ...")
    E = to_structs(StorageUnit, joinpath(inputs_dir, filename))
    println(" ok, loaded ", length(E), " storage units.")

    return E
end

function stochastic_capex_model!(sys, mod:: Model)

    N = @views sys.N
    T = @views sys.T
    S = @views sys.S
    E = @views sys.E

    # Filter energy storage units by bus
    E_AT_BUS = [filter(e -> e.bus == n.bus, E) for n in N]

    # Define generation variables
    # vCHARGE: power discharged from the storage system (MW)
    # vDISCHA: power discharged from the storage system (MW)
    # vSOC: state of charge of the storage system (MWh)
    # vPCAP: power capacity of the storage system (MW)
    # vECAP: energy capacity of the storage system (MWh)
    @variables(mod, begin
            vDISCHA[E, S, T] ≥ 0
            vCHARGE[E, S, T] ≥ 0       
            vSOC[E, S, T] ≥ 0
            vPCAP[E, S] ≥ 0  
            vECAP[E, S] ≥ 0  
    end)

    @constraint(mod, cMinEnerCapStor[e ∈ E, s ∈ S],     
                        vECAP[e, s] ≥ e.cap_existing_energy_MWh)
    
    @constraint(mod, cMinPowerCapStor[e ∈ E, s ∈ S], 
                        vPCAP[e, s] ≥ e.cap_existing_power_MW)

    @constraint(mod, cMaxPowerCapStor[e ∈ E, s ∈ S], 
                        vPCAP[e, s] ≤ e.cap_max_power_MW)

    # Define constraints for energy storage systems
    @constraint(mod, cMaxCharge[e ∈ E, s ∈ S, t ∈ T], 
                        vCHARGE[e, s, t] ≤ vPCAP[e, s])
    
    @constraint(mod, cMaxDischa[e ∈ E, s ∈ S, t ∈ T], 
                        vDISCHA[e, s, t] ≤ vPCAP[e, s])
    
    E_fixduration = filter(e -> e.duration_hr > 0, E)

    if !isempty(E_fixduration)
    @constraint(mod, cFixEnergyPowerRatio[e ∈ E_fixduration, s ∈ S], 
                        vECAP[e, s] ==  e.duration_hr * vPCAP[e, s] )
    end

    @constraint(mod, cMaxSOC[e ∈ E, s ∈ S, t ∈ T], 
                        vSOC[e, s, t] ≤ vECAP[e, s])

    # SOC in the next time is a function of SOC in the previous time
    # with circular wrapping for the first and last timepoints within a timeseries
    @constraint(mod, cStateOfCharge[e ∈ E, s ∈ S, t ∈ T],
                        vSOC[e, s, t] == vSOC[e, s, T[t.prev_timepoint_id]] +
                                        t.duration_hr*(vCHARGE[e, s, t]*e.efficiency_charge 
                                                        - vDISCHA[e, s, t]*1/e.efficiency_discharge) )

    # Power generation by bus
    @expression(mod, eNetDischargeAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    - sum(vCHARGE[e, s, t] for e ∈ E_AT_BUS[n.id]) 
                    + sum(vDISCHA[e, s, t] for e ∈ E_AT_BUS[n.id]) )

    # Storage cost per timepoint
    @expression(mod, eStorCostPerTp[t ∈ T],
                     1/length(S)*(sum(s.probability * e.cost_variable_USDperMWh * vCHARGE[e, s, t] for e ∈ E, s ∈ S) ) )
    
    eCostPerTp =  @views mod[:eCostPerTp]
    unregister(mod, :eCostPerTp)
    @expression(mod, eCostPerTp[t ∈ T], eCostPerTp[t] + eStorCostPerTp[t])

                     
    # Storage cost per period
	@expression(mod, eStorCostPerPeriod,
                    1/length(S)*sum( s.probability * (e.cost_fixed_power_USDperkW * vPCAP[e, s] * 1000 
                                                + e.cost_fixed_energy_USDperkWh * vECAP[e, s] * 1000) for e ∈ E, s ∈ S ))

    eCostPerPeriod =  @views mod[:eCostPerPeriod]
    unregister(mod, :eCostPerPeriod)
    @expression(mod, eCostPerPeriod, eCostPerPeriod + eStorCostPerPeriod)
                

    # Total costs
    @expression(mod, eStorTotalCost, sum(eStorCostPerTp[t] * t.weight for t in T) + eStorCostPerPeriod)
    
end

function toCSV_stochastic_capex(sys, mod:: Model, outputs_dir:: String)

    # Print vCHARGE AND vDISCHARGE 
    df1 = to_df(mod[:vCHARGE], [:storage, :scenario, :timepoint, :Charge_MW])

    df2 = to_df(mod[:vDISCHA], [:storage, :scenario, :timepoint, :Discharge_MW])
    
    df_mix1 = outerjoin(df1, df2, on=[:storage, :scenario, :timepoint])
    CSV.write(joinpath(outputs_dir, "storage_dispatch.csv"), df_mix1)

    to_df(mod[:vSOC], [:storage, :scenario, :timepoint, :SOC_MWh];
                    csv_dir = joinpath(outputs_dir, "storage_soc.csv"))
    
    # Print vGENV variable solution
    df1 = to_df(mod[:vPCAP], [:storage, :scenario, :PowerCap_MW])

    # Print vCAPV variable solution
    df2 = to_df(mod[:vECAP], [:storage, :scenario, :EnergyCap_MWh])

    df_mix1 = outerjoin(df1, df2, on=[:storage, :scenario])
    CSV.write(joinpath(outputs_dir, "storage_capacity.csv"), df_mix1)

    # Print cost expressions
    filename = "storage_costs_itemized.csv"
    costs =  DataFrame(component  = ["CostPerTimepoint", "CostPerPeriod", "TotalCost"], 
                            cost  = [   value(sum(mod[:eStorCostPerTp][t] * t.weight for t in sys.T)), 
                                        value(mod[:eStorCostPerPeriod]), 
                                        value(mod[:eStorTotalCost])]) 
    CSV.write(joinpath(outputs_dir, filename), costs)
    println(" > $filename printed.")

end


end # module EnergyStorage