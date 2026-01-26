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
mutable struct StorageUnit
    id:: Int64
    name:: String
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
    expand_capacity:: Bool
end

function StorageUnit(id, storage, technology, bus, cap_existing_energy_MWh, cap_existing_power_MW, cap_max_power_MW, cost_fixed_energy_USDperkWh, cost_fixed_power_USDperkW, cost_variable_USDperMWh, duration_hr, efficiency_charge, efficiency_discharge, c0_USD, c1_USDperMWh, c2_USDperMWh2; expand_capacity = true)
    return StorageUnit(id, storage, technology, bus, cap_existing_energy_MWh, cap_existing_power_MW, cap_max_power_MW, cost_fixed_energy_USDperkWh, cost_fixed_power_USDperkW, cost_variable_USDperMWh, duration_hr, efficiency_charge, efficiency_discharge, c0_USD, c1_USDperMWh, c2_USDperMWh2, expand_capacity)
    
end

function load_data(inputs_dir::String)

    # Load energy storage units using CSVs
    start_time = time() 
    filename = "storage.csv"
    println(" > $filename ...")
    E = to_structs(StorageUnit, joinpath(inputs_dir, filename))
    println(" └ Completed, loaded ", length(E), " storage units. Elapsed time ", round(time() - start_time, digits=2), " seconds.")

    # Extra calculations or checks can be added here
    for e in E
        if e.cap_existing_power_MW >= e.cap_max_power_MW
            e.expand_capacity = false
        end
    end

    return E
end

function stochastic_capex_model!(sys, mod:: Model, model_settings::Dict)

    N = @views sys.N
    T = @views sys.T
    S = @views sys.S
    E = @views sys.E

    # Filter energy storage units by bus
    E_AT_BUS = [filter(e -> e.bus == n.name, E) for n in N]

    # Define generation variables
    # vCHARGE: power discharged from the storage system (MW)
    # vDISCHA: power discharged from the storage system (MW)
    # vSOC: state of charge of the storage system (MWh)
    # vPCAP: power capacity of the storage system (MW)
    # vECAP: energy capacity of the storage system (MWh)    
    @variable(mod, vSOC[E, S, T] ≥ 0)
    @variable(mod, vPCAP[E, S] ≥ 0)  
    @variable(mod, vECAP[E, S] ≥ 0)  

    if model_settings["consider_single_storage_injection"] == true
            @variable(mod, vDISCHA[E, S, T])
    else
            @variable(mod, vDISCHA[E, S, T] ≥ 0)
            @variable(mod, vCHARGE[E, S, T] ≥ 0)       
    end

    for e in E
        if e.expand_capacity == false
            for s in S
                fix(vPCAP[e, s], 0.0; force=true)
                fix(vECAP[e, s], 0.0; force=true)
            end
        end
    end

    println(" - Constraints of power capacity expansion")
    @constraint(mod, cMaxPowerCapStor[e ∈ E, s ∈ S; e.expand_capacity == true], 
                        vPCAP[e, s] ≤ e.cap_max_power_MW - e.cap_existing_power_MW)


    println(" - Maximum charge and discharge constraints")
    if model_settings["single_storage_injection"] == true
        @constraint(mod, cMaxCharge[e ∈ E, s ∈ S, t ∈ T], 
                                          vDISCHA[e, s, t] >= - (vPCAP[e, s] + e.cap_existing_power_MW))
        @constraint(mod, cMaxDischa[e ∈ E, s ∈ S, t ∈ T], 
                                          vDISCHA[e, s, t] <= vPCAP[e, s] + e.cap_existing_power_MW)
    else
        # Define constraints for energy storage systems
        @constraint(mod, cMaxCharge[e ∈ E, s ∈ S, t ∈ T; e.expand_capacity == true], 
                        vCHARGE[e, s, t] ≤ vPCAP[e, s] + e.cap_existing_power_MW)
    
        @constraint(mod, cMaxDischa[e ∈ E, s ∈ S, t ∈ T; e.expand_capacity == true], 
                        vDISCHA[e, s, t] ≤ vPCAP[e, s] + e.cap_existing_power_MW)
    end
    
    E_fixduration_expandable = filter(e -> ((e.duration_hr > 0)  && (e.expand_capacity == true)), E)

    if !isempty(E_fixduration_expandable)
        println(" - Energy-power ratio constraints for expandable storage with fixed duration")
    @constraint(mod, cFixEnergyPowerRatio[e ∈ E_fixduration_expandable, s ∈ S], 
                        (vECAP[e, s] + e.cap_existing_energy_MWh) ==  e.duration_hr * (vPCAP[e, s] + e.cap_existing_power_MW) )
    end

    println(" - Maximum state of charge constraints")
    @constraint(mod, cMaxSOC[e ∈ E, s ∈ S, t ∈ T], 
                        vSOC[e, s, t] ≤ (vECAP[e, s] + e.cap_existing_energy_MWh))

    # SOC in the next time is a function of SOC in the previous time
    # with circular wrapping for the first and last timepoints within a timeseries
    println(" - State of charge constraints")
    if model_settings["single_storage_injection"] == true
        @constraint(mod, cStateOfCharge_SingleInjection[e ∈ E, s ∈ S, t ∈ T],
                        vSOC[e, s, t] == vSOC[e, s, T[t.prev_timepoint_id]] -
                                        t.duration_hr*(vDISCHA[e, s, t]) )
        # Power generation by bus
        @expression(mod, eNetDischargeAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    sum(vDISCHA[e, s, t] for e ∈ E_AT_BUS[n.id]) )

    else
        @constraint(mod, cStateOfCharge[e ∈ E, s ∈ S, t ∈ T],
                        vSOC[e, s, t] == vSOC[e, s, T[t.prev_timepoint_id]] +
                                        t.duration_hr*(vCHARGE[e, s, t]*e.efficiency_charge 
                                                        - vDISCHA[e, s, t]*1/e.efficiency_discharge) )
        # Power generation by bus
        @expression(mod, eNetDischargeAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    - sum(vCHARGE[e, s, t] for e ∈ E_AT_BUS[n.id]) 
                    + sum(vDISCHA[e, s, t] for e ∈ E_AT_BUS[n.id]) )
    end
    
    println(" - Storage cost per timepoint expressions")
    if model_settings["single_storage_injection"] == true
        @expression(mod, eStorCostPerTp[t ∈ T], 0.0)  # No variable cost in single injection model
    else
        @expression(mod, eStorCostPerTp[t ∈ T],
                     1/length(S)*(sum(s.probability * e.cost_variable_USDperMWh * (vCHARGE[e, s, t] + vDISCHA[e, s, t]) for e ∈ E, s ∈ S) ) )
    end

    eCostPerTp =  @views mod[:eCostPerTp]
    unregister(mod, :eCostPerTp)
    @expression(mod, eCostPerTp[t ∈ T], eCostPerTp[t] + eStorCostPerTp[t])
                 
    # Storage cost per period
    println(" - Storage cost per period expression")
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
    df1 = to_df(mod[:vDISCHA], [:storage, :scenario, :timepoint, :discharge_MW]; struct_fields = [:name, :name, :name])
    
    if haskey(mod, :vCHARGE)
        df2 = to_df(mod[:vCHARGE], [:storage, :scenario, :timepoint, :charge_MW]; struct_fields = [:name, :name, :name])
        df1 = outerjoin(df1, df2, on=[:storage, :scenario, :timepoint])
    end
   
    df3 = to_df(mod[:vSOC], [:storage, :scenario, :timepoint, :state_of_charge_MWh]; struct_fields = [:name, :name, :name])

    df_mix1 = outerjoin(df1, df3, on=[:storage, :scenario, :timepoint])

    CSV.write(joinpath(outputs_dir, "storage_dispatch.csv"), df_mix1)
    println("   - storage_dispatch.csv printed.")
    
    # Print vGENV variable solution
    df1 = to_df(mod[:vPCAP], [:storage, :scenario, :built_power_capacity_MW]; struct_fields = [:name, :name])

    # Print vCAPV variable solution
    df2 = to_df(mod[:vECAP], [:storage, :scenario, :built_energy_capacity_MWh]; struct_fields = [:name, :name])

    df_mix1 = outerjoin(df1, df2, on=[:storage, :scenario])
    CSV.write(joinpath(outputs_dir, "storage_capacity.csv"), df_mix1)
    println("   - storage_capacity.csv printed.")

    # Print cost expressions
    filename = "storage_costs_summary.csv"
    costs =  DataFrame(component  = ["CostPerTimepoint_USD", "CostPerPeriod_USD", "TotalCost_USD"], 
                            cost  = [   value(sum(mod[:eStorCostPerTp][t] * t.weight for t in sys.T)), 
                                        value(mod[:eStorCostPerPeriod]), 
                                        value(mod[:eStorTotalCost])]) 
    CSV.write(joinpath(outputs_dir, filename), costs)
    println("   - $filename printed.")

end


end # module EnergyStorage