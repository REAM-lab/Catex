module EnergyStorage

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP, DataFrames, CSV

# Use internal modules
using ..Utils

# Export variables and functions
export EnergyStorageUnit, stochastic_capex_model!, toCSV_stochastic_capex

"""
Energy storage unit represents an energy storage system in the power system.
# Fields:
- es_id: ID of the storage system
- es_tech: technology type of the storage system, for example, "battery", "pumped_hydro". It could be any string.
- bus_id: ID of the bus where the storage system is connected to.
- invest_cost: investment cost per MW of power capacity (USD/MW)
- exist_power_cap: pre-existing power capacity of the storage system (MW)
- exist_energy_cap: pre-existing energy capacity of the storage system (MWh)
- var_om_cost: variable operation and maintenance cost (USD/MW)
- efficiency: round-trip efficiency of the storage system (between 0 and 1)
- duration: duration of the storage system at full power (hours)
"""
struct EnergyStorageUnit
    id:: Int64
    name:: String
    tech:: String
    bus_name:: String
    invest_cost:: Float64
    exist_power_cap:: Float64
    exist_energy_cap:: Float64
    power_cap_limit:: Float64
    var_om_cost:: Float64
    charge_effic:: Float64
    discha_effic:: Float64
    duration:: Float64
end

function stochastic_capex_model!(mod:: Model, sys, pol)

    N = @views sys.N
    T = @views sys.T
    S = @views sys.S
    E = @views sys.E

    E_AT_BUS = [filter(e -> e.bus_name == n.name, E) for n in N]

    # Define generation variables
    @variables(mod, begin
            vDISCHA[E, S, T] ≥ 0
            vCHARGE[E, S, T] ≥ 0       
            vSOC[E, S, T] ≥ 0
            vPCAP[E, S] ≥ 0  
            vECAP[E, S] ≥ 0  
    end)

    @constraint(mod, cMinEnerCapStor[e ∈ E, s ∈ S], 
                        vECAP[e, s] ≥ e.exist_energy_cap)
    
    @constraint(mod, cMinPowerCapStor[e ∈ E, s ∈ S], 
                        vPCAP[e, s] ≥ e.exist_power_cap)

    @constraint(mod, cMaxPowerCapStor[e ∈ E, s ∈ S], 
                        vPCAP[e, s] ≤ e.power_cap_limit)

    # Define constraints for energy storage systems
    @constraint(mod, cMaxCharge[e ∈ E, s ∈ S, t ∈ T], 
                        vCHARGE[e, s, t] ≤ vPCAP[e, s])
    
    @constraint(mod, cMaxDischa[e ∈ E, s ∈ S, t ∈ T], 
                        vDISCHA[e, s, t] ≤ vPCAP[e, s])
    
    E_fixdur = filter(e -> e.duration > 0, E)

    if !isempty(E_fixdur)
    @constraint(mod, cFixEnergyPowerRatio[e ∈ E, s ∈ S, t ∈ T], 
                        vECAP[e, s] ==  e.duration * vPCAP[e, s] )
    end

    # Energy capacity must be less ×4 the power capacity
    @constraint(mod, cMaxSOC[e ∈ E, s ∈ S, t ∈ T], 
                        vSOC[e, s, t] ≤ vECAP[e, s])

    # SOC in the next time is a function of SOC in the pervious time
    # with circular wrapping for the first and last t ∈ P_i
    
    prev = [id == 1 ? T[end] : T[id - 1] for (id,_) in enumerate(T)]

    @constraint(mod, cStateOfCharge[e ∈ E, s ∈ S, t ∈ T],
                        vSOC[e, s, t] == vSOC[e, s, prev[t.id]] 
                                        + vCHARGE[e, s, t]*e.charge_effic 
                                        - vDISCHA[e, s, t]*1/e.discha_effic)

    # Power generation by bus
    @expression(mod, eNetDischargeAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    - sum(vCHARGE[e, s, t] for e ∈ E_AT_BUS[n.id]) 
                    + sum(vDISCHA[e, s, t] for e ∈ E_AT_BUS[n.id]) )

    # 
    @expression(mod, eStorVariableCosts,
                     1/length(S)*(sum(s.prob * t.duration * e.var_om_cost * vCHARGE[e, s, t] for e ∈ E, s ∈ S, t ∈ T) ) ) 

    # Fixed costs 
	@expression(mod, eStorFixedCosts,
                    1/length(S)*sum( (s.prob * e.invest_cost * vPCAP[e, s]) for e ∈ E, s ∈ S ))

    # Total costs
    @expression(mod, eStorTotalCosts,
                        eStorVariableCosts + eStorFixedCosts)
    
end

function toCSV_stochastic_capex(sys, pol, mod:: Model, outputs_dir:: String)

    # Print vCHARGE AND vDISCHARGE 
    df1 = to_df(mod[:vCHARGE], [:es_name, :sc_name, :tp_name, :Charge_MW]; 
                                struct_fields=[:name, :name, :name])
    df2 = to_df(mod[:vDISCHA], [:es_name, :sc_name, :tp_name, :Discharge_MW];
                                struct_fields=[:name, :name, :name])
    
    df_mix1 = outerjoin(df1, df2, on=[:es_name, :sc_name, :tp_name])
    CSV.write(joinpath(outputs_dir, "storage_dispatch.csv"), df_mix1)

    to_df(mod[:vSOC], [:es_name, :sc_name, :tp_name, :SOC_MWh];
                    struct_fields=[:name, :name, :name], csv_dir = joinpath(outputs_dir, "storage_soc.csv"))
    
    # Print vGENV variable solution
    df1 = to_df(mod[:vPCAP], [:es_name, :sc_name, :PowerCap_MW]; 
                        struct_fields=[:name, :name])

    # Print vCAPV variable solution
    df2 = to_df(mod[:vECAP], [:es_name, :sc_name, :EnergyCap_MWh]; 
                        struct_fields=[:name, :name])

    df_mix1 = outerjoin(df1, df2, on=[:es_name, :sc_name])
    CSV.write(joinpath(outputs_dir, "storage_capacity.csv"), df_mix1)

    # Print cost expressions
    filename = "storage_costs_itemized.csv"
    costs =  DataFrame(component  = ["variable_costs", "fixed_costs", "total_costs"], 
                            cost  = [value(mod[:eStorVariableCosts]), value(mod[:eStorFixedCosts]), value(mod[:eStorTotalCosts])]) 
    CSV.write(joinpath(outputs_dir, filename), costs)
    println(" > $filename printed.")
    

end


end # module EnergyStorage