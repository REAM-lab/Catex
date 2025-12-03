module EnergyStorages

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP

# Use internal modules
using ..Utils

# Export variables and functions
export EnergyStorage, load_data

"""
Energy storage represents an energy storage system in the power system.
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
struct EnergyStorage 
    es_id:: String
    es_tech:: String
    bus_id:: String
    invest_cost:: Float64
    exist_power_cap:: Float64
    exist_energy_cap:: Float64
    var_om_cost:: Float64
    charge_effic:: Float64
    discha_effic:: Float64
    duration:: Float64
end

function stochastic_capex_model!(mod:: Model, sys, pol)

    N = @views sys.N
    T = @views sys.T
    S = @views sys.S

    Nids = getfield.(N, :bus_id)

    E_AT_BUS = [filter(e -> e.bus_id == n, E) for n in Nids]

    # Define generation variables
    @variables(mod, begin
            vDISCHA[E, S, T] ≥ 0
            vCHARGE[E, S, T] ≥ 0       
            vSOC[E, S, T] ≥ 0
    end)

    # Define constraints for energy storage systems
    @constraint(mod, cMaxCharge[e ∈ E, s ∈ S, t ∈ T], 
                        vCHARGE[e, s, t] ≤ e.max)
    
    # Energy capacity must be less ×4 the power capacity
    #@constraint(mod, cMaxStateOfCharge[e ∈ E, s ∈ S, t ∈ T], 
    #                    vSOC[e, s, t] ≤ e.duration * e.exist_power_cap)
    
    # SOC in the next time is a function of SOC in the pervious time
    # with circular wrapping for the first and last t ∈ P_i
    
    prev = [idx == 1 ? T[end] : T[idx - 1] for (idx,_) in enumerate(T)]

    @constraint(mod, cStateOfCharge[e ∈ E, s ∈ S, t ∈ T],
                        vSOC[e, s, t] == vSOC[e, s, prev[t.idx]] + vCHARGE[e, s, t]*e.charge_effic - vDISCHA[e, s, t]*1/e.discha_effic)
end
end # module EnergyStorage