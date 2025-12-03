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
    id:: Int64
    name:: String
    tech:: String
    bus_name:: String
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

    E_AT_BUS = [filter(e -> e.bus_name == n.bus_name, E) for n in N]

    # Define generation variables
    @variables(mod, begin
            vDISCHA[E, S, T] ≥ 0
            vCHARGE[E, S, T] ≥ 0       
            vSOC[E, S, T] ≥ 0
            vPCAP[E, S] ≥ 0  
            vECAP[E, S] ≥ 0  
    end)

    # Define constraints for energy storage systems
    @constraint(mod, cMaxCharge[e ∈ E, s ∈ S, t ∈ T], 
                        vCHARGE[e, s, t] ≤ vPCAP[e, s])
    
    @constraint(mod, cMaxDischarge[e ∈ E, s ∈ S, t ∈ T], 
                        vDISCHA[e, s, t] ≤ vPCAP[e, s])
    
    E_fixdur = filter(e -> e.duration > 0, E)

    if !isempty(E_fixdur)
    @constraint(mod, cFixEnergyPowerRatio[e ∈ E, s ∈ S, t ∈ T], 
                        vECAP[e, s] ==  e.duration * vPCAP[e, s] )
    end

    # Minimum power capacity of storage
    @constraint(mod, cFixPowerCapStor[e ∈ E, s ∈ S], 
                    vPCAP[e, s] ≥ e.exist_power_cap)

    # Minimum energy capacity of storage
    @constraint(mod, cMaxEnergyCapStor[e ∈ E, s ∈ S], 
                    vECAP[e, s] ≥ g.exist_energy_cap)

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
    @expression(mod, eNetChargeAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    sum(vCHARGE[e, s, t] for e ∈ E_AT_BUS[n.idx]) 
                    - sum(vDISCHA[e, s, t] for e ∈ E_AT_BUS[n.idx]) )

    
end
end # module EnergyStorage