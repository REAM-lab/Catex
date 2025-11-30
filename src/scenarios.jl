module Scenarios   

    using NamedArrays

    """
    Scenario represents a scenario in the stochastic capacity expansion problem.
    
    # Fields:
    - sc_id: ID of the scenario
    - prob: probability of the scenario
    """
    struct Scenario
        sc_id:: String
        prob:: Float64
    end

    function load_scenarios(inputs_dir:: String)

    # Get a list of Scenario structures
    scen = to_Structs(Scenario, inputs_dir, "scenarios.csv")

    # Get a list of the scenario IDs
    S = getfield.(scen, :sc_id)

    # Transform scen into NamedArray, so we can access scenarios by their IDs
    scen = NamedArray(scen, (S))

    return scen

    end

end # module Scenarios