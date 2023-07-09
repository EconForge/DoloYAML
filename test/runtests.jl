module DoloTests

using Test

    using DoloYAML, DataStructures

    tests = length(ARGS) > 0 ? ARGS : [
                                   "test_import"
    ]

    for t in tests
        include("$(t).jl")
    end

end  # module
~
~
