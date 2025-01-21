module DoloYAML

    import Dolo:ndims

    import Dolo
    import Dolo: getprecision
    import Dolo: ⫫
    import Dolo: get_domain
    import Dolo: MarkovChain, MvNormal, VAR1, UNormal
    import Dolo: CartesianSpace, ProductSpace, GridSpace
    import Dolo: transition, arbitrage,complementarities
    import Dolo: AbstractModel, AModel
    using Dolo: YModel 

    const Normal = MvNormal
    
    # using Dolo: ModelCalibration
    using Dolang: sanitize
    using Dolang: Tree, Token
    import Dolang: solve_triangular_system
    import Dolang
    import Dolang: yaml_node_from_string, yaml_node_from_file
    using Dolang: SymExpr
    import Dolang: LTree
    import Dolang

    import YAML


    using DataStructures: OrderedDict


    using Printf

    using StaticArrays

    import Dolang: SymExpr, list_syms

    # model import utils

    using LinearAlgebra
    
    import YAML; using YAML: load_file, load

    import HTTP


    # Dolang
    using Dolang
    using Dolang: _to_expr, inf_to_Inf, solution_order, solve_triangular_system, _get_oorders
    import Dolang: Language, add_language_elements!, ToGreek

    using MacroTools  # used for eval_with

    using StringDistances


    import Base.size
    import Base.eltype
    import Base.*
    using Format
    using LinearAlgebra


    Expression = Union{Expr, Symbol, Float64, Int64}
    

    # recursively make all keys at any layer of nesting a symbol
    # included here instead of util.jl so we can call it on RECIPES below
    _symbol_dict(x) = x
    _symbol_dict(d::AbstractDict) =
        Dict{Symbol,Any}([(Symbol(k), _symbol_dict(v)) for (k, v) in d])

    const src_path = dirname(@__FILE__)
    const pkg_path = dirname(src_path)
    const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))

    using StaticArrays

    minilang = Language(Dict())
    add_language_elements!(minilang, Dict(
        "!Normal"=>UNormal,
        "!UNormal"=>UNormal,
        "!MvNormal"=>MvNormal,
        "!MarkovChain"=>MarkovChain,
        "!VAR1"=>VAR1,
        # "!Product"=>Product,
        # "!PoissonProcess"=>PoissonProcess,
        # "!ConstantProcess"=>ConstantProcess,
        # "!DeathProcess"=>DeathProcess,
        # "!AgingProcess"=>AgingProcess,
        # "!Mixture"=>Mixture,
        # "!Bernouilli"=>Bernouilli,
    ))
    # ))

    # include("linter.jl")
    include("calibration.jl")
    include("minilang.jl")
    include("model.jl")
    include("domodel.jl")
    include("misc.jl")

end # module DoloYAML
