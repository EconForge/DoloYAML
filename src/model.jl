# defind language understood in yaml files
language = minilang

mutable struct ModelSource{ID, Dom, P} <: AbstractModel{Dom}

    data::YAML.MappingNode
    symbols::Dict{Symbol, Vector{Symbol}}
    calibration::ModelCalibration
    parameters::P
    exogenous
    domain::Dom
    factories
    definitions
    filename::String

end

struct DummyModel{ID, Dom, P} <: AbstractModel{Dom}
    parameters::P
end

id(model::ModelSource{ID, Dom}) where ID where  Dom = ID

function check_exogenous_domain_type(data)
    if !("exogenous" in keys(data))
        return EmptyProcess
    end
    exo = data["exogenous"]
    exo_types = [exo[k].tag for k in keys(exo)]

    iid_types = ("!UNormal", "!ULogNormal", "!Normal", "!MvNormal", "!ConstantProcess")
    mc_types = ("!MarkovChain", "!ConstantProcess")
    acorr_types = ("!AR1", "!VAR1", "!ConstantProcess")

    if exo_types ⊆ iid_types
        return  EmptyDomain
    elseif exo_types ⊆ mc_types
        return DiscreteDomain
    else
        return CartesianDomain
    end
end


function ModelSource(url::AbstractString; print_code=false)
    
    fname = basename(url)

    # typ = check_exogenous_domain_type(data)
    # it looks like it would be cool to use the super constructor ;-)
    if match(r"(http|https):.*", url) != nothing
        res = HTTP.request("GET", url)
        txt = String(res.body)
    else
        txt = read(open(url), String)
    end
    data = Dolang.yaml_node_from_string(txt)
    # domain = det_domain(data, states, calib)
    # fspace = ProductDomain(typ, typeof())
    return ModelSource(data, fname)

end


function ModelSource(data::YAML.Node, filename)

    idt = gensym()

    calibration = get_calibration(data)
    symbols = get_symbols(data)
    exogenous = get_exogenous(data, symbols[:exogenous], calibration.flat)

    endo_domain = get_domain(data, symbols[:states], calibration.flat)
    exo_domain = get_domain(exogenous)
    
    if exogenous isa MvNormal
        domain = endo_domain
    else
        domain = ProductSpace(exo_domain, endo_domain)
    end

    dom_typ = typeof(domain)

    p = SVector(calibration[:parameters]...)

    model = ModelSource{idt, dom_typ, typeof(p)}(
        data,
        symbols,
        calibration,
        p,
        exogenous,
        domain,  ### how is this computed?
        nothing,
        nothing,
        filename
    )

    create_factories!(model)

    return model

end

function create_factories!(model)

    factories = get_factories(model)
    idt = id(model)
        
    model_type = Dolo.YModel{idt}

    for (k,fact) in factories
        # code = Dolang.gen_generated_gufun(fact; dispatch=typeof(model))
        # code = Dolang.gen_generated_gufun(fact; dispatch=model_type)
        code = Dolang.gen_kernel2(fact, 0; dispatch=model_type)
        # print_code && println("equation '", eq_type, "'", code)
        println(code)
        Core.eval(DoloYAML, code)
    end

    # Create definitions
    defs = get_definitions(model; stringify=true)
    model.definitions = defs
    definitions = OrderedDict{Symbol, SymExpr}( [(stringify(k),v) for (k,v) in defs])
    vars = cat(model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls]; dims=1)
    args = OrderedDict(
        :past => [Dolang.stringify(v,-1) for v in vars],
        :present => [Dolang.stringify(v,0) for v in vars],
        :future => [Dolang.stringify(v,1) for v in vars],
        :params => [Dolang.stringify(e) for e in model.symbols[:parameters]]
    )
    ff = Dolang.FunctionFactory(definitions, args, OrderedDict{Symbol, SymExpr}(), :definitions)
    # println("Compiling Definitions")
    # @show ff
    # code = Dolang.gen_kernel2(ff, 0; dispatch=model_type, funname=:evaluate_definitions)
    # Core.eval(DoloYAML,code)


    # # add compatibility
    # TODO: replace
    symbols = get_symbols(model)
    if (model.exogenous isa DoloYAML.VAR1) | (model.exogenous isa DoloYAML.MarkovChain)
        symbols[:states] = cat(symbols[:exogenous], symbols[:states]; dims=1)
    end

    n_p = length(symbols[:parameters])
    n_x = length(symbols[:controls])
    n_s = length(symbols[:states])
    n_e = length(symbols[:exogenous])

    code = quote

        function transition(model::$model_type, s::SVector{$n_s}, x::SVector{$n_x}, E::SVector{$n_e})
            # m_,s_ = get_ms(model,s)
            p = SVector((model.calibration[k] for k=1:$n_p)...)
            transition(model, s, x, E, p)
        end

        function arbitrage(model::$model_type, s::SVector{$n_s}, x::SVector{$n_x}, S::SVector{$n_s}, X::SVector{$n_x})
            # m_,s_ = get_ms(model,s)
            # M_,S_ = get_ms(model,S)
            p = SVector((model.calibration[k] for k=1:$n_p)...)
            arbitrage(model, s, x, S, X, p)
        end

        function controls_lb(model::$model_type, s::SVector{$n_s})
            p = SVector((model.calibration[k] for k=1:$n_p)...)
            controls_lb(model, s, p)
        end

        function controls_ub(model::$model_type, s::SVector{$n_s})
            p = SVector((model.calibration[k] for k=1:$n_p)...)
            controls_ub(model, s, p)
        end

        function complementarities(model::$model_type, s_::Dolo.QP, x::SVector{$n_x}, Fv::SVector{$n_x})
            s_ = s_.val
            # p = SVector((model.calibration[k] for k=1:$n_p)...)
            a = controls_lb(model, s_)
            b = controls_ub(model, s_)
            r = -( -(Fv .⫫ (x-a)) .⫫ (b-x) )
            return r
        end
    end

    println(code)
    Core.eval(DoloYAML, code)

    # end
    

    # if model.exogenous isa MvNormal
    #     println("Specializing for MvNormal")
    #     n_p = length(model.symbols[:parameters])
    #     n_x = length(model.symbols[:controls])
    #     n_s = length(model.symbols[:states])
    #     n_e = length(model.symbols[:exogenous])
    #     code = quote
            
    #         function transition(model::$model_type, s::SVector{$n_s}, x::SVector{$n_x}, E::SVector{$n_e})
    #             p = SVector((model.calibration[k] for k=1:$n_p)...)
    #             # here we should be allowed to ignore te first appearance of E
    #             T = getprecision(model)
    #             __E = zero(SVector{$n_e, T})
    #             transition(model, s, x, E, p)
    #         end

    #         function arbitrage(model::$model_type, s::SVector{$n_s}, x::SVector{$n_x}, S::SVector{$n_s}, X::SVector{$n_x})
    #             # m_,s_ = get_ms(model,s)
    #             # M_,S_ = get_ms(model,S)
    #             p = SVector((model.calibration[k] for k=1:$n_p)...)
    #             # here we should be allowed to ignore `m_` and `M_`
    #             T = getprecision(model)
    #             __E = zero(SVector{$n_e, T})

    #             arbitrage(model, s, x, S, X, p)
    #         end

    #         function controls_lb(model::$model_type, s::SVector{$n_s})
    #             p = SVector((model.calibration[k] for k=1:$n_p)...)
    #             controls_lb(model, s, p)
    #         end

    #         function controls_ub(model::$model_type, s::SVector{$n_s})
    #             p = SVector((model.calibration[k] for k=1:$n_p)...)
    #             controls_ub(model, s, p)
    #         end

    #         function complementarities(model::$model_type, s_::Dolo.QP, x::SVector{$n_x}, Fv::SVector{$n_x})
    #             s_ = s_.val
    #             # p = SVector((model.calibration[k] for k=1:$n_p)...)
    #             a = controls_lb(model, s_)
    #             b = controls_ub(model, s_)
    #             r = -( -(Fv .⫫ (x-a)) .⫫ (b-x) )
    #             return r
    #         end
    #     end
    #     println(code)
    #     Core.eval(DoloYAML, code)

    # end

    factories[:definitions] = ff

    model.factories = factories

end



yaml_import_old(filename::AbstractString) = ModelSource(filename)

function Base.show(io::IO, model::ModelSource)
    print(io, "Model")
end

function get_symbols(model::ModelSource)
    return get_symbols(model.data)
end

function get_symbols(data)
    syms = data[:symbols]
    symbols = OrderedDict{Symbol,Vector{Symbol}}()
    for k in keys(syms)
        symbols[Symbol(k)] = [Symbol(e.value) for e in syms[k]]
    end
    return symbols
end

function get_variables(model::ModelSource)
    return get_variables(model.data)
end

function get_variables(data)
    symbols = get_symbols(data)
    vars = cat(values(symbols)...; dims=1)
    dynvars = setdiff(vars, symbols[:parameters] )
    dynvars = union( dynvars, get_defined_variables(data) )
end

function get_defined_variables(model::ModelSource)
    return get_defined_variables(model.data)
end

function get_defined_variables(data)
    if !("definitions" in keys(data))
        return Symbol[]
    else
        # TODO: this is awfully redundant. This should be done "after" definitions are parsed
        try
            # OBSOLETE
            return [Symbol(e) for e in keys(data["definitions"])]
        catch
            defeqs = Dolang.parse_assignment_block(data["definitions"].value).children
            return [Symbol(eq.children[1].children[1].children[1].value) for eq in defeqs]

        end
    end
end

function get_name(model::ModelSource)
    if "name" in keys(model.data)
        name = model.data[:name].value
    else
        name = "anonymous"
    end
    return name
end

function get_exogenous(model::AModel)
    data = model.data
    exosyms = model.symbols[:exogenous]
    fcalib = model.calibration.flat
    return get_exogenous(data, exosyms, fcalib)
end

function get_exogenous(data, exosyms, fcalib)

    calibration = fcalib

    if !("exogenous" in keys(data))
        return nothing
    end

    exo_dict = data[:exogenous]
    cond = !(exo_dict.tag=="tag:yaml.org,2002:map")

    syms = cat( [[Symbol(strip(e)) for e in split(k, ",")] for k in keys(exo_dict)]..., dims=1)
    
    expected = exosyms
    if (syms != expected)
        msg = string("In 'exogenous' section, shocks must be declared in the same order as shock symbol. Found: ", syms, ". Expected: ", expected, ".")
        throw(ErrorException(msg))
    end
    # processes = []
    # for k in keys(exo_dict)
    #     v = exo_dict[k]
    #     p = Dolang.eval_node(v, calibration, minilang, ToGreek())
    #     push!(processes, p)
    # end
    processes = Dolang.eval_node_with_keys(exo_dict, calibration, minilang, ToGreek())
    
    if length(processes) > 1
        return ProductProcess(Tuple(processes))
    else 
        return processes[1]
    end

end


function get_calibration(model::ModelSource; kwargs...)

    return get_calibration(model.data; kwargs...)
end

function get_calibration(data; kwargs...)
    symbols = get_symbols(data)
    calib = data[:calibration]
    # prep calib: parse to Expr, Symbol, or Number
    _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
    # add calibration for definitions
    for (k,v) in kwargs
        calib[k].value = string(v)
    end
    for k in keys(calib)
        x = Dolang.parse_equation(calib[k].value)
        e = Dolang.convert(Expr, x, stringify=false)
        _calib[Symbol(k)] = e
    end
    # so far _calib is a symbolic calibration
    calibration = solve_triangular_system(_calib)
    return ModelCalibration(calibration, symbols)
end

function set_calibration!(model::ModelSource, key::Symbol, value)
    # TODO: set proper type for ScalarNode
    # this will fail is parameter wasn't defined before
    data = model.data
    model.data[:calibration][key].value = string(value)
    calibration = get_calibration(model)
    symbols = get_symbols(data)

    exogenous = get_exogenous(data, symbols[:exogenous], calibration.flat)

    endo_domain = get_domain(data, symbols[:states], calibration.flat)
    exo_domain = get_domain(exogenous)

    domain = ProductDomain(exo_domain, endo_domain)

    model.calibration = calibration
    model.exogenous = exogenous
    model.domain = domain

end

function set_calibration!(model::ModelSource; values...)
    calib = model.data[:calibration]
    data = model.data
    for (k,v) in values
        calib[k].value = string(v)
    end
    calibration = get_calibration(model)
    model.calibration = calibration
    
    exogenous = get_exogenous(model)
    symbols = get_symbols(data)
    endo_domain = get_domain(data, symbols[:states], calibration.flat)
    exo_domain = get_domain(model.exogenous)
    domain = ProductDomain(exo_domain, endo_domain)
    model.exogenous = exogenous
    model.domain = domain;
end



function get_equation_block(model, eqname; stringify=true)
    variables = get_variables(model)
    variables_str = String[string(e) for e in variables]
    
    yml_node = model.data["equations"][eqname]

    if "tag:yaml.org,2002:str" == yml_node.tag
        block = Dolang.parse_equation_block(yml_node.value; variables=variables_str)
        equations = block.children
    else
        equations = [Dolang.parse_equation(line.value; variables=variables_str) for line in yml_node]
    end

    eqs = Union{Expr, Number, Symbol}[]
    eqs_lb = Union{Expr, Number, Symbol}[]
    eqs_ub = Union{Expr, Number, Symbol}[]    

    for eq in equations
        if eq.data == "double_complementarity"
            inequality = eq.children[2]
            # @assert inegality.data == "double_inequality"
            eq_lb = Dolang.convert(Expr, inequality.children[1]; stringify=stringify)
            eq_cc = Dolang.convert(Expr, inequality.children[2]; stringify=stringify)
            # TODO: check complementarities order
            eq_ub = Dolang.convert(Expr, inequality.children[3]; stringify=stringify)
            eq = eq.children[1]
            push!(eqs_lb, eq_lb)
            push!(eqs_ub, eq_ub)
        else
            push!(eqs_lb, :(-Inf) )
            push!(eqs_ub, :(Inf) )
        end
        if eq.data == "equality"
            eq = Tree("sub", [eq.children[2], eq.children[1]])          
        end
        # eq = Dolang.stringify(eq)
        push!(eqs, Dolang.convert(Expr, eq; stringify=stringify))
    end
    return (eqs, eqs_lb, eqs_ub)
end

function get_assignment_block(model, eqname; stringify=true)
    variables = get_variables(model)
    variables_str = String[string(e) for e in variables]
    
    yml_node = model.data["equations"][eqname]

    if "tag:yaml.org,2002:str" == yml_node.tag
        block = Dolang.parse_assignment_block(yml_node.value; variables=variables_str)
        equations = block.children
    else
        equations = [Dolang.parse_assignment(line.value; variables=variables_str) for line in yml_node]
    end


    # eqs = OrderedDict{Tuple{Symbol, Int64}, SymExpr}()
    eqs = OrderedDict()
    for eq in equations

        lhs = eq.children[1]
        rhs = eq.children[2]
        dest_s = Symbol(lhs.children[1].children[1].value)
        dest_d = parse(Int, lhs.children[2].children[1].value)
        
        dest = Dolang.stringify(dest_s, dest_d)
        
        # TODO : remove type instability
        if stringify
            dest = Dolang.stringify(dest_s, dest_d)
        else
            dest = (dest_s, dest_d)
        end
        
        expr = Dolang.convert(Expr, rhs; stringify=stringify)

        eqs[dest] = expr

    end
    return eqs
end


function get_domain(model::AModel)::AbstractDomain
    states = model.symbols[:states]
    fcalib = model.calibration.flat
    data = model.data
    return get_domain(data, states, fcalib)
end

function get_domain(data, states, calib)

    if !("domain" in keys(data))
        return EmptyDomain(states)
    end

    domain = data["domain"]
    
    kk = tuple([Symbol(k) for k in keys(domain)]...)
    ss = tuple(states...)
    if kk != ss
        error("The declaration order in the domain section ($kk) must match the symbols orders ($ss)")
    end
    vals = [Dolang.eval_node(domain[k], calib) for k in kk]
    return CartesianSpace(;(k=>tuple(v...) for (k,v) in zip(states,vals))...)

end

function get_definitions(model::ModelSource; tshift=0, stringify=false) # ::OrderedDict{Tuple{Symbol,Int}}
    
    # return OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()

    # parse defs so values are Expr
    if !("definitions" in keys(model.data))
        return OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()
    end

    if model.data["definitions"].tag == "tag:yaml.org,2002:map"

        # LEGACY
        defs = Dict()

        dynvars = string.( cat(get_variables(model), keys(defs)...; dims=1) )

        _defs =  OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()  # {Symbol,Int}()

        for (kkk,v) in defs
            kk = Dolang.parse_equation(kkk)
            if kk.data == "symbol"
                k = kk.children[1].value
            elseif kk.data=="variable"
                k = kk.children[1].children[1].value
            else
                error("Invalid key.")
            end
            vv = Dolang.parse_equation(v.value; variables=(dynvars))::LTree
            if tshift != 0
                vv = Dolang.time_shift(vv, tshift)::LTree
            end
            _defs[(Symbol(k),tshift)] = Dolang.convert(Expr, vv; stringify=stringify)::SymExpr
        end

        return _defs

    else

        @assert model.data["definitions"].tag == "tag:yaml.org,2002:str"
        _defs =  OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()  # {Symbol,Int}()
        equations = Dolang.parse_assignment_block(model.data["definitions"].value).children
        # TODO : check that equations are correct
        for eq in equations
            dest, val = eq.children
            # TODO : create conveniences functions in dolang to avoid the follwowing
            name = dest.children[1].children[1].value
            # TODO : check that name is already defined at this stage
            t = parse(Int, dest.children[2].children[1].value)

            @assert t == 0

            vv = Dolang.time_shift(val, tshift)

            _defs[(Symbol(name),tshift)] = Dolang.convert(Expr, vv; stringify=stringify)::SymExpr
        end
        return _defs
    end

end


function get_options(model::ModelSource)
    cc = model.calibration.flat
    if "options" in keys(model.data)
        return Dolang.eval_node(model.data["options"], cc, language)
    else
        return Dict()
    end
end

import Dolang: FunctionFactory
import Dolang: stringify

function get_factories(model::ModelSource)

    facts = Dict()
    for eq_type in keys(model.data["equations"])
        if eq_type != "arbitrage"
            factories = (get_factory(model, eq_type), )
        else
            factories = get_factory(model, eq_type)
        end
        for f in factories
            facts[f.funname] = f
        end
    end

    return facts

end

function get_factory(model::ModelSource, eq_type::String)


    # TODO: replace
    symbols = get_symbols(model)
    if (model.exogenous isa DoloYAML.VAR1) | (model.exogenous isa DoloYAML.MarkovChain)
        symbols[:states] = cat(symbols[:exogenous], symbols[:states]; dims=1)
    end

    if eq_type == "arbitrage"

        defs_0 = get_definitions(model; stringify=true)
        defs_1 = get_definitions(model; tshift=1, stringify=true)
        definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in merge(defs_0, defs_1)])

        eqs, eq_lb, eq_ub = get_equation_block(model, eq_type)


        #TODO : fix crazy bug: it doesn't work without the trailing underscore !
        equations = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
            Symbol(string("out_", i, "_")) => eqs[i]
            for i =1:length(eqs)
        )
        arguments = OrderedDict(
            :s => [stringify(e,0) for e in symbols[:states]],
            :x => [stringify(e,0) for e in symbols[:controls]],
            :S => [stringify(e,1) for e in symbols[:states]],
            :X => [stringify(e,1) for e in symbols[:controls]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )

            
        
        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))

        equations_lb = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i, "_")) => eq_lb[i]
                for i =1:length(eqs)
            )
        equations_ub = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i,"_")) => eq_ub[i]
                for i =1:length(eqs)
            )
        arguments_bounds = OrderedDict(
            # :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :p => [stringify(e) for e in symbols[:parameters]]
        )

        # definitions: we should remove definitions depending on controls
        # definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in defs_0] )
        definitions = OrderedDict{Symbol, SymExpr}()
        ff_lb = FunctionFactory(equations_lb, arguments_bounds, definitions, :controls_lb)
        ff_ub = FunctionFactory(equations_ub, arguments_bounds, definitions, :controls_ub)

        return (ff, ff_lb, ff_ub)
           
    elseif eq_type=="transition"
        # defs_0  = get_definitions(model; stringify=true)
        defs_m1 = get_definitions(model; tshift=-1, stringify=true)

        definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

        equations = get_assignment_block(model, eq_type)

        arguments = OrderedDict(
            :s => [stringify(e,-1) for e in symbols[:states]],
            :x => [stringify(e,-1) for e in symbols[:controls]],
            :M => [stringify(e,0) for e in symbols[:exogenous]],
            :p => [stringify(e) for e in symbols[:parameters]]

            
        )
        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
        return ff
    else
        specs = RECIPES[:dtu][:specs]
        if !(Symbol(eq_type) in keys(specs))
            error("No spec for equation type '$eq_type'")
        end
        recipe = specs[Symbol(eq_type)]
        if !(:target in keys(recipe))
            error("Not implemented")
        end
        equations = get_assignment_block(model, eq_type)
        symbols = get_symbols(model)
        times = unique([e[2] for e in recipe[:eqs]])
        defs = [get_definitions(model; tshift=t, stringify=true) for t in times]
        definitions = OrderedDict{Symbol, SymExpr}(
            [ (Dolang.stringify(k), v) for (k,v) in merge(defs...)]
        )
        arguments = OrderedDict(
            Symbol(l[3]) => [stringify(e,l[2]) for e in symbols[Symbol(l[1])]]
            for l in recipe[:eqs] if !(l[1]=="parameters")
        )
        arguments[:p] = [stringify(e) for e in symbols[:parameters]]

        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
    end

end


function get_options(model::AModel; options=Dict())
    if "options" in keys(model.data)
        return model.data["options"]
        # opts = model.data["options"]
        # return rmerge(opts, options)
    else
        return Dict()
    end
end

# using Dolo: CartesianGrid

function get_discretization_options(model::AModel)

    domain = get_domain(model)
    d = length(domain.states)

    c = model.calibration.flat

    if ("options" in keys(model.data)) && ( "discretization" in keys(model.data["options"]) )
        gopt = model.data["options"]["discretization"]
        d = Dolang.eval_node(gopt, c)
        if !(:exo in keys(d))
            d[:exo] = Dict()
        end
        return d
    elseif ("options" in keys(model.data)) && ( "grid" in keys(model.data["options"]) )
        return Dict{Any,Any}(:exo => Dict{Any, Any}(), :endo => Dict{Any, Any}(:n=>DoloYAML.get_options(model)[:grid].orders))
    else
        return Dict(
            :endo=>Dict(),
            :exo=>Dict()
        )
    end


        
end


function discretize(model::ModelSource; kwargs...)

    opts = get_discretization_options(model; kwargs...)
    
    opts_endo = merge(opts[:endo], get(kwargs, :endo, Dict()) )
    opts_exo = merge(opts[:exo], get(kwargs, :exo, Dict()) )

    endo_domain = model.domain.endo
    grid_endo = discretize(endo_domain;  opts_endo...) 
    dprocess = discretize(model.exogenous;  opts_exo...) 

    grid_exo = dprocess.grid
    grid = ProductGrid(grid_exo, grid_endo)
    return grid, dprocess


    # dprocess = discretize(model.exogenous)
#     return ProductGrid(dprocess.grid, grid), dprocess
    # exo_grid = dprocess.
    # return ProductGrid(grid, exo_grid)
    # return dprocess.grid, grid, dprocess
end

