function DoloModel(filename::AbstractString)
    source = yaml_import_old(filename)
    DoloModel(source)
end

yaml_import(filename::AbstractString) = DoloModel(filename)


function DoloModel(source)

    m = NamedTuple( v=>source.calibration.flat[v] for v in source.symbols[:exogenous] )
    s = NamedTuple( v=>source.calibration.flat[v] for v in source.symbols[:states] ) 
    x = NamedTuple( v=>source.calibration.flat[v] for v in source.symbols[:controls] )
    p = NamedTuple( v=>source.calibration.flat[v] for v in source.symbols[:parameters] )

    calibration = merge(p, m, s, x)

    # if typeof(source.exogenous) <: MarkovChain

    #     (;states, min, max) = source.domain.endo

    #     aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))
    #     cspace = Dolo.CartesianSpace(; aargs... )
    #     exo_vars = tuple(source.symbols[:exogenous]...)
    #     states = get_domain(source.exogenous) × cspace
    #     exogenous = source.exogenous

    # elseif typeof(source.exogenous) <: MvNormal

    #     (;states, min, max) =  source.domain.endo
    #     aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))
    #     states = Dolo.CartesianSpace(; aargs... )
    #     exo_vars = tuple(source.symbols[:exogenous]...)
    #     exogenous = Dolo.MvNormal(exo_vars, source.exogenous.μ, source.exogenous.Σ)
    #     # exogenous = Dolo.MvNormal(source.exogenous.μ, source.exogenous.Σ) # TODO

    # elseif typeof(source.exogenous) <: VAR1
        
    #     exo_vars = tuple(source.symbols[:exogenous]...)
    #     # TODO: assert R diag(r)
    #     # ρ = source.exogenous.ρ
    #     Σ = source.exogenous.Σ
    #     p = size(Σ,1)


    #     (;states, min, max) = source.domain.endo

    #     aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))
    #     cspace = Dolo.CartesianSpace(aargs... )

    #     exo_domain = get_domain(source.exogenous)

    #     min = tuple(exo_domain.min..., cspace.min...)
    #     max = tuple(exo_domain.max..., cspace.max...)
    #     names = tuple(Dolo.variables(exo_domain)..., Dolo.variables(cspace)...)
        
    #     states = Dolo.CartesianSpace(; ( n=>(a,b) for (n,a,b) in zip(names, min, max))... )
    #     exogenous = source.exogenous

    # end

    states = source.domain
    exogenous = source.exogenous

    controls = Dolo.CartesianSpace(
        (k=>(-Inf, Inf) for k in source.symbols[:controls])...
    )

    name = try
        Symbol(source.data[:name].value)
    catch e
        :anonymous
    end

    return Dolo.YModel(id(source), states, controls, exogenous, calibration, source)
    # return Dolo.YModel(name, states, controls, exogenous, calibration, source)

end


# function discretize(ym::Dolo.YModel{<:Dolo.MarkovChain,B,C,D,E}) where A where B where C where D where E<:AbstractModel

function recalibrate(ym::Dolo.YModel{A,B,C,D,E, F}; kwargs...) where A where B where C where D where E where F<:AbstractModel
    
    source = (ym.source)

    DoloYAML.set_calibration!(source; kwargs...)
    m = (; (v=>source.calibration.flat[v] for v in source.symbols[:exogenous])...) 
    s = (; (v=>source.calibration.flat[v] for v in source.symbols[:states])...) 
    x = (; (v=>source.calibration.flat[v] for v in source.symbols[:controls])...) 
    p = (; (v=>source.calibration.flat[v] for v in source.symbols[:parameters])...) 
    calibration = merge(m, s, x, p,)

    Dolo.YModel(Dolo.name(ym), ym.states, ym.controls, ym.exogenous, calibration, source)
    
end

function discretize(ym::Dolo.YModel{<:Dolo.MvNormal,B,C,D,E, F}; endo=nothing, exo=nothing) where B where C where D where E where F<:AbstractModel

    options = DoloYAML.get_options(ym.source)
    discopts = get(options, :discretization, Dict())
    if typeof(endo)<:Nothing
        endo_opts = get(discopts, :endo, Dict())
    else
        endo_opts = endo
    end
    if typeof(exo)<:Nothing
        exo_opts = get(discopts, :exo, Dict())
    else
        exo_opts = exo
    end

    dvar = discretize(ym.exogenous, exo_opts)
    grid = discretize(ym.states, endo_opts)

    return Dolo.DYModel(ym, grid, dvar)
    
end
# only for VAR and MC


function discretize(ym::Dolo.YModel{<:Dolo.MarkovChain,B,C,D,E, F}; endo=nothing, exo=nothing) where B where C where D where E where F<:AbstractModel

    options = DoloYAML.get_options(ym.source)
    discopts = get(options, :discretization, Dict())
    if typeof(endo)<:Nothing
        endo_opts = get(discopts, :endo, Dict())
    else
        endo_opts = endo
    end
    if typeof(exo)<:Nothing
        exo_opts = get(discopts, :exo, Dict())
    else
        exo_opts = exo
    end

    dvar = discretize(ym.exogenous, exo_opts)
    exo_grid = SGrid(dvar.Q)
    endo_grid = discretize(ym.states.spaces[2], endo_opts)
    grid = exo_grid × endo_grid
    return Dolo.DYModel(ym, grid, dvar)
    
end


function discretize(ym::Dolo.YModel{<:Dolo.VAR1,B,C,D,E, F}; endo=nothing, exo=nothing) where B where C where D where E where F<:AbstractModel

    options = get_options(ym.source)
    discopts = get(options, :discretization, Dict())
    if typeof(endo)<:Nothing
        endo_opts = get(discopts, :endo, Dict())
    else
        endo_opts = endo
    end
    if typeof(exo)<:Nothing
        exo_opts = get(discopts, :exo, Dict())
    else
        exo_opts = exo
    end

    dvar = discretize(ym.exogenous, exo_opts)
    exo_grid = SGrid(dvar.Q)
    
    d = Dolo.ndims(exo_grid)
    min = ym.states.min[d+1:end]
    max = ym.states.max[d+1:end]
    endo_grid = discretize(
        Dolo.CSpace(min, max),
        endo_opts
    )
    grid = exo_grid × endo_grid
    return Dolo.DYModel(ym, grid, dvar)
    
end

# only for VAR and MC

function transition(model::Dolo.YModel{MOD,B,C,D,E,sS}, s::NamedTuple, x::NamedTuple, M::NamedTuple) where B where C where D where E where sS<:AbstractModel where MOD<:Union{<:Dolo.MarkovChain,<:Dolo.VAR1}

    svars = (Dolo.variables(model.states))
    exo = (Dolo.variables(model.exogenous))
    n_e = length(exo)
    endo = svars[n_e+1:length(svars)]

    controls = Dolo.variables(model.controls)

    mm = SVector( (s[i] for i in exo)... )
    ss = SVector( (s[i] for i in endo)... )
    xx = SVector( (x[i] for i in controls)... )
    MM = SVector( (M[i] for i in exo)... )

    # p = SVector(model.source.calibration[:parameters]...)
    p = model.source.parameters
    S = transition(model.source, mm, ss, xx, MM, p)

    return NamedTuple{endo}(S)

end

# function transition(model::Dolo.YModel{MOD,B,C,D,E,sS}, s::SVector, x::SVector, M::SVector) where B where C where D where E where sS<:AbstractModel where MOD<:Union{<:Dolo.MarkovChain,<:Dolo.VAR1}

#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))

#     mm = SVector( (s[i] for i=1:d)... )
#     ss = SVector( (s[i] for i=(d+1):n)... )
    
#     # p = SVector(model.source.calibration[:parameters]...)
#     p = model.source.parameters

#     S = Dolo.transition(model.source, mm, ss, x, M, p)

#     return S

# end

function get_ms(model::Dolo.YModel, s::SVector)
    d = length(Dolo.variables(model.exogenous))
    n = length(Dolo.variables(model.states))
    mm = SVector( (s[i] for i=1:d)...)
    ss = SVector( (s[i] for i=(d+1):n)...)
    return (mm, ss)
end

# function get_ms(model::Dolo.YModel{<:Dolo.MvNormal,B,C,D,E,sS}, s::SVector) where B where C where D where E where sS<:AbstractModel
#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))
#     m = SVector( (NaN64 for i=1:d)... )
#     return (m,s)
# end


# function arbitrage(model::Dolo.YModel{A,B,C,D,E,sS}, s::SVector, x::SVector, S::SVector, X::SVector) where A where B where C where D where E where sS<:AbstractModel

#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))


    
#     mm,ss = get_ms(model, s)
#     # mm = SVector( (s[i] for i=1:d)...)
#     # ss = SVector( (s[i] for i=(d+1):n)...)
    
#     MM = SVector( (S[i] for i=1:d)...)
#     SS = SVector( (S[i] for i=(d+1):n)...)
    
#     p = model.source.parameters

#     res = arbitrage(model.source, mm, ss, x, MM, SS, X, p)

#     return SVector(res...)

# end


# function transition(model::Dolo.YModel{<:Dolo.MvNormal,B,C,D,E,sS}, s::SVector, x::SVector, M::SVector) where B where C where D where E where sS<:AbstractModel

#     p = model.source.parameters
#     M_ = M # this argument is ignored
#     S = transition(model.source, M_, s, x, M, p)

#     return S

# end


# function arbitrage(model::Dolo.YModel{<:Dolo.MvNormal,B,C,D,E,sS}, s::SVector, x::SVector, S::SVector, X::SVector)  where B where C where D where E where sS<:AbstractModel

#     # TODO: currently we don't allow for iid shocks in Euler equations
#     # TODO: warning

#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))
#     m0 = zero(SVector{d,Float64})

#     p = model.source.parameters

#     res = arbitrage(model.source, m0, s, x, m0, S, X, p)

#     return SVector(res...)

# end

# function bounds(model::Dolo.YModel{A, B, C, D, E, sS}, ss::Dolo.QP) where A where B where C where D where E where sS<:AbstractModel
    
#     s = ss.val

#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))

#     mm, ss_ = get_ms(model, ss.val)

#     p = model.source.parameters

#     lb = controls_lb(model.source, mm, ss_,  p)
#     ub = controls_ub(model.source, mm, ss_,  p)

#     return (lb, ub)


# end


# function reward(model::Dolo.YModel{A, B, C, D, E, sS}, ss::Dolo.QP, x::SVector) where A where B where C where D where E where sS<:AbstractModel
    

#     mm, ss_ = get_ms(model, ss.val)

#     p = model.source.parameters

#     rew = felicity(model.source, mm, ss_, x, p)
    

#     return rew[1]


# end



# function complementarities(model::Dolo.YModel{<:Dolo.MarkovChain, B, C, D, E, sS}, ss::Dolo.QP, x::SVector, Ef::SVector) where B where C where D where E where sS<:AbstractModel
    
#     s = ss.val

#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))
    
#     mm = SVector( (s[i] for i=1:d)... )
#     ss_ = SVector( (s[i] for i=(d+1):n)... )

#     p = model.source.parameters
#     lb = controls_lb(model.source, mm, ss_,  p)

#     r = (Ef .⫫ (x-lb)) 
#     # r = (-r) .⫫ (ub-x)
#     # -r

# end


# function complementarities(model::Dolo.YModel{<:Dolo.VAR1, B, C, D, E, sS}, ss::Dolo.QP, x::SVector, Ef::SVector) where B where C where D where E where sS<:AbstractModel
    
#     s = ss.val

#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))

#     mm = SVector( (s[i] for i=1:d)... )
#     ss_ = SVector( (s[i] for i=(d+1):n)... )

#     p = model.source.parameters
#     lb = controls_lb(model.source, mm, ss_,  p)
#     ub = controls_ub(model.source, mm, ss_,  p)

#     r = Ef .⫫ (x-lb)
#     # r = (-Ef) .⫫ (ub-x)
#     # -r
    
# end


# function complementarities(model::Dolo.YModel{<:Dolo.MvNormal, B, C, D, E, sS}, ss::Dolo.QP, x::SVector, Ef::SVector) where B where C where D where E where sS<:AbstractModel
    

#     d = length(Dolo.variables(model.exogenous))
#     n = length(Dolo.variables(model.states))
#     m0 = zero(SVector{d,Float64})
    
#     s = ss.val

#     p = model.source.parameters
#     lb = controls_lb(model.source, m0, s,  p)
#     ub = controls_ub(model.source, m0, s,  p)

#     # Ef
#     r = Ef .⫫ (x-lb)
#     r = (-Ef) .⫫ (ub-x)
#     -r

#     # r = Ef .⟂ (x-lb)
#     # r
#     # r = (-Ef) .⟂ (ub-x)
#     # -r

# end
