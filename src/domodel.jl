using StaticArrays
    

function DoloModel(filename::AbstractString)
    source = yaml_import(filename)
    DoloModel(source)
end

function DoloModel(source)

    m = (; (v=>source.calibration.flat[v] for v in source.symbols[:exogenous])...) 
    s = (; (v=>source.calibration.flat[v] for v in source.symbols[:states])...) 
    x = (; (v=>source.calibration.flat[v] for v in source.symbols[:controls])...) 
    p = (; (v=>source.calibration.flat[v] for v in source.symbols[:parameters])...) 

    calibration = merge(m, s, x, p,)

    if typeof(source.exogenous) <: DiscreteMarkovProcess
        p1,q1=size(source.exogenous.values)
        p2,q2=size(source.exogenous.transitions)

        Q = SMatrix{p1,q1}(source.exogenous.values)
        points = SVector( tuple((Q[i,:] for i=1:size(Q,1))...) )

        P = SMatrix{p2,q2}(source.exogenous.transitions)

        (;states, min, max) = source.domain.endo

        aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))
        cspace = NoLib.CartesianSpace(; aargs... )
        exo_vars = tuple(source.symbols[:exogenous]...)
        states = NoLib.GridSpace(exo_vars, points  ) × cspace
        exogenous = NoLib.MarkovChain(exo_vars, P, points)

    elseif typeof(source.exogenous) <: MvNormal

        (;states, min, max) = source.domain.endo
        aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))
        states = NoLib.CartesianSpace(; aargs... )
        exo_vars = tuple(source.symbols[:exogenous]...)
        exogenous = NoLib.MvNormal(exo_vars, source.exogenous.μ, source.exogenous.Σ)
        # exogenous = NoLib.MvNormal(source.exogenous.μ, source.exogenous.Σ) # TODO

    elseif typeof(source.exogenous) <: VAR1
        
        exo_vars = tuple(source.symbols[:exogenous]...)
        # TODO: assert R diag(r)
        ρ = source.exogenous.R[1,1]
        Σ = source.exogenous.Σ
        p = size(Σ,1)


        (;states, min, max) = source.domain.endo

        aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))
        cspace = NoLib.CartesianSpace(; aargs... )
        exo_vars = tuple(source.symbols[:exogenous]...)
        keys =  Dict(k=>(-Inf,Inf) for k in exo_vars )
        exo_domain = NoLib.CartesianSpace(; keys... )

        min = tuple(exo_domain.min..., cspace.min...)
        max = tuple(exo_domain.max..., cspace.max...)
        names = tuple(NoLib.variables(exo_domain)..., NoLib.variables(cspace)...)
        states = NoLib.CartesianSpace{length(names), names}(min, max)
        exogenous = NoLib.VAR1(exo_vars, ρ, SMatrix{p,p,Float64,p*p}(Σ))


    end

        controls = NoLib.CartesianSpace(;
            Dict(k=>[-Inf, Inf] for k in source.symbols[:controls])...
        )

    name = try
        Symbol(source.data[:name].value)
    catch e
        :anonymous
    end

    return NoLib.YModel(name, states, controls, exogenous, calibration, source)

end


# function discretize(ym::NoLib.YModel{<:NoLib.MarkovChain,B,C,D,E}) where A where B where C where D where E<:Model

function recalibrate(ym::NoLib.YModel{A,B,C,D,E, F}; kwargs...) where A where B where C where D where E where F<:Model
    
    source = (ym.source)

    DoloYAML.set_calibration!(source; kwargs...)
    m = (; (v=>source.calibration.flat[v] for v in source.symbols[:exogenous])...) 
    s = (; (v=>source.calibration.flat[v] for v in source.symbols[:states])...) 
    x = (; (v=>source.calibration.flat[v] for v in source.symbols[:controls])...) 
    p = (; (v=>source.calibration.flat[v] for v in source.symbols[:parameters])...) 
    calibration = merge(m, s, x, p,)

    NoLib.YModel(NoLib.name(ym), ym.states, ym.controls, ym.exogenous, calibration, source)
    
end

function discretize(ym::NoLib.YModel{<:NoLib.MvNormal,B,C,D,E, F}; endo=nothing, exo=nothing) where B where C where D where E where F<:Model

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

    return NoLib.DYModel(ym, grid, dvar)
    
end
# only for VAR and MC


function discretize(ym::NoLib.YModel{<:NoLib.MarkovChain,B,C,D,E, F}; endo=nothing, exo=nothing) where B where C where D where E where F<:Model

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
    return NoLib.DYModel(ym, grid, dvar)
    
end


function discretize(ym::NoLib.YModel{<:NoLib.VAR1,B,C,D,E, F}; endo=nothing, exo=nothing) where B where C where D where E where F<:Model

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
    
    d = NoLib.ndims(exo_grid)
    min = ym.states.min[d+1:end]
    max = ym.states.max[d+1:end]
    endo_grid = discretize(
        NoLib.CSpace(min, max),
        endo_opts
    )
    grid = exo_grid × endo_grid
    return NoLib.DYModel(ym, grid, dvar)
    
end

# only for VAR and MC


function transition(model::NoLib.YModel{MOD,B,C,D,E,sS}, s::NamedTuple, x::NamedTuple, M::NamedTuple) where B where C where D where E where sS<:Model where MOD<:Union{<:NoLib.MarkovChain,<:NoLib.VAR1}

    svars = (NoLib.variables(model.states))
    exo = (NoLib.variables(model.exogenous))
    endo = tuple( (v for v in svars if !(v in exo))...)
    controls = NoLib.variables(model.controls)

    mm = SVector( (s[i] for i in exo)... )
    ss = SVector( (s[i] for i in endo)... )
    xx = SVector( (x[i] for i in controls)... )
    MM = SVector( (M[i] for i in exo)... )

    # p = SVector(model.source.calibration[:parameters]...)
    p = model.source.parameters
    S = transition(model.source, mm, ss, xx, MM, p)

    return NamedTuple{endo}(S...)

end

# function transition(model::NoLib.YModel{MOD,B,C,D,E,sS}, s::SVector, x::SVector, M::SVector) where B where C where D where E where sS<:Model where MOD<:Union{<:NoLib.MarkovChain,<:NoLib.VAR1}

#     d = length(NoLib.variables(model.exogenous))
#     n = length(NoLib.variables(model.states))

#     mm = SVector( (s[i] for i=1:d)... )
#     ss = SVector( (s[i] for i=(d+1):n)... )
    
#     # p = SVector(model.source.calibration[:parameters]...)
#     p = model.source.parameters

#     S = Dolo.transition(model.source, mm, ss, x, M, p)

#     return S

# end

function get_ms(model::NoLib.YModel{A,B,C,D,E,sS}, s::SVector) where A where B where C where D where E where sS<:Model
    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))
    m = SVector( (s[i] for i=1:d)... )
    s = SVector( (s[i] for i=(d+1):n)... )
    return (m,s)
end

function get_ms(model::NoLib.YModel{<:NoLib.MvNormal,B,C,D,E,sS}, s::SVector) where B where C where D where E where sS<:Model
    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))
    m = SVector( (NaN64 for i=1:d)... )
    return (m,s)
end


function arbitrage(model::NoLib.YModel{A,B,C,D,E,sS}, s::SVector, x::SVector, S::SVector, X::SVector) where A where B where C where D where E where sS<:Model

    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))
    mm,ss = get_ms(model, s)
    
    MM = SVector( (S[i] for i=1:d)...)
    SS = SVector( (S[i] for i=(d+1):n)...)
    
    p = model.source.parameters

    res = arbitrage(model.source, mm, ss, x, MM, SS, X, p)

    return SVector(res...)

end


function transition(model::NoLib.YModel{<:NoLib.MvNormal,B,C,D,E,sS}, s::SVector, x::SVector, M::SVector) where B where C where D where E where sS<:Model

    p = model.source.parameters
    M_ = M # this argument is ignored
    S = transition(model.source, M_, s, x, M, p)

    return S

end


function arbitrage(model::NoLib.YModel{<:NoLib.MvNormal,B,C,D,E,sS}, s::SVector, x::SVector, S::SVector, X::SVector)  where B where C where D where E where sS<:Model

    # TODO: currently we don't allow for iid shocks in Euler equations
    # TODO: warning

    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))
    m0 = zero(SVector{d,Float64})

    p = model.source.parameters

    res = arbitrage(model.source, m0, s, x, m0, S, X, p)

    return SVector(res...)

end

function bounds(model::NoLib.YModel{A, B, C, D, E, sS}, ss::NoLib.QP) where A where B where C where D where E where sS<:Model
    
    s = ss.val

    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))

    mm, ss_ = get_ms(model, ss.val)

    p = model.source.parameters

    lb = controls_lb(model.source, mm, ss_,  p)
    ub = controls_ub(model.source, mm, ss_,  p)

    return (lb, ub)


end


function reward(model::NoLib.YModel{A, B, C, D, E, sS}, ss::NoLib.QP, x::SVector) where A where B where C where D where E where sS<:Model
    

    mm, ss_ = get_ms(model, ss.val)

    p = model.source.parameters

    rew = felicity(model.source, mm, ss_, x, p)
    

    return rew[1]


end



function complementarities(model::NoLib.YModel{<:NoLib.MarkovChain, B, C, D, E, sS}, ss::NoLib.QP, x::SVector, Ef::SVector) where B where C where D where E where sS<:Model
    
    s = ss.val

    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))
    
    mm = SVector( (s[i] for i=1:d)... )
    ss_ = SVector( (s[i] for i=(d+1):n)... )

    p = model.source.parameters
    lb = controls_lb(model.source, mm, ss_,  p)

    r = (Ef .⫫ (x-lb)) 
    # r = (-r) .⫫ (ub-x)
    # -r

end


function complementarities(model::NoLib.YModel{<:NoLib.VAR1, B, C, D, E, sS}, ss::NoLib.QP, x::SVector, Ef::SVector) where B where C where D where E where sS<:Model
    
    s = ss.val

    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))

    mm = SVector( (s[i] for i=1:d)... )
    ss_ = SVector( (s[i] for i=(d+1):n)... )

    p = model.source.parameters
    lb = controls_lb(model.source, mm, ss_,  p)
    ub = controls_ub(model.source, mm, ss_,  p)

    r = Ef .⫫ (x-lb)
    # r = (-Ef) .⫫ (ub-x)
    # -r
    
end


function complementarities(model::NoLib.YModel{<:NoLib.MvNormal, B, C, D, E, sS}, ss::NoLib.QP, x::SVector, Ef::SVector) where B where C where D where E where sS<:Model
    

    d = length(NoLib.variables(model.exogenous))
    n = length(NoLib.variables(model.states))
    m0 = zero(SVector{d,Float64})
    
    s = ss.val

    p = model.source.parameters
    lb = controls_lb(model.source, m0, s,  p)
    ub = controls_ub(model.source, m0, s,  p)

    # Ef
    r = Ef .⫫ (x-lb)
    r = (-Ef) .⫫ (ub-x)
    -r

    # r = Ef .⟂ (x-lb)
    # r
    # r = (-Ef) .⟂ (ub-x)
    # -r

end
