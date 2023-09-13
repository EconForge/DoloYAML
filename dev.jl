import DoloYAML
import Dolo

model = DoloYAML.yaml_import("integration_A.yaml")
dmodel = Dolo.discretize(model)

using StaticArrays
# Dolo.time_iteration(model; improve=true)

# s0 = Dolo.QP((5, 200), Dolo.calibrated(model, :states))
# x0 = Dolo.calibrated(model, :controls)
# ϕ0  = Dolo.Policy(model.states, model.controls, s->x0)
# Dolo.F(dmodel, s0, x0, ϕ0)

function testalloc(dmodel)

    model = dmodel.model

    # s0 = Dolo.QP((5, 200), Dolo.calibrated(model, :states))
    s0 = Dolo.calibrated(model, :states)
    x0 = Dolo.calibrated(model, :controls)
    m0 = Dolo.calibrated(model, :exogenous)
    p0 = Dolo.calibrated(model, :parameters)
    print(s0, x0, m0, p0)

    r = Dolo.arbitrage(model, m0, s0, x0, m0, s0, x0, p0)

    return r
    
    # ϕ0  = Dolo.Policy(model.states, model.controls, s->x0)

    # r = Dolo.F(dmodel, s0, x0, ϕ0)
    
    # return sum(r)

    # r = tuple( (S for (S) in  Dolo.τ(model, s0, x0))...)

    # M = SVector(0.0, 0.0, 0.0)
    # S = Dolo.transition(model, s0, x0)

    # return sum(S.val)

    # r = tuple( (w for (w, S) in  Dolo.τ(dmodel, s0, x0))...)

    # return sum(r)

end

testalloc(dmodel)

@time testalloc(dmodel)



model = DoloYAML.yaml_import("integration_B.yaml")
dmodel = Dolo.discretize(model)

Dolo.time_iteration(model)

s0 = Dolo.QP((5, 4001), Dolo.calibrated(model, :states))
x0 = Dolo.calibrated(model, :controls)
ϕ0  = Dolo.Policy(model.states, model.controls, s->x0)
Dolo.F(dmodel, s0, x0, ϕ0)




function testalloc(dmodel, s0, x0,ϕ0)
   r = Dolo.F(dmodel, s0, x0, ϕ0)
   return sum(r)
end

@time testalloc(dmodel, s0, x0, ϕ0)

