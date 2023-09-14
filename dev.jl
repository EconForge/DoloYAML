import DoloYAML
import Dolo

model = DoloYAML.yaml_import("integration_A.yaml")
dmodel = Dolo.discretize(model)

@time Dolo.time_iteration(model)
@time Dolo.time_iteration(model, improve=true)



using StaticArrays
# Dolo.time_iteration(model; improve=true)
s0 = Dolo.calibrated(model, :states)
x0 = Dolo.calibrated(model, :controls)
m0 = Dolo.calibrated(model, :exogenous)
p0 = Dolo.calibrated(model, :parameters)
vals = (;s0, x0, m0, p0)
ss = Dolo.QP((5,200), s0)


ss0 = Dolo.QP((5, 200), Dolo.calibrated(model, :states))
x0 = Dolo.calibrated(model, :controls)
ϕ0  = Dolo.Policy(model.states, model.controls, u->x0)

Dolo.F(dmodel, ss0, x0, ϕ0)

xvec = Dolo.GVector(dmodel.grid, [x0 for i=1:length(dmodel.grid)])
φ1 = Dolo.DFun(model, xvec)

res = Dolo.F(dmodel, xvec, φ1);
dres = Dolo.dF_1(dmodel, xvec, φ1);

import Dolo: CPU

res1 = res*0.0
@time begin res1 = res*0.0;(Dolo.F!(res1, dmodel, xvec, φ1);println(sum(sum(res1)))) end
@time begin res2 = res*0.0;Dolo.F!(res2, dmodel, xvec, φ1, CPU());println(sum(sum(res2))) end


function test(dmodel,ss0, x0, ϕ0)
    return sum(Dolo.dF_1(dmodel,ss0, x0, ϕ0))
end
test(dmodel,ss0, x0, ϕ0)

@time Dolo.dF_1!(dres1, dmodel, xvec, φ1);

@time begin dres1 = dres*0.0;(Dolo.dF_1!(dres1, dmodel, xvec, φ1);println(sum(sum(dres1)))) end
@time begin dres2 = dres*0.0;Dolo.dF_1!(dres2, dmodel, xvec, φ1, CPU());println(sum(sum(dres2))) end

@time Dolo.time_iteration(model;)

@time Dolo.time_iteration(model; engine=:cpu)



# x0 = Dolo.GVector(dmodel.grid, [x0 for i=1:length(dmodel.grid)])
(;φ) = Dolo.time_iteration_workspace(dmodel; )
# ϕ0  = Dolo.Policy(model.states, model.controls, s->x0)


function test(dmodel, vals, φ)

    model = dmodel.model

    (;m0, s0, x0, p0) = vals

    # r = Dolo.arbitrage(model.source, m0, s0, x0, m0, s0, x0, p0) # ok
    # s = Dolo.transition(model.source, m0, s0, x0, m0, p0) # ok
    
    ss = Dolo.QP((5,200), s0)
    # SS = Dolo.QP((5,SVector(s0[3],s0[4])), s0)

    Dolo.F(dmodel, ss, x0, φ) # ok


    # if abs(r)>-0.1
    #     return nothing
    # else
    #     return "No"
    # end
    # return sum(r)

end

test(dmodel, vals, φ)

@time test(dmodel, vals, φ);





function test2(dmodel, φ, xvec, rvec) # OK!

    model = dmodel.model
    Dolo.F!(rvec, dmodel, xvec, φ) 

end

xvec = Dolo.GVector(dmodel.grid, [x0 for i=1:length(dmodel.grid)])
rvec = deepcopy(xvec)

test2(dmodel, φ, xvec, rvec)
@time test2(dmodel, φ, xvec, rvec)

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

