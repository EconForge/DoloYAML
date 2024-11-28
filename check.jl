import DoloYAML

model = DoloYAML.yaml_import("examples/consumption_savings_iid.yaml");
mmodel = DoloYAML.orphan(model)
isbits(mmodel)
dmodel = Dolo.discretize(model)
sol = Dolo.time_iteration(model);


model = DoloYAML.yaml_import("examples/rbc_iid.yaml");
mmodel = DoloYAML.orphan(model)
isbits(mmodel)
dmodel = Dolo.discretize(model)
sol = Dolo.time_iteration(model);




model = DoloYAML.yaml_import("examples/rbc_mc.yaml");
mmodel = DoloYAML.orphan(model)
isbits(mmodel)
dmodel = Dolo.discretize(model)
sol = Dolo.time_iteration(model);

model = DoloYAML.yaml_import("examples/sudden_stop.yaml");
mmodel = DoloYAML.orphan(model)
isbits(mmodel)
dmodel = Dolo.discretize(model)
sol = Dolo.time_iteration(model);


model = DoloYAML.yaml_import("examples/neoclassical.yaml");
mmodel = DoloYAML.orphan(model)
isbits(mmodel)
dmodel = Dolo.discretize(model)
sol = Dolo.time_iteration(model,verbose=false);

model = DoloYAML.yaml_import("examples/rbc_ar1.yaml");
mmodel = DoloYAML.orphan(model)
isbits(mmodel)
dmodel = Dolo.discretize(model)
sol = Dolo.time_iteration(model);