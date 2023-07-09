# DoloYAML

Import Dolo models from a YAML file.

To import one of the example models:

```
using DoloYAML

model = DoloYAML.yaml_import("examples/rbc_mc.yaml")

# recover steadystate vectors
m,s,x,p = model.calibration[:exogenous, :states, :controls, :parameters]

# compute state transition according to equation defined in file
S = DoloYAML.transition(m,s,x,m,p)
```

For the accepted syntax, see the Dolo [documentation](http://www.econforge.org/Dolo.jl/latest/index.html).

## Developer note

DoloYAML.jl was formerly part of the Dolo.jl package. It is packaged out in the hope of reducing dependency complexity.