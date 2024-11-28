using DoloYAML
import Dolo

@testset "Testing yaml_import" begin
    path = joinpath(DoloYAML.pkg_path, "examples")
    files = [f for f in readdir(path) if occursin("yaml", f)]
    for fname in files
        println("Importing ", fname)
        model = DoloYAML.yaml_import(joinpath(path, fname))
        mmodel = DoloYAML.orphan(model)
        @assert isbits(mmodel)
        dmodel = Dolo.discretize(model)
        # sol = Dolo.time_iteration(model);
    end
end
