using DoloYAML

@testset "Testing yaml_import" begin
    path = joinpath(DoloYAML.pkg_path, "examples")
    files = [f for f in readdir(path) if occursin("yaml", f)]
    for fname in files
        println("Importing ", fname)
        model = DoloYAML.yaml_import(joinpath(path, fname))
    end
end
