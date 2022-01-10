push!(LOAD_PATH, "../src/")

using Documenter, EmpiricalOrthogonalFunctions

pages = ["Home" => "index.md", "API" => "api.md"]

makedocs(;
    modules = [EmpiricalOrthogonalFunctions],
    authors = "Kel Markert",
    repo = "https://github.com/KMarkert/EmpiricalOrthogonalFunctions.jl/blob/{commit}{path}#L{line}",
    sitename = "EmpiricalOrthogonalFunctions.jl",
    pages = pages,
)

deploydocs(;
    repo = "github.com/KMarkert/EmpiricalOrthogonalFunctions.jl.git",
    devbranch = "main",
)
