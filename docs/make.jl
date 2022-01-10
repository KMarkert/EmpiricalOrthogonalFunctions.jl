push!(LOAD_PATH, "../src/")

using Documenter, EmpiricalOrthoFuncs

pages = ["Home" => "index.md", "API" => "api.md"]

makedocs(;
    modules = [EmpiricalOrthoFuncs],
    authors = "Kel Markert",
    repo = "https://github.com/KMarkert/EmpiricalOrthoFuncs.jl/blob/{commit}{path}#L{line}",
    sitename = "EmpiricalOrthoFuncs.jl",
    pages = pages,
)

deploydocs(;
    repo = "github.com/KMarkert/EmpiricalOrthoFuncs.jl.git",
    devbranch = "main",
)
