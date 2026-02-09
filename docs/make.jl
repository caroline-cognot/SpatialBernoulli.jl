using SpatialBernoulli
using Documenter , Example
DocMeta.setdocmeta!(
    SpatialBernoulli, :DocTestSetup, :(using SpatialBernoulli); recursive=true
)

makedocs(;
    modules=[SpatialBernoulli],
    authors="Caroline Cognot <caroline.cognot@agroparistech.fr>",
    sitename="SpatialBernoulli.jl",
    format=Documenter.HTML(;
        canonical="https://caroline-cognot
.github.io/SpatialBernoulli.jl", edit_link="main", assets=String[]
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/caroline-cognot/SpatialBernoulli.jl", devbranch="main")


# Only deploy on GitHub Actions
if get(ENV, "GITHUB_ACTIONS", "false") == "true"
    deploydocs(
        repo="github.com/caroline-cognot/SpatialBernoulli.jl",
        devbranch="main"
    )
end