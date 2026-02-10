using SpatialBernoulli
using Documenter 

DocMeta.setdocmeta!(
    SpatialBernoulli, :DocTestSetup, :(using SpatialBernoulli); recursive=true
)

makedocs(;
    modules=[SpatialBernoulli],
    authors="Caroline Cognot <caroline.cognot@agroparistech.fr>",
    sitename="SpatialBernoulli.jl",
    format=Documenter.HTML(;
        canonical="https://caroline-cognot.github.io/SpatialBernoulli.jl", edit_link="main"
    ),
    pages=["Home" => "index.md" ,
    "Example" => "example.md",
    "Functions"=> "functions.md"],

)

deploydocs(; repo="github.com/caroline-cognot/SpatialBernoulli.jl.git")

