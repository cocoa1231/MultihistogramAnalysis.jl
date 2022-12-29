using MultihistogramAnalysis
using Documenter

DocMeta.setdocmeta!(MultihistogramAnalysis, :DocTestSetup, :(using MultihistogramAnalysis); recursive=true)

makedocs(;
    modules=[MultihistogramAnalysis],
    authors="Cocoa <cocoathepenguin@protonmail.com> and contributors",
    repo="https://github.com/cocoa1231/MultihistogramAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="MultihistogramAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cocoa1231.github.io/MultihistogramAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cocoa1231/MultihistogramAnalysis.jl",
    devbranch="main",
)
