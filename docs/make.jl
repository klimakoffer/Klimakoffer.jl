using Documenter
import Pkg
using Klimakoffer

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(Klimakoffer, :DocTestSetup, :(using Klimakoffer); recursive=true)

# Make documentation
makedocs(
    # Specify modules for which docstrings should be shown
    modules = [Klimakoffer],
    # Set sitename to Klimakoffer
    sitename="Klimakoffer.jl",
    # Provide additional formatting options
    format = Documenter.HTML(
        # Disable pretty URLs during manual testing
        prettyurls = get(ENV, "CI", nothing) == "true",
        # Explicitly add favicon as asset
        # assets = ["assets/favicon.ico"],
        # Set canonical URL to GitHub pages URL
        canonical = "https://klimakoffer.github.io/Klimakoffer.jl/dev"
    ),
    # Explicitly specify documentation structure
    pages = [
        "Home" => "index.md",
        "Reference" => "reference.md",
        "License" => "license.md"
    ],
    strict = true # to make the GitHub action fail when doctests fail, see https://github.com/neuropsychology/Psycho.jl/issues/34
)

deploydocs(
    repo = "github.com/klimakoffer/Klimakoffer.jl",
    devbranch = "main",
    push_preview = true
)
