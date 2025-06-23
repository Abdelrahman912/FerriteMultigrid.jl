using TimerOutputs

dto = TimerOutput()
reset_timer!(dto)

const liveserver = "liveserver" in ARGS

if liveserver
    using Revise
    @timeit dto "Revise.revise()" Revise.revise()
end

using Documenter, DocumenterCitations, FerriteMultigrid

const is_ci = haskey(ENV, "GITHUB_ACTIONS")

# bibtex_plugin = CitationBibliography(
#     joinpath(@__DIR__, "src", "assets", "references.bib"),
#     style=:numeric
# )

# Build documentation.
@timeit dto "makedocs" makedocs(
    format = Documenter.HTML(
        assets = [
            "assets/custom.css",
            "assets/citations.css",
            # "assets/favicon.ico"
        ],
        # canonical = "https://localhost/",
        collapselevel = 1,
    ),
    sitename = "FerriteMultigrid.jl",
    doctest = false,
    warnonly = true,
    draft = liveserver,
    pages = Any[
        "Home" => "index.md",
        
        ],
    # plugins = [
    #     bibtex_plugin,
    # ]
)

print_timer(dto)
