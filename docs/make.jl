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
    format = Documenter.HTML(),
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

# Deploy built documentation (only if not liveserver)
if !liveserver
    @timeit dto "deploydocs" deploydocs(
        repo = "github.com/Abdelrahman912/FerriteMultigrid.jl.git",
        devbranch = "main",  
        push_preview = true
    )
end


print_timer(dto)
