# Run `julia --project=docs docs/make.jl` from the top directory to execute this script.
# Run `julia --project=docs docs/serve.jl` to build and serve the docs locally.

using Scri
using Documenter
using DocumenterCitations
using DocumenterInterLinks
using Revise

# When run via `docs/serve.jl`, this ensures that any changes in the docstrings are picked
# up without needing to restart the server.
Revise.revise()

DocMeta.setdocmeta!(Scri, :DocTestSetup, :(using Scri); recursive=true)

include("local_notes.jl")

bibliography = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib"); style=:numeric
)

# Add ability to link to other packages' docs.  The keys are the names of packages used in
# external references by this package's docs, and the values are the base URLs of the
# documentation for those packages.  Use these as, e.g., [`Quaternionic.rotor`](@extref).
links = InterLinks(
    "Quaternionic" => "https://moble.github.io/Quaternionic.jl/stable/",
    "SphericalFunctions" => "https://moble.github.io/SphericalFunctions.jl/stable/",
    "Julia" => "https://docs.julialang.org/en/v1/",
)

# Add titles of sections and overrides page titles
const titles = Dict(
    # "10-tutorials" => "Tutorials", # example folder title
    "80-details" => "Details",
    "91-developer.md" => "Developer docs",
    "65-api.md" => "API",
    "99-local_notes" => "Notes",
)

function recursively_list_pages(folder; path_prefix="")
    pages_list = Any[]
    for file ∈ readdir(folder)
        if file == "index.md"
            # We add index.md separately to make sure it is the first in the list
            continue
        end
        # this is the relative path according to our prefix, not @__DIR__, i.e., relative to `src`
        relpath = joinpath(path_prefix, file)
        # full path of the file
        fullpath = joinpath(folder, relpath)

        if isdir(fullpath)
            # If this is a folder, enter the recursion case
            subsection = recursively_list_pages(fullpath; path_prefix=relpath)

            # Ignore empty folders
            if length(subsection) > 0
                title = if haskey(titles, relpath)
                    titles[relpath]
                else
                    @error "Bad usage: '$relpath' does not have a title set. Fix in 'docs/make.jl'"
                    relpath
                end
                push!(pages_list, title => subsection)
            end

            continue
        end

        if splitext(file)[2] != ".md" # non .md files are ignored
            continue
        elseif haskey(titles, relpath) # case 'title => path'
            push!(pages_list, titles[relpath] => relpath)
        else # case 'title'
            push!(pages_list, relpath)
        end
    end

    return pages_list
end

function list_pages()
    root_dir = joinpath(@__DIR__, "src")
    pages_list = recursively_list_pages(root_dir)

    return ["index.md"; pages_list]
end

makedocs(;
    modules=[Scri],
    authors="Michael Boyle <michael.oliver.boyle@gmail.com>",
    # repo="https://github.com/moble/Scri.jl/blob/{commit}{path}#{line}",
    repo=Documenter.Remotes.GitHub("moble", "Scri.jl"),
    sitename="Scri.jl",
    format=Documenter.HTML(;
        prettyurls=(!("local" in ARGS)),  # Use clean URLs, unless built as a "local" build
        edit_link="main",  # Link out to "main" branch on github
        canonical="https://moble.github.io/Scri.jl",
        assets=String["assets/citations.css", "assets/custom.css"],
    ),
    remotes=notes_remotes,
    pages=[list_pages()..., notes_pages...],
    plugins=[bibliography, links],
)

deploydocs(; repo="github.com/moble/Scri.jl")
