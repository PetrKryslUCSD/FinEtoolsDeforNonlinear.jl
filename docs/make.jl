using Documenter, FinEtools, FinEtoolsDeforLinear, FinEtoolsDeforNonlinear

makedocs(
	modules = [FinEtoolsDeforNonlinear],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsDeforNonlinear.jl",
	pages = Any[
	"Home" => "index.md",
	"Guide" => "guide/guide.md",
	"Types and Functions" => Any[
		"man/types.md",
		"man/functions.md"]
		]
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl.git",
)
