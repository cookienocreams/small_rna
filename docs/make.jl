using Documenter, small_rna, DocumenterTools

makedocs(sitename="NEXTFLEX Small RNA Anlaysis Documentation")

deploydocs(
    repo = "github.com/cookienocreams/small_rna.git"
    , branch = "docs"
    , devbranch = "main"
)
