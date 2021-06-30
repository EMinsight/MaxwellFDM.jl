module MaxwellFDFD

using Reexport
@reexport using MaxwellBase
using AbbreviatedTypes

# The order of inclusion matters: if types or functions in file A are used in file B, file A
# must be included first.
include("source/source.jl")
# include("equation.jl")
include("model/model.jl")

end # module MaxwellFDFD
