# Source attribution

This is a reduced version of the original codebase, for public consumption. (What remains private is messy/unnecessary code.)

The numerical library is an **adapted subset** of [ConsistentKoopman.jl](https://github.com/clairevalva/ConsistentKoopman.jl), based on commit `9c179ee81fe581d032371655ebfc791e4cc5a219`. Its [MIT license](vendor/ConsistentKoopman/LICENSE) is retained. Unused modules, methods, exports, imports, and experimental comments were removed; the active numerical formulas are preserved.

The public drivers use shared paths and paper settings, the supported projection/reconstruction functions, and consistent OSE file conventions. Plotting uses end dates, matches climate fields by physical date, and filters Hovmöller observations along time. Forecast reconstruction retains the paper's full-delay average from a single issue date. 

See [DATA.md](DATA.md) for scientific conventions, input sources, and optional snapshot metadata. This repository illustrates replication of the experiments and does not promise bit-identical archived results.
