# Tree2DoubleStrangeness

Single-threaded programs, ideally executed in a job scheduler like Slurm or HTCondor.

`t2ds` parses output from [**E2T**](https://github.com/HD-ALICE-Sexaquark/esd2tree). It has two modes:
- `Finder`: given a collection of tracks, form anti-sexaquark candidates
- `Verifier`: given a collection of previously found (anti)lambdas, form resonance-like (anti)h-dibaryons
The output is stored in event-based RNTuples. It can parse real data or MC.

`skim` flattens the candidate vectors of one or more `t2ds` RNTuples into a scalar-only cache, applying the config's
baseline preselection on the way. Only works for MC productions. Its only argument is a skim config file (`configs/`),
which names the input files, the variables to cache and where the cache goes.

## Requirements

- CMake (v3.28 or higher)
- C++ compiler compatible with C++23
- [**ROOT**](https://root.cern.ch) v6.40.02+
- Internet connection to fetch:
  - [**CLIUtils/CLI11**](https://github.com/CLIUtils/CLI11)
  - [**Eigen**](https://gitlab.com/libeigen/eigen)
  - [**nlohmann/json**](https://github.com/nlohmann/json) (only needed by `skim`, and only fetched when the ROOT
  installation wasn't built with `-Dbuiltin_nlohmannjson=ON` and no system copy is found)

## Building

```bash
git clone --recurse-submodules git@github.com:HD-ALICE-Sexaquark/tree2secondaries.git <repo dir>
# alternatively:
#     git clone <placeholder for remote url> <repo dir>
#     git submodule update --init
cd <repo dir>
mkdir -p build && cd build
cmake ../ <options>
cmake --build .
# alternatively, build split targets:
#     cmake --build . --target t2ds
#     cmake --build . --target skim
```

Additional `<options>`:

- `-DCMAKE_BUILD_TYPE=` -- choose between `Release` (default), `DebWithRelInfo`, `Debug` (enables directive
`T2DS_DEBUG`)
- `-DENABLE_PROFILING=ON` -- (default: `OFF`) enable profiling

## Usage: T2DS

```
./t2ds [OPTIONS] SUBCOMMAND [SUBCOMMAND OPTIONS]

OPTIONS:

  -h, --help                  Print this help message and exit
  -i, --input FILE1 FILE2 ... [REQUIRED] Path(s) of input file(s). Faulty files get skipped.
  -o, --output TEXT           Path of output file
  -n, --nevents NUMBER        Limit to N events, counted across all input files

SUBCOMMANDS:

  find data
  find mc   : Read "AnalysisResults.root" to search for antisexaquark reactions via three channels and their charged-
              conjugates as reference background.
              Additionally, `find mc` requires:
              -m, --mass : {1.73,1.8,1.87,1.94,2.01} Assign injected antisexaquark mass

  verify data
  verify mc   : Read "AnalysisResults.root" to verify the existence of (anti)h-dibaryons.
                It also produces artificial lambda+antilambda background.
                Doesn't require further options.
```

More information: [docs/FINDER.md](./docs/FINDER.md) / [docs/VERIFIER.md](./docs/VERIFIER.md)

## Usage: Skimmer

```
./skim [OPTIONS] <config>

POSITIONALS:
  <config>      [REQUIRED] Path of JSON config file

OPTIONS:
  -h, --help    Print this help message and exit
```

More information: [docs/SKIMMER.md](./docs/SKIMMER.md)
