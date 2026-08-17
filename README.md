# tree2secondaries

Single-threaded program, ideally executed in a job scheduler like Slurm or HTCondor.

## Requirements

- CMake (v3.28 or higher)
- C++ compiler compatible with C++23
- Internet connection to fetch **[CLIUtils/CLI11] (<https://github.com/CLIUtils/CLI11>)** and
**[Eigen](https://gitlab.com/libeigen/eigen)**
- **[ROOT](https://root.cern.ch)** v6.40.02+

## Building

```bash
mkdir -p build && cd build
cmake ../ <options>
cmake --build .
```

Additional `<options>`:

- `-DCMAKE_BUILD_TYPE=` -- (`Debug`, `Release`, `DebWithRelInfo`) if not specified, it defaults to `Release`
- `-DENABLE_PROFILING=ON` -- (default: `OFF`) enable profiling (see below)

## Usage

```
./t2ds [OPTIONS] SUBCOMMAND [SUBCOMMAND OPTIONS]

OPTIONS:

  -h, --help                  Print this help message and exit
  -i, --input FILE1 FILE2 ... [REQUIRED] Path(s) of input file(s). Faulty files get skipped.
  -o, --output TEXT           Path of output file
  -n, --nevents NUMBER        Limit to N events, counted across all input files

SUBCOMMANDS:

  pack data
  pack mc   : Read "AnalysisResults.root" to package injected anti-sexaquarks, tracks and V0s.
              Additionally, `pack mc` requires:
              -c, --channel : {A,D,H} Process a standard reaction channel
              -m, --mass : {1.73,1.8,1.87,1.94,2.01} Assign injected anti-sexaquark mass

  search data
  search mc   : Read "PackedRNT.root" files to search for anti-sexaquark reactions.
                Both cases require option:
                -c, --channel : {A,D,H} Process a standard reaction channel

  verify data
  verify mc   : Read "AnalysisResults.root" to verify the existence of (anti)h-dibaryons.
                Doesn't require further options.
```
