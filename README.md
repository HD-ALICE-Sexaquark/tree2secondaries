# tree2secondaries

## Requirements

- CMake (v3.28 or higher)
- C++ compiler compatible with C++23
- Internet connection to fetch **[CLIUtils/CLI11] (<https://github.com/CLIUtils/CLI11>)** and
**[Eigen](https://gitlab.com/libeigen/eigen)**
- **[ROOT](https://root.cern.ch)** v6.40.00+

## Building

```bash
mkdir -p build && cd build
cmake ../ <options>
cmake --build .
```

Additional `<option>`:

- `-DCMAKE_BUILD_TYPE=` -- (`Debug`, `Release`, `DebWithRelInfo`) if not specified, it defaults to `Release`
- `-DENABLE_PROFILING=ON` -- (default: `OFF`) enable profiling (see below)

## Usage

```
./t2ds [OPTIONS] SUBCOMMAND

OPTIONS:

  -h, --help                  Print this help message and exit
  -i, --input TEXT REQUIRED   Path of input file
  -o, --output TEXT           Path of output file
  -n, --nevents NUMBER        Limit to N events

SUBCOMMANDS:

  pack data
  pack mc   : Read "AnalysisResults.root" to package injected anti-sexaquarks, tracks and V0s.
              Both cases require option:
              -c, --channel : {A,D,H} Process a standard reaction channel
              And `pack mc` additionaly requires:
              -m, --mass : {1.73,1.8,1.87,1.94,2.01} Assign injected anti-sexaquark mass

  search data
  search mc   : Read "PackedRNT.root" files to search for anti-sexaquark reactions.
                Both cases require option:
                -c, --channel : {A,D,H} Process a standard reaction channel

  verify data
  verify mc   : Read "AnalysisResults.root" to verify the existence of (anti)h-dibaryons
```
