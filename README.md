# tree2secondaries

Second step.
ROOT-based standalone application that groups ROOT trees in a simple vector format into V0s.

## Requirements

- [ROOT](https://root.cern.ch)

## Building

```
mkdir build && cd build
cmake ../ <options>
cmake --build .
```

Additional `<options>`:

* `-DCMAKE_EXPORT_COMPILE_COMMANDS=1` -- export `compile_commands.json` file
* `-DT2S_DEBUG` -- enable debug messages
* `-DENABLE_PROFILING=ON` -- enable profiling (see below)

## Usage

```
./src/App [OPTIONS] SUBCOMMAND [SUBOPTIONS]

OPTIONS:
  -h,--help                 Print help message and exit
  -i,--input [REQUIRED]     Path(s) of input file(s)
  -o,--output               Path of output file
  -n,--nevents              Limit to N events
SUBCOMMANDS:
  data                      Process data
  mc                        Process MC
    SUBOPTIONS:
    -s,--signal             Process Signal MC
    -c,--channel {A,D,E,H}  Process a standard reaction channel
    -a,--anti    {A,D,E,H}  Process an anti-reaction channel
```

## Debugging

```
root -l -b -q caca.root -e 'ProEvents->Scan("TMath::Sqrt(V0_E*V0_E-V0_Px*V0_Px-V0_Py*V0_Py-V0_Pz*V0_Pz):TMath::Sqrt(V0_Xv*V0_Xv+V0_Yv*V0_Yv)", "V0_PID == -3122")'
root caca.root -e 'ProEvents->Draw("TMath::Sqrt(V0_E*V0_E-V0_Px*V0_Px-V0_Py*V0_Py-V0_Pz*V0_Pz)", "V0_PID == -3122")'
root -l -b -q SexaquarkResults_SignalMC.root -e 'Events->Scan("TMath::Sqrt(K0S_E*K0S_E-K0S_Px*K0S_Px-K0S_Py*K0S_Py-K0S_Pz*K0S_Pz):TMath::Sqrt(K0S_Xv*K0S_Xv+K0S_Yv*K0S_Yv):K0S_IsSignal")'
```

```
 > current_D.txt ; \
root -l -b -q SexaquarkResults_SignalMC.root -e \
'Events->Scan("L_Neg_EsdIdx:L_Pos_EsdIdx:L_Mass:L_Radius:L_IsSignal:L_Px:L_Py:L_Pz:L_Xv:L_Yv:L_Zv:L_Neg_Px:L_Neg_Py:L_Neg_Pz:L_Pos_Px:L_Pos_Py:L_Pos_Pz")' >> current_D.txt ; \
root -l -b -q SexaquarkResults_SignalMC.root -e \
'Events->Scan("K0S_Neg_EsdIdx:K0S_Pos_EsdIdx:K0S_Mass:K0S_Radius:K0S_IsSignal:K0S_Px:K0S_Py:K0S_Pz:K0S_Xv:K0S_Yv:K0S_Zv:K0S_Neg_Px:K0S_Neg_Py:K0S_Neg_Pz:K0S_Pos_Px:K0S_Pos_Py:K0S_Pos_Pz")' >> current_D.txt ; \
root -l -b -q SexaquarkResults_SignalMC.root -e \
'Events->Scan("ASA_V0a_Idx:ASA_V0b_Idx:ASA_Xv:ASA_Yv:ASA_Mass:ASA_Pt:ASA_Radius:ASA_IsSignal:ASA_ReactionID:ASA_IsHybrid")' >> current_D.txt ; \
root -l -b -q SexaquarkResults_SignalMC.root -e \
'Events->Scan("BA_V0a_Idx:BA_V0b_Idx:BA_Xv:BA_Yv:BA_Mass:BA_Pt:BA_Radius:BA_IsSignal:BA_ReactionID:BA_IsHybrid")' >> current_D.txt ; \
root -l -b -q SexaquarkResults_SignalMC.root -e \
'Events->Scan("ASD_V0_Idx:ASD_Ba_Entry:ASD_Xv:ASD_Yv:ASD_Mass:ASD_Pt:ASD_Radius")' >> current_D.txt ; \
root -l -b -q SexaquarkResults_SignalMC.root -e \
'Events->Scan("BD_V0_Idx:BD_Ba_Entry:BD_Xv:BD_Yv:BD_Mass:BD_Pt:BD_Radius")' >> current_D.txt
```

## Memory Leaks

```
valgrind --leak-check=full --suppressions=$(root-config --etcdir)/valgrind-root.supp APP &> LOG
```

## Profiling

```
mkdir profile && cd build
cmake ../ -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DENABLE_PROFILING=ON
cmake --build .
cd apps/
./App [options] # will generate gmon.out
gprof App > gprof.log
```
