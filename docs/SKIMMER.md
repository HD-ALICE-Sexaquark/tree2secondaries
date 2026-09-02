# `skim`

## The Skim Config File

One JSON file per analysis channel, in `configs/`. It holds everything that changes the skim's output: which files to
read, which variables to cache, which preselection to apply, and where the cache goes. It is also the only argument
`skim` takes.

| key                    | required       | default | possible values                                  |
|------------------------|----------------|---------|--------------------------------------------------|
| `channel`              | yes            |  (none) | `ChannelA`, `ChannelD`, `ChannelH`, `LambdaPair` |
| `signal_mass`          | yes            |  (none) | `1.73`,`1.8`,`1.87`,`1.94`,`2.01`,`2.234`        |
| `output`               | yes            |  (none) | path of the cache file, must end in `.root`      |
| `n_events_limit`       | no             |     `0` | events per sample to read (`0` = all)            |
| `weights_pt`           | no             |    `""` | path to a blast-wave weights file                |
| `weights_radius`       | no             |    `""` | path to a material-budget weights file           |
| `samples[]`            | yes, non-empty |  (none) |                                                  |
| `samples[].path`       | yes, full path |  (none) |                                                  |
| `samples[].role`       | no             |  `both` | `both`, `signal`, `background`                   |
| `variables[]`          | yes, non-empty |  (none) | any available in `Skimmer/VariableRegistry.hxx`  |
| `baseline[]`           | no             |    `[]` |                                                  |
| `baseline[].variable`  | yes            |  (none) | any available in `Skimmer/VariableRegistry.hxx`  |
| `baseline[].direction` | yes            |  (none) | `lower`, `upper`, `window`                       |
| `baseline[].value`     | yes            |  (none) |                                                  |

`channel` decides which schema is read, and the names of both input and output RNTuples.

`n_events_limit` caps how many events are read *per sample*. Anything above `0` makes the skim partial, which no
consumer may use for normalization -- see `NInjected` below.

The tuple names come from the schema family, not from the channel, because the three sexaquark channels share a schema
and never end up in the same file. `Hdibaryon` survives as that family name (`Schema::FoundHdibaryon`,
`Cached::Hdibaryon`) even though the channel it serves is called `LambdaPair`.

`variables[]` is a plain array of names. Each must exist in the compiled registry of `Skimmer/VariableRegistry.hxx`
for the declared channel, and none may collide with a reserved cache field (see below). The observable is just one
more of these and has to be listed like any other -- how finely to scan a variable, from where, and in which direction
is the optimizer's business, not the skim's.

Every `baseline[].variable` is cached too, whether or not it also appears in `variables[]`. The cut is applied here
and cannot be loosened downstream, but its surviving distribution stays inspectable -- and a baseline cut nobody can
see the shape of is one nobody can check.

Cut directions are inclusive: `lower` keeps `x >= value`, `upper` keeps `x <= value`, `window` requires `low` and `high`
and keeps `low <= x <= high`.

`samples[].role` has an extra rule: either every sample is `both` or no sample is `both`, and there is at least one
`signal` and at least one `background`.

## The Shape Weight

The dedicated MC injects antisexaquarks flat in pt over `[0, 5)` GeV/c and flat in 2D radius over `[5, 180)` cm.
Neither is physical, so `Weight` reshapes them, using two files named by the config:

| key              | file                                         | holds                                                            |
|------------------|----------------------------------------------|------------------------------------------------------------------|
| `weights_pt`     | `extra/Weights_BlastWave_Pt.root`            | 50 `TH1D` named `"<mass>_<centrality class>"`, e.g. `1.80_10-20` |
| `weights_radius` | `extra/Weights_MaterialBudget_Radius2D.root` | one `TH1D` named `Radius2D`, from real-data photon conversions   |

Both are optional and independent.

The lookup is `template(x)` divided by the flat prior density over the injected range, so the **mean weight over
that range is 1**. That is all `Weight` is: a shape.

Only rows with `Classification == kSignal` are weighted, and only in the three sexaquark channels; every other
row carries exactly `1`. The pt and radius fed to the lookup are MC truth, read off the injected antisexaquark
linked to the candidate.

Two things stop the run rather than being papered over, since either would leave a cache that is wrong without
looking wrong:

- a true-signal candidate with no linked injection;
- a truth pt or radius outside the injected range, which means the templates and the production disagree about
  how the signal was injected.

The centrality classes (`common/DB_Centrality.hpp`) cover `[0, 90)`, which is a real-data acceptance rather than
anything MC respects. An event outside them still has an injection to reweight, so it draws a class uniformly in
percentile, from a fixed seed. The end-of-run summary reports how many rows were weighted and how many of them
drew.

The yield weights, for how many antisexaquarks or (anti)h-dibaryons a real run would actually hold, should be applied
downstream, by the consumer of the cache, on top of this column. The consumer reads the skim config anyway, and both
paths are keys in it, so it can tell a reweighted cache from a flat one.

## Output 1: `RNTuple "Cached{Sexaquark|Hdibaryon}"`

`CachedSexaquark` / `CachedHdibaryon` -- one row per surviving candidate: `SampleIndex`, `Classification`, `RunNumber`,
`DirNumber`, `DirNumberB`, `EventNumber`, `Weight` (see above), then one `float` per cached variable.

`RunNumber`, `DirNumber`, `DirNumberB` and `EventNumber` are copied straight off `POD::Event`; `DirNumberB` is only
meaningful for real data.

`Classification` is the MC-truth label (see `Skimmer/Classification.hxx` and `common/docs/MC_LABELS.md`). Possible
values:
- `kSignal = 0` signal
- `kRealBkg = 1` real background
- `kHybrid = 2` hybrid (none of the above)
- `kReferenceReal = 3` wrong-sign or mixed channel, where no component descends from an injection
- `kReferenceHybrid = 4` wrong-sign or mixed channel, where at least one component descends from injection
- `kUnknown = 5` unknown

`SampleIndex` is the row's position in the config's `samples[]`, and is the join back to `Meta`.

## Output 2: `RNTuple "Meta"`

One row per input file. Fields: `SampleIndex`, `Path`, `Role`, `NInjected`, `NEvents`, `NEventsRead`, `NCandidatesRead`,
`NCandidatesWritten`.

`Role` is stored as the underlying value of `ERole`: `0` signal, `1` background, `2` both.

`NInjected` is the sum of `Injected` signal over every read event. That is exact for a full run, because `t2ds` keeps
every MC event whose `Injected` vector is non-empty -- the events it dropped held no injection to miss. A run with a
non-zero `n_events_limit` stops early and undercounts, which is why such a run must not be used for normalization.

The count is empirical, so it follows the file rather than the declared role: a signal production reused as a background
sample reports its real injections. Sum `NInjected` over the `Role == signal` rows.

`NEvents` comes from the `N_Events` histogram `t2ds` fills once per event, before any candidate is built. It is the
event denominator, and the way to see that an input was truncated -- the RNTuple's entry count is always smaller, since
candidate-less events are dropped.
