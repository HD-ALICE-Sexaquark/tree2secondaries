# `skim`

## The Skim Config File

One JSON file per analysis channel, in `configs/skim/`. It holds only what changes the skim's output: which files to
read, which variables to cache, and which preselection to apply.

| key                   | required       | default | possible values                                  |
|-----------------------|----------------|---------|--------------------------------------------------|
| `channel`             | yes            |  (none) | `ChannelA`, `ChannelD`, `ChannelH`, `LambdaPair` |
| `signal_mass`         | yes            |   (none)| `1.73`,`1.8`,`1.87`,`1.94`,`2.01`,`2.234`        |
| `samples[]`           | yes, non-empty |  (none) |                                                  |
| `samples[].path`      | yes, full path |  (none) |                                                  |
| `samples[].role`      | no             |  `both` | `both`, `signal`, `background`                   |
| `variables[]`         | yes, non-empty |  (none) | any available in `Skimmer/VariableRegistry.hxx`  |
| `baseline[]`          | no             |    `[]` |                                                  |
| `baseline[].variable` | yes            |  (none) | any available in `Skimmer/VariableRegistry.hxx`  |
| `baseline[].direction` | yes           |  (none) | `lower`, `upper`, `window`                       |
| `baseline[].value`    | yes            |  (none) |                                                  |

`channel` decides everything else about the file layout: which schema is read, which RNTuple that schema lives in, and
which one the cache is written to.

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

## Output 1: `RNTuple "Cached{Sexaquark|Hdibaryon}"`

`CachedSexaquark` / `CachedHdibaryon` -- one row per surviving candidate: `SampleIndex`, `Classification`, `RunNumber`,
`DirNumber`, `EventNumber`, `Weight`, then one `float` per cached variable.

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
every MC event whose `Injected` vector is non-empty -- the events it dropped held no injection to miss. A run with `-n`
stops early and undercounts, which is why such a run is flagged partial and must not be used for normalization.

The count is empirical, so it follows the file rather than the declared role: a signal production reused as a background
sample reports its real injections. Sum `NInjected` over the `Role == signal` rows.

`NEvents` comes from the `N_Events` histogram `t2ds` fills once per event, before any candidate is built. It is the
event denominator, and the way to see that an input was truncated -- the RNTuple's entry count is always smaller, since
candidate-less events are dropped.

## Output 3: `TObjString "CacheSource"`

It hold one JSON object:

```json
{
  "channel": "LambdaPair",
  "signal_mass": 2.234,
  "ntuple": "CachedHdibaryon",
  "fields": ["Mass", "Rapidity", "Pt", "..."],
  "config_path": "../configs/hdibaryon.json",
  "is_partial": false,
  "n_events_limit": 0
}
```
