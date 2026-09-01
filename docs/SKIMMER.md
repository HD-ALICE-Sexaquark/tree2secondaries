# `skim`

The limit is per sample, not across all of them: counted globally it would consume the samples in config order and stop,
which for a config with roles can mean "all of the signal, none of the background". A run with `-n` is flagged as
partial in `Provenance` and must not be used for normalization.

### The Config File

One JSON file per analysis channel, in `configs/`, contains cut values, a scan range, or normalization assumptions.

| key                                        | required            | default |
|--------------------------------------------|---------------------|---------|
| `channel`                                  | yes                 |  (none) |
| `hypothesis`                               | yes                 |  (none) |
| `signal_mass`                              | no                  |   `0.0` |
| `samples[]`                                | yes, non-empty      |  (none) |
| `samples[].path`                           | yes, full path      |  (none) |
| `samples[].ntuple`                         | yes                 |  (none) |
| `samples[].run_number`                     | no                  |     `0` |
| `samples[].role`                           | no                  |  `both` |
| `samples[].n_injected_per_event`           | no                  |     `0` |
| `observable.variable`                      | yes                 |  (none) |
| `observable.bins`                          | no                  |   `200` |
| `observable.range`                         | yes                 |  (none) |
| `baseline[]`                               | no                  |    `[]` |
| `baseline[].variable`                      | yes                 |  (none) |
| `baseline[].direction`                     | yes                 |  (none) |
| `baseline[].value`                         | yes                 |  (none) |
| `variables[]`                              | yes, non-empty      |  (none) |
| `variables[].name`                         | yes                 |  (none) |
| `variables[].direction`                    | yes                 |  (none) |
| `variables[].range`                        | yes                 |  (none) |
| `variables[].steps`                        | no, at least `2`    |   `100` |
| `variables[].initial`                      | yes                 |  (none) |
| `variables[].in_grid`                      | no                  | `false` |
| `fom.formula`                              | yes                 |  (none) |
| `fom.a`                                    | no                  |   `3.0` |
| `fom.f_syst`                               | no                  |   `0.2` |
| `fom.n_signal_expected`                    | see below           |  (none) |
| `fom.signal_yield`                         | see below           |  (none) |
| `fom.signal_yield.dndy`                    | yes, within block   |  (none) |
| `fom.signal_yield.delta_y`                 | yes, within block   |  (none) |
| `fom.signal_yield.source`                  | yes, within block   |  (none) |
| `fom.signal_yield.branching_ratio`         | no                  |   `1.0` |
| `fom.signal_yield.n_injected_species`      | no                  |     `1` |
| `fom.signal_yield.interaction_probability` | no                  |   `1.0` |
| `guards.min_raw_signal`                    | no                  |    `20` |
| `guards.min_raw_background`                | no                  |    `20` |
| `sentinel_ok[]`                            | no                  |    `[]` |

Cut directions are inclusive: `lower` keeps `x >= value`, `upper` keeps `x <= value`, `window` keeps `low <= x <= high`.

A `window` cut must give its value as a two-element array. This is enforced rather than defaulted, so that a window
written with a single number fails loudly instead of silently acquiring an upper edge of zero.

`fom.formula` is one of `punzi`, `asimov` or `poisson`. The latter two are built on an absolute signal yield and so
require either `n_signal_expected` (the number) or a `signal_yield` block (the recipe for it) -- exactly one of the two,
never both, so they cannot disagree. Punzi is the default precisely because it needs neither: only the signal
efficiency, which is measurable without assuming a production rate.

`sentinel_ok[]` names variables whose dummy population is acknowledged rather than accidental. Every name must be a
scanned variable or a baseline cut.

`samples[].role` is one of `signal`, `background` or `both`, and the whole `samples[]` array must be consistent:
- either every sample is `both` -- one production serving as its own background,
- or no sample is `both`, and there is at least one `signal` and at least one `background`.

`samples[].n_injected_per_event` overrides the channel's compile-time count (20 reactions per event for the sexaquark
channels, 100 (anti)h-dibaryons for `Hdibaryon`). `0` means "use the channel's", except on a `background` sample, which
by definition carries no injected signal and so resolves to `0`. The resolved value is what lands in `Meta`.

### The Output Cache

One file holding three RNTuples:

- **`<Channel>`** -- one row per surviving candidate: `SampleIndex`, `Process`, `RunNumber`, `EventNumber`, `Weight`,
  then one `float` per cached variable.
- **`Meta`** -- one row per input file: `SampleIndex`, `Path`, `RunNumber`, `Role`, `NInjectedPerEvent`, `NEvents`,
  `NEventsRead`, `NCandidatesRead`, `NCandidatesWritten`.
- **`Provenance`** -- one row: `Channel`, `Hypothesis`, `ConfigPath`, `Observable`, `SignalMass`,
  `NExpectedEventsInRealData`, `IsPartial`.

`SampleIndex` is the row's position in the config's `samples[]`, and is the join back to `Meta` -- two samples may share
a run number, or omit it, so `RunNumber` alone cannot attribute a row to its file. `Role` is stored as the underlying
value of `ERole`: `0` signal, `1` background, `2` both.

`Process` is the MC-truth label (`Skimmer/Process.hxx`) and is independent of the sample's role: `0` signal, `1` hybrid,
`2` combinatorial, `3` reference (wrong-sign or mixed channel), `4` unknown. The normalization denominator for a signal
efficiency is `sum(NEvents * NInjectedPerEvent)` over the signal-side `Meta` rows -- `NEvents` comes from the `N_Events`
histogram `t2ds` fills once per event, before any candidate is built, which is why it is the only correct denominator.
