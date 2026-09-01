# Seeder + KalmanFitter

A rewrite of [**KFParticle**](https://github.com/alisw/KFParticle), where the algebra machinery is done with
[**Eigen**](https://gitlab.com/libeigen/eigen).

Deliberate numerical differences:

- `double` instead of `float`, everywhere.
- `InvertSym3` uses Eigen's closed-form 3x3 cofactor inverse, not the original's `InvertCholetsky3`. Every fit
  solves against the same matrix for a whole batch of right-hand sides at once -- 8, 15 and 8 columns for
  `VertexUpdate`, `VertexUpdateMC` and `AtProductionVertex` respectively -- so one inverse spends a single division
  where an LDLT solve spends three per column. The input is always a sum of position covariances, SPD and better
  conditioned than either term, so the cofactor form also keeps far more precision than a `float` Cholesky
  inverse. The result is forced symmetric afterwards.
- `var(S)`'s middle term is restored in `AtProductionVertex`. With `j = ds_dr` and `j3` its position part,
  `var(S) = j' C j - 2 j' G j3 + j3' V_decay j3`; the original drops the middle term, treating the decay point as
  independent of the PV-constrained state. That is exact before the update -- the transport lands on the PCA, so
  `j3` is parallel to `p` and the state's dependence on the decay point projects out -- but the update mixes in a
  K-weighted piece of the position measurement and breaks it. The dropped term vanishes again as `K -> 0` and as
  `K -> 1`, so it is a correction rather than a factor of two. It is restored because `DecayLengthErr()` feeds a
  significance that gets cut on.
- `PinToMassShell` seeds lambda with the conjugate root. The textbook smaller root of the quadratic part is
  `(E² + p² - sqrt(d)) / a`, which subtracts two nearly equal numbers in exactly the near-on-shell case that
  dominates here. Since the roots multiply to `c0/a`, the same root is `c0 / (E² + p² + sqrt(d))` -- no
  cancellation, and it degrades gracefully to the linear root `-c0/b` as `a -> 0`, so no separate fallback is
  needed.

Tricky but exact reductions:

- `PropagateCov`'s 6x6 + 6x2 split is the old 8x8 `J C J'`. Rows and columns 6-7 (`E` and `S`) really are
  transported by the identity, so only three of the four blocks move; the full triple product would spend most of
  its work multiplying by zeros.
- `corr` is never materialised. The only way one leg's state reaches the other is through the scalar `ds`, so
  `corr = d(r')/d(ds) · d(ds)/d(r_other)'` is exactly rank one. Keeping the two factors instead turns every
  `corr V corr'` into a variance of `ds` times an outer product, and every `corr X jacob'` into a row or column
  times one.
- `mDf` narrowed 7x7 -> 4x7. The original rotates all seven rows through `mJ2 · mDf · mJ1'` but only ever
  folds rows 3-6 back into `fC`; and both Jacobians are the identity outside their `(px,py,pz,E)` corner, so the
  position columns pick up only the left factor and a pin that didn't fire skips its side outright.
- `cross` reproduces the original's `D = F3·C1·F1' + F4·C2·F2'` with the same index order. `mD.transpose()`
  going into `ApplyCrossCorrection` mirrors `SetProductionVertex`'s opposite `D[k][i]` ordering -- the two
  originals genuinely disagree with each other: `AddDaughterWithEnergyFit` has `D[i][k]  K2[k][j]`,
  `SetProductionVertex` has `D[k][i]  K2[k][j]`.
- `PinToMassShell`'s 4x4 `j` is the old 7x7 restricted to its non-identity corner. The untouched position
  block comes out bit-identical, not merely equal to roundoff -- which is what keeps the later cross correction in
  `VertexUpdateMC` valid.
- `mK = mCHt  mS_inv`. `S` is symmetric (and forced so), hence so is `S^-1`, so both transposes in
  `(S^-1 CHt')'` drop out and what is left is a plain gemm.
- `VertexUpdateMC` mirrors one 4x3 block instead of re-symmetrising the whole 7x7. The one-sided
  `fC[3:7, 0:3] += mDf[:, 0:3]` fold is the only structural asymmetry either update introduces: everything
  else `fC`'s 7x7 corner picks up arrives as `M + M'` (`mDf`'s momentum-energy block, `ApplyCrossCorrection`)
  or through a `selfadjointView` (`PinDaughterToMassShell`). Mirroring that block into `fC[0:3, 3:7]` costs 12
  stores where the old full `selfadjointView<Lower>` cost a 49-element temporary and two full passes, and it
  leaves `VertexUpdateMC` on exactly the invariant `VertexUpdate` already held.

Known omissions and latent inconsistencies:

- `PinToMassShell` leaves row and column 7 (`S`) alone, matching the original. That leaves `S`'s cross-covariances
  with the rescaled momenta stale -- harmless only because they are zero until `AtProductionVertex` fills them,
  but a real inconsistency rather than a mere omission.
- `fC(7,6)`, the `S`-`E` covariance, is never filled. The original's loop stops at `fC[33]` too.
- `AtProductionVertex` assumes the vertex is independent of the candidate. If the granddaughter tracks entered the
  primary-vertex fit upstream, that correlation is of the same order as the restored `var(S)` term, and nothing at
  this level can recover it.
- A singular input to `InvertSym3` yields inf/nan, unguarded, and nothing downstream checks. An LDLT solve would
  do the same, so this is no worse -- but neither is it a check.
- Both `VertexUpdate` and `VertexUpdateMC` are kept as monolithic as possible, despite their similarities.
- `Particle::fC` is symmetric up to roundoff, not bitwise. `PropagateCov`'s `J C J'` and the vertex updates'
  `mK · mCHt'` are symmetric in exact arithmetic only: Eigen evaluates `(i,j)` and `(j,i)` as different summation orders
  over differently-rounded factors, so they can disagree in the last ulp. Nothing consumes `fC` at that precision --
  `Cov<N>()` and every named accessor read the lower triangle -- but the upper triangle is read (`PropagateCov`'s 6x2
  corner, `PinMotherMass`'s 4x4 block), so this is an invariant to state rather than one to lean on.
