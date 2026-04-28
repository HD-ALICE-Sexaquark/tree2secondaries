tree2secondaries -- Finder
==========================

## File layout

* `Injected` -- (`TTree`) One entry per input event.
* `FoundChannel{CH}` -- (`TTree`) One entry per input event.
* `N_Events` -- (`TH1D`) Single-bin event counter, filled once per `ProcessEvent()`.
* `CutFlow_{}` -- (`TH1D`) 20-bin cut-flow histograms, one per species considered in the chosen channel.

## TTree `Injected`

All branches below are scalars.

* `RunNumber` -- (`unsigned int`)
* `DirNumber` -- (`unsigned int`)
* `EventNumber` -- (`unsigned int`)
* `ReactionID` -- (`unsigned int`)
* `SV_X` -- (`float`)
* `SV_Y` -- (`float`)
* `SV_Z` -- (`float`)
* `Sexa_Px` -- (`float`)
* `Sexa_Py` -- (`float`)
* `Sexa_Pz` -- (`float`)
* `Sexa_E` -- (`float`)
* `Nucleon_Px` -- (`float`)
* `Nucleon_Py` -- (`float`)
* `Nucleon_Pz` -- (`float`)
* `Nucleon_E` -- (`float`)

## TTree `FoundChannel{CH}`

All branches below are scalars.

Event properties.

* `RunNumber` -- (`unsigned int`)
* `DirNumber` -- (`unsigned int`)
* `DirNumberB` -- (`unsigned int`) **(RD only)**
* `EventNumber` -- (`unsigned int`)
* `Centrality` -- (`float`)
* `MagneticField` -- (`float`)
* `PV_X` -- (`float`)
* `PV_Y` -- (`float`)
* `PV_Z` -- (`float`)
* `MC_PV_X` -- (`float`) **(MC only)**
* `MC_PV_Y` -- (`float`) **(MC only)**
* `MC_PV_Z` -- (`float`) **(MC only)**

Reconstructed properties of anti-sexaquark candidate.

* `SV_X` -- (`float`)
* `SV_Y` -- (`float`)
* `SV_Z` -- (`float`)
* `Sexa_Px` -- (`float`)
* `Sexa_Py` -- (`float`)
* `Sexa_Pz` -- (`float`)
* `Sexa_E` -- (`float`)
* `Sexa_Chi2NDF` -- (`float`)
* `Sexa_E_MinusNucleon` -- (`float`)
* `Sexa_IsAntiChannel` -- (`bool`)

**(MC only)** MC properties of linked injected anti-sexaquark. When no matching MC anti-sexaquark is found, values are
dummy (`-999.` for `float`, `-1` for `int`, `false` for `bool`).

* `MC_Before_Px` -- (`float`)
* `MC_Before_Py` -- (`float`)
* `MC_Before_Pz` -- (`float`)
* `MC_Before_E` -- (`float`)
* `MC_After_Px` -- (`float`)
* `MC_After_Py` -- (`float`)
* `MC_After_Pz` -- (`float`)
* `MC_After_E` -- (`float`)
* `MC_Nucleon_Px` -- (`float`)
* `MC_Nucleon_Py` -- (`float`)
* `MC_Nucleon_Pz` -- (`float`)
* `MC_Nucleon_E` -- (`float`)
* `MC_SV_X` -- (`float`)
* `MC_SV_Y` -- (`float`)
* `MC_SV_Z` -- (`float`)
* `MC_ReactionID` -- (`int`)
* `MC_IsSignal` -- (`bool`)
* `MC_IsHybrid` -- (`bool`)

Properties of V0 product and their negative/positive daughters.

* `{V0}_Decay_X` -- (`float`)
* `{V0}_Decay_Y` -- (`float`)
* `{V0}_Decay_Z` -- (`float`)
* `{V0}_Px` -- (`float`)
* `{V0}_Py` -- (`float`)
* `{V0}_Pz` -- (`float`)
* `{V0}_E` -- (`float`)
* `{V0}_Chi2NDF` -- (`float`)
* `{V0}_PCAwrtSV_X` -- (`float`)
* `{V0}_PCAwrtSV_Y` -- (`float`)
* `{V0}_PCAwrtSV_Z` -- (`float`)
* `{V0}_PCAwrtSV_Px` -- (`float`)
* `{V0}_PCAwrtSV_Py` -- (`float`)
* `{V0}_PCAwrtSV_Pz` -- (`float`)
* `{V0}_{Neg/Pos}_EsdEntry` -- (`unsigned int`)
* `{V0}_{Neg/Pos}_Charge` -- (`int`)
* `{V0}_{Neg/Pos}_X` -- (`float`)
* `{V0}_{Neg/Pos}_Y` -- (`float`)
* `{V0}_{Neg/Pos}_Z` -- (`float`)
* `{V0}_{Neg/Pos}_Px` -- (`float`)
* `{V0}_{Neg/Pos}_Py` -- (`float`)
* `{V0}_{Neg/Pos}_Pz` -- (`float`)
* `{V0}_{Neg/Pos}_PreDCAxy` -- (`float`)
* `{V0}_{Neg/Pos}_PreDCAz` -- (`float`)
* `{V0}_{Neg/Pos}_TPC_Signal` -- (`float`)
* `{V0}_{Neg/Pos}_NSigmaPion` -- (`float`)
* `{V0}_{Neg/Pos}_NSigmaKaon` -- (`float`)
* `{V0}_{Neg/Pos}_NSigmaProton` -- (`float`)
* `{V0}_{Neg/Pos}_PCAwrtV0_X` -- (`float`)
* `{V0}_{Neg/Pos}_PCAwrtV0_Y` -- (`float`)
* `{V0}_{Neg/Pos}_PCAwrtV0_Z` -- (`float`)
* `{V0}_{Neg/Pos}_PCAwrtV0_Px` -- (`float`)
* `{V0}_{Neg/Pos}_PCAwrtV0_Py` -- (`float`)
* `{V0}_{Neg/Pos}_PCAwrtV0_Pz` -- (`float`)

**(MC only)** MC properties of V0 product and their negative/positive daughters. When no matching MC particle is found,
values are dummy (`-999.` for `float`, `-1` for `int`, `false` for `bool`).

* `MC_{V0}_McEntry` -- (`int`)
* `MC_{V0}_PdgCode` -- (`int`)
* `MC_{V0}_Origin_X` -- (`float`)
* `MC_{V0}_Origin_Y` -- (`float`)
* `MC_{V0}_Origin_Z` -- (`float`)
* `MC_{V0}_Decay_X` -- (`float`)
* `MC_{V0}_Decay_Y` -- (`float`)
* `MC_{V0}_Decay_Z` -- (`float`)
* `MC_{V0}_Px` -- (`float`)
* `MC_{V0}_Py` -- (`float`)
* `MC_{V0}_Pz` -- (`float`)
* `MC_{V0}_E` -- (`float`)
* `MC_{V0}_IsTrue` -- (`bool`)
* `MC_{V0}_IsSecondary` -- (`bool`)
* `MC_{V0}_IsSignal` -- (`bool`)
* `MC_{V0}_IsHybrid` -- (`bool`)
* `MC_{V0}_ReactionID` -- (`int`)
* `MC_{V0}_Mother_McEntry` -- (`int`)
* `MC_{V0}_Mother_PdgCode` -- (`int`)
* `MC_{V0}_{Neg/Pos}_McEntry` -- (`int`)
* `MC_{V0}_{Neg/Pos}_PdgCode` -- (`int`)
* `MC_{V0}_{Neg/Pos}_Px` -- (`float`)
* `MC_{V0}_{Neg/Pos}_Py` -- (`float`)
* `MC_{V0}_{Neg/Pos}_Pz` -- (`float`)
* `MC_{V0}_{Neg/Pos}_E` -- (`float`)
* `MC_{V0}_{Neg/Pos}_IsTrue` -- (`bool`)
* `MC_{V0}_{Neg/Pos}_IsSecondary` -- (`bool`)
* `MC_{V0}_{Neg/Pos}_IsSignal` -- (`bool`)
* `MC_{V0}_{Neg/Pos}_ReactionID` -- (`int`)

Reconstructed properties of stable charged particle product.

* `{CT}_EsdEntry` -- (`unsigned int`)
* `{CT}_Charge` -- (`int`)
* `{CT}_X` -- (`float`)
* `{CT}_Y` -- (`float`)
* `{CT}_Z` -- (`float`)
* `{CT}_Px` -- (`float`)
* `{CT}_Py` -- (`float`)
* `{CT}_Pz` -- (`float`)
* `{CT}_PreDCAxy` -- (`float`)
* `{CT}_PreDCAz` -- (`float`)
* `{CT}_TPC_Signal` -- (`float`)
* `{CT}_NSigmaPion` -- (`float`)
* `{CT}_NSigmaKaon` -- (`float`)
* `{CT}_NSigmaProton` -- (`float`)
* `{CT}_PCAwrtSV_X` -- (`float`)
* `{CT}_PCAwrtSV_Y` -- (`float`)
* `{CT}_PCAwrtSV_Z` -- (`float`)
* `{CT}_PCAwrtSV_Px` -- (`float`)
* `{CT}_PCAwrtSV_Py` -- (`float`)
* `{CT}_PCAwrtSV_Pz` -- (`float`)

**(MC only)** MC properties of stable charged particle product. When no matching MC particle is found, values are dummy
(`-999.` for `float`, `-1` for `int`, `false` for `bool`).

* `MC_{CT}_McEntry` -- (`int`)
* `MC_{CT}_PdgCode` -- (`int`)
* `MC_{CT}_Origin_X` -- (`float`)
* `MC_{CT}_Origin_Y` -- (`float`)
* `MC_{CT}_Origin_Z` -- (`float`)
* `MC_{CT}_Px` -- (`float`)
* `MC_{CT}_Py` -- (`float`)
* `MC_{CT}_Pz` -- (`float`)
* `MC_{CT}_E` -- (`float`)
* `MC_{CT}_IsTrue` -- (`bool`)
* `MC_{CT}_IsSecondary` -- (`bool`)
* `MC_{CT}_IsSignal` -- (`bool`)
* `MC_{CT}_ReactionID` -- (`int`)
* `MC_{CT}_Mother_McEntry` -- (`int`)
* `MC_{CT}_Mother_PdgCode` -- (`int`)
* `MC_{CT}_GM_McEntry` -- (`int`)
* `MC_{CT}_GM_PdgCode` -- (`int`)
