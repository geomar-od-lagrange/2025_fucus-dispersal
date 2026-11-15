# Beaching statistics for bottom particles

Status: **plan, not started (2026-11-15).** Touches 002, 024d–g and the
beaching consumers; start only after the coastline-grid-fix rerun has landed
([coastline_grid_fix.md](coastline_grid_fix.md)), which is editing and
rerunning the same notebooks.

## Premise

`bottom` / `bottom_ub49` drifters feel no Stokes drift, but in water shallow
against the wave length the orbital motion reaches the bed and the
surface/bed distinction fades. Bottom drifters then strand under **the same
rate model** ([beaching.md](../docs/beaching.md)) with onshore Stokes `w_on`
replaced by its value at the bed, `w_bed = F·w_on`. Surface drifters sit at
`z = 0`: `F ≡ 1`. The double-counting objection (Daily et al.) holds a
fortiori: the bottom runs carry no Stokes at all.

## Step 0: band width, decided first

The production band `band_m = 2000` is questioned: under strong onshore Stokes
surface drifters reach the wall anyway (ε rest point 0.01–0.07 cell ≈ 10–60 m
fine, 50–350 m coarse; [wall_velocity.md](../docs/wall_velocity.md)), so a
~500 m band may be the right one, for **both** `surface_stokes` and bottom,
with the same absolute `w_c = 0.10` m/s. **Author decision; this test is its
basis**, on the existing `surface_stokes` sidecars, no 024d rebuild:

- 024e/024f members at `band_m = 750` (see below) against 2000, production
  trap/wave point, all release years: stranded %, afloat-by-age curve (029),
  stranding map and receiving-hex count (029/031). Expectation if the argument
  holds: a small change for surface. Repeat on the bottom prototype
  (Validation): there a **much larger** reduction is expected, since bottom
  drifters are not pushed onto the wall and spread across the band.
- **Tag collision:** `band_m` is in the parquet attrs but not in the member
  tag, so a band member would overwrite production. Add `_b<band_m>` to the
  tag whenever `band_m ≠ 2000` (or always, if the band changes production).
- **Raster quantisation:** `dist` is the EDT of the 500 m raster
  (`raster_dx_m`) in whole-cell units, and the reducer gates `dist·25 <
  band_m`, so `band_m = 500` selects *nothing*; 501–700 m takes the
  edge-adjacent ring, 701–1000 m adds the diagonals (707 m → 700 stored). "500 m band" = the first
  raster ring, `band_m = 750` (water within ~0–500 m of the raster land
  edge, which itself sits within ±250 m of the model coast). A finer raster
  (250 m) is a full 024d rebuild per regime (~50 min on the ~140-node
  allocation, 4× the raster arrays, sidecar size unchanged) and only sharpens
  the fine nest — the coarse model coast is 5 km anyway.
- If the band changes: beaching.md must state `w_c` as an **absolute**
  threshold (today: "the p93 of in-band forcing", a 2 km-band statistic),
  and every number in its Production/Sensitivity tables is rebuilt.

## The condition: Stokes at the bed

Runs are 2D; a bottom drifter sits on the bed, depth `h = H0` of its model
cell ([h0_semantics.md](../docs/h0_semantics.md)), fine over coarse — the H0
024d samples for its land raster. Monochromatic linear waves in depth `h`:

```
u_s(z) = ω k a² cosh(2k(z+h)) / (2 sinh²(kh))
F(kh)  = u_s(−h) / u_s(0) = 1 / cosh(2kh)        (deep water: 2·e^(−2kh))
ω² = g k tanh(kh),  ω = 2π/T                      (3 Newton steps from Eckart)
```

`F` scales the WAM surface Stokes (a spectral integral); only the decay takes
one period. `H0 ≤ 0` (tidal flat) → `kh = 0`, `F = 1`. `F > 0.5` needs
`kh < 0.66`: `h ≲ 2.5 / 5 / 8 m` at `T = 4 / 6 / 8 s`.

| option | needs | trade-off |
|---|---|---|
| **`F = sech(2kh)`** (recommended) | `h`, one period | continuous, no free knob; one period stands in for the spectrum |
| Deep-water spectral profile (Breivik et al. 2014) | `Hs`, `Tm01` | better spectral shape, biased low exactly where `kh` is small |
| `h < h_c` | `h` only | no download; ignores that long storm waves reach deeper; `h_c` a free knob |
| `kh < κ` | `h`, `T` | same data as `sech`, adds a threshold, drops the continuous decay |
| `F ≡ 1` in band | nothing | upper bound (bottom as if at surface); keep as a bracket member |

**Period.** Surface Stokes weights the spectrum by ω³, so its effective period
is shorter than the peak. The catalogue (`cmems_mod_bal_wav_my_PT1H-i`) has
`VTPK`, `VTM10`, `VTM02` (no total `VTM01`); `Tm02 ≤ Tm10 ≤ Tp`. Recommend
**`VTM02`** (closest to the Stokes-weighted period, smallest `F`:
conservative); `VTPK` as the prototype's liberal bracket.

**Download.** `002` pulls only `VSDX`/`VSDY` (2019-10-18 file: those two,
float32, uncompressed, 118 MB/day). `VTM02` adds 59 MB/day → ~230 GB for
2016-01 … 2026-09 (`baltic_highres` only; 024d reads no WAVERYS), less with
`netcdf_compression_level`; ~6 h over ten parallel year-windows. Own day file
(`stokes/baltic_highres/<year>/tm02_<YYYYMMDD>.nc`) so 002's skip-if-exists
resume works per variable. Check the NRT dataset carries `VTM02`.

**Depth caveats.** H0 is a cell mean (fine ~0.9 km, coarse ~5 km), deeper than
the true nearshore depth, so `F` is **underestimated** near the coast —
conservative. Sea level is ignored: Baltic SSH is tens of cm (surge ~1–2 m in
the west), small against the cell-mean bias, and `H0 + η` costs an hourly 2D
read of the BSH z files; set-up coincides with onshore wind, so it would lower
`F` in the strong-wave hours. A nearshore slope `h = min(H0, β·dist)` would
fix the bias but adds a free `β` — open question.

**Optional diagnostic only** (not a gate): share of bottom in-band hours with
orbital bed speed `u_b = πHs/(T sinh kh)` above `U_b = 0.049 m/s`
([bottom_friction.md](../docs/bottom_friction.md)). Needs `VHM0` (~230 GB).

## Pipeline changes

- **002**: `VTM02` day files as above.
- **024d**: `_ZARR_STEM_RE` accepts `bottom_ub\d+(p\d+)?` (today it rejects
  `bottom_ub49`; same in 024e/f/g). New parameter `particle_level`
  (`"surface"` | `"bed"`), set per regime by the job script and written as a
  sidecar attr — explicit plumbing, no regime-name parsing (the zarrs carry no
  level). New column **`kh`**: uint8 in 0.01 steps, in band only, 255 = out of
  band / no data / `kh ≥ 2.55` (`F < 0.012`, read as 0); `h` from the per-block
  H0 raster lookup, `T` sampled like `VSDX/VSDY` (nearest hour/cell, same donor
  map). `kh` rather than `F` keeps the profile form a reducer choice and the
  `kh < κ` variant free; `ΔF ≤ 0.01`. Written for every regime (one schema),
  attr `period_var`; `surface_stokes` sidecars need no rebuild for the numbers.
- **024e/f/g**: each `beaching_exponent` copy uses `w = w_on·sech(2·kh)` when
  `particle_level == "bed"`, else `w_on`. **Members unchanged** (same absolute
  `w_c`, `τ`, trap, band; `τ_calm` unscaled in band). Optional parameter
  `bed_profile` (`"sech"` | `"none"` = `F ≡ 1`, tag `_fnone`); `_b<band>` tag
  as in Step 0. Filenames carry the regime: no collisions across regimes.
- **Consumers**: 029/030/031/033 take `regime` already as a parameter; 035
  `stranded_regime` → one stranded curve per regime of `regimes_csv`; 032
  table A per regime; 034 optional stranded row `bottom` vs `bottom_ub49`.
- **Scripts**: 024d and 024efg job scripts take `regime` instead of pinning
  `surface_stokes` (submit_024d/024efg pass it); `submit_downstream_chain.sh`
  stages B/C loop over the chain regimes; `025-031_ConsumersAll_job.sh` gates
  029–031 per regime and builds the beaching-member glob per regime.
- **Docs**: beaching.md, survival_occupancy.md, bottom_friction.md,
  job_scripts.md ("beaching is `surface_stokes`-only" goes), stokes_drift.md,
  AGENTS.md bullets for 002/024d/035.

## Cost and disk

- 024d: bottom zarrs 2.8 GB/release vs 1.8 GB, same 87 200 × 5280 shape;
  ~1.5–2× the surface per-zarr cost → ~2–3 h for both regimes (~1460 zarrs)
  on the ~140-node allocation; +59 MB day cube, peak stays ≲ 12 GB of 18 GB.
- Disk: ~93 MB/sidecar + `kh` (5–15 MB) → ~160 GB for both regimes; ~230 GB
  download. 024e/f/g: seconds per zarr per member, ×2 regimes.

## Validation (`tmp_bottom_beaching/`, release 2019-10-18)

1. `VTM02` + `VTPK` for the 220-d window (~26 GB, tmp only).
2. Prototype sidecars (script copy of 024d + `kh`) for `surface_stokes`,
   `bottom`, `bottom_ub49` of `Fucus_BSH_20191018T001500_*`.
3. Report: H0, `kh`, `F` distributions in band (fine vs coarse nest); share of
   bottom in-band hours with `F > 0.5` and with `F·w_on ≥ w_c` (vs surface);
   stranded % at production and inert members for the three regimes ×
   {`Tm02`, `Tp`, `F ≡ 1`} × `band_m` {750, 2000}.
4. Sanity: new 024d on the `surface_stokes` zarr reproduces the production
   sidecar's five arrays bit for bit; new 024e/f/g reproduce that release's
   production parquets exactly (`F ≡ 1` leaves `w_on` untouched).
5. **Gate:** if bottom strong-wave hours are < ~10 % of surface ones, bottom
   stranding is near zero and mostly set by the H0 bias — numbers to the
   authors before the full build.

## Open questions for the authors

1. Band: 2000 m or the first raster ring (`band_m = 750`), for all regimes
   (Step 0)? Finer raster wanted?
2. Period for `F`: `Tm02` (recommended), `Tp`, or both as a bracket?
3. Depth: H0 only (recommended), `H0 + η`, or a nearshore-slope correction?
4. Same absolute `w_c = 0.10` m/s for bottom (recommended)?
5. Report `F ≡ 1` as the bottom upper bound? Pull `VHM0` for the `u_b` share?
6. Which bottom products: 029/031/032/035 only, or also 030, 024f/g-weighted
   025–028 and a stranded row in 034?
