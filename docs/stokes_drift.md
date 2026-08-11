# Stokes drift forcing

The `surface_stokes` regime adds wave-induced Stokes drift to BSH
surface Eulerian currents. Combined field is written by
`notebooks/003_prepare_2d_fields.py` (variant `surface_stokes`); from
Parcels' perspective it's just another 2D fieldset, no kernel changes.
Wave-orbit asymmetry is comparable to the Eulerian current in storm
conditions; the contribution decays as `exp(-2 k z)`, so it is added
to the surface variant only.

## Source datasets (layered)

The Baltic high-resolution hindcast covers most of the BSH domain at
2 km / hourly; the global WAVERYS hindcast fills the German Bight
strip (BSH water west of ~9 °E) at coarser resolution.

| Property       | Baltic high-res                              | WAVERYS global                          |
|----------------|----------------------------------------------|-----------------------------------------|
| Dataset ID     | `cmems_mod_bal_wav_my_PT1H-i`                | `cmems_mod_glo_wav_my_0.2deg_PT3H-i`    |
| Resolution     | 2 km, hourly                                 | 0.2°, 3-hourly                          |
| Coverage       | 9–30 °E, 53–66 °N (Baltic + Danish Straits + Kattegat) | global                       |
| Period         | 1980 onwards                                 | 1980 onwards                            |
| Variables used | `VSDX`, `VSDY` (surface Stokes drift, m/s)   | same                                    |

`interpolate_stokes` interpolates both onto BSH U/V faces, then layers
them: Baltic high-res where defined, WAVERYS where Baltic is NaN,
zero where neither covers (e.g. far west of the BSH coarse-grid edge,
which has no Fucus releases anyway). Without the layering, the strip
between ~6.2 and 9 °E (German Bight, SH North Sea coast, Heligoland
Bight) would silently fall back to BSH-only physics under the
`surface_stokes` regime.

## Download

`notebooks/002_download_stokes.py` pulls each product day-by-day under
`output_root/stokes/<product>/<year>/stokes_<YYYYMMDD>.nc`, resumable:

```bash
pixi run python notebooks/002_download_stokes.py \
    --output-root /path/to/outputs --year 2020
```

## Regrid onto BSH C-grid

CMEMS delivers Stokes on a regular A-grid (cell centres); BSH lives
on a NEMO-style C-grid (U east face, V south face). `interpolate_stokes`
flips lat to descending, linearly interpolates `VSDX`/`VSDY` onto the
matching U/V face coordinates, and nearest-neighbours in time onto the
BSH 15-min step (sub-hourly refinement of the wave field would
over-resolve — Baltic is hourly, WAVERYS is 3-hourly).

The combined field is summed in BSH convention before the edge mask
and V-roll, so masking after summation ensures no flux across the BSH
coastline even where wave-model and BSH disagreed about land. Hygiene
chain in [2d_field_extraction.md](2d_field_extraction.md). Considered
in-kernel runtime sum; rejected for I/O cost and harder face-alignment
validation.

Stokes is also shut off per-timestep where BSH says the face is
blocked (`u_surf == 0` or NaN) so tidal flats correctly receive Stokes
when wet and zero Stokes when dry, and the Stokes contribution doesn't
push particles through no-slip walls.

## Open concern: spread bridges thin land barriers

The 3×3 rolling-mean spread (N=5 iterations, in
`_spread_into_nan` / `_interp_stokes_to_bsh`) bridges land barriers
narrower than the spread reach. With Baltic high-res at 2 km grid and
N=5, the spread crosses ~10 km of wave-model land — wide enough to
push open-ocean Stokes across thin barriers like the Curonian Spit
(~1–3 km wide) into the sheltered Curonian Lagoon, where actual fetch
resets and real Stokes is small.

Real wave physics: Stokes drift depends on wind-driven wave growth
along an unobstructed fetch. Thin land barriers reset fetch on their
leeward side; the sheltered water body sees small Stokes regardless
of conditions on the open side.

The Baltic high-res wave model encodes this physics — it returns NaN
inside such enclosed lagoons. Our current spread overrides that
silence and propagates open-ocean values inward.

Possible mitigations (none implemented):

- **Reduce N to 1** so the spread only catches the BSH-vs-wave-model
  coastline 1-cell drift, minimising barrier crossing. Costs ~12 % of
  near-coast BSH cells losing Stokes entirely (those farther than one
  wave-model cell from the wave-model coast).
- **Treat Baltic-bbox interior NaN as zero** in the layering step
  (rather than falling back to WAVERYS, which is too coarse to resolve
  the barriers anyway). Pairs naturally with reduced N.
- **Geodesic spread** that respects the wave-model land/water
  topology, not crossing through cells that were originally land in
  the wave-model file. More complex; requires a connected-components
  analysis.

Pending further thought; kept at N=5 with the per-timestep face mask
as the production setting for now.

## Cross-references

- [2d_field_extraction.md](2d_field_extraction.md) — preprocessor;
  Stokes is one of three variants.
- [seeding.md](seeding.md) — how the sweep targets `surface_stokes`.
- [wam_extrapolation.md](wam_extrapolation.md) — how the beaching
  diagnostic re-reads this field at the coast, where the blocked-face
  mask above suppresses it.
