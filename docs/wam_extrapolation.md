# Extrapolating the WAM wave field onto the BSH domain

The beaching rate is driven by onshore Stokes drift sampled from the CMEMS
Baltic wave hindcast (FMI-WAM, `BALTICSEA_MULTIYEAR_WAV_003_015`, 1 nmi
≈ 1.6 × 1.9 km, hourly). **WAM's water mask is not BSH's**, so a bare
nearest-cell lookup leaves much of the BSH coastline with no wave data — and
in [beaching.md](beaching.md) no wave data means *rate zero*, a coastline
that cannot strand at any parameter. This doc covers how the field is
extended to full BSH coverage.

## Scale of the gap

WAM masks its own coastal cells, and its bbox starts at 9.01 °E, excluding
the German Bight — while the beaching band is the 2 km strip hugging the
shore. So the overlap is worst exactly where it matters: **78.7 % of in-band
samples need filling**. Left alone that is a structural bias, not a parameter
choice, and it falls hardest on sheltered fjords and archipelago, i.e. prime
*Fucus* habitat.

## Method

The fill is a breadth-first expansion of the **donor index** over the static
masks — one 4-neighbour dilation per round, front advancing ~1.6 km, BSH
land blocking it. Propagating indices rather than values keeps it a one-off
precompute; per-hour sampling stays a single gather. Three details carry
weight:

- **Geodesic, not Euclidean.** A nearest-wet lookup by straight-line
  distance can draw from across a headland or outside a fjord mouth,
  importing an open-water wave climate into sheltered water that WAM
  excluded *because* it is not open water. Geodesic donors instead share
  the particle's water body. **Measured effect: none.** Rebuilding the full
  sweep on the geodesic fill reproduced the Euclidean numbers to within
  0.15 points on every member (e.g. 60.1 % → 60.2 % at `w_tau = 0.4`), and
  fill distance barely moved (2.74 → 2.67 km). So Euclidean nearest-wet
  rarely crossed land in a way that mattered. Geodesic is kept because it
  cannot do so *by construction* and because it distinguishes unreachable
  water (capped and counted) from reachable, which the Euclidean version
  could not — not because it changes the answer.
- **4-neighbour, not 8.** A 3×3 dilation squeezes between diagonally
  touching land cells — the same thin-barrier bridging the `surface_stokes`
  N=5 Stokes spread is faulted for.
- **Capped and counted.** Cells not reached within `stokes_fill_max_cells`
  rounds keep `w_onshore = 0` and are reported (currently 0.09 % of
  samples), rather than silently importing from an absurd distance as the
  Euclidean version did.

## Land versus ice

The mask that makes this safe is that **WAM NaN means land *or* ice** — the
wet-cell count varies hour to hour, ~4.4 % of the grid seasonally
ice-blanked in the Bothnian Bay and Gulf of Finland. Only *static* land is
filled, identified by an `ever_wet` mask (finite in any hour of a monthly
sample spanning the seasonal cycle). **Ice keeps `w_onshore = 0`**, the
physical answer: ice suppresses waves, so no-beaching-under-ice falls out
for free, consistent with the currents side where Fucus rides the upper-cell
velocity.

## Masks and known bias

The BSH water mask on the WAM grid uses a **3×3 footprint** rule, not centre
sampling: WAM cells are ~1.6 km and coastal ones are part land, so a centre
landing on BSH land would exclude cells particles legitimately occupy, which
would then never receive a donor and be silently forced to zero.

Known bias: filled values are read 1–3 km offshore, where Stokes drift is
stronger than at the shoreline. Zero was a large negative bias; this is a
smaller positive one. A depth taper is *not* the obvious remedy — Baltic
wind seas are short-period (Tp 3–6 s), so depth-limited breaking is confined
to a strip far narrower than the 500 m raster or the 1.6 km WAM cell.

## Cross-references

- [beaching.md](beaching.md) — the rate model this feeds.
- [stokes_drift.md](stokes_drift.md) — the wave field itself and the
  blocked-face mask this reverses at the coast.
