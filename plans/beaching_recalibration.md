# Beaching recalibration: move the rate scale into `τ0`

**Status:** proposed, not implemented. Current production is `τ0 = 24 h`,
`w_half = 1.5 m/s` — see [../docs/beaching.md](../docs/beaching.md).

## Problem

The rate is `τ = τ0 / s(w_onshore)` with `s(w) = w/(w + w_half)`, a saturating
ramp. `w_half` is meant to be the half-saturation point: the forcing at which
the rate reaches half its ceiling. **Production sets it far outside the range
the forcing ever takes**, so the saturation is never exercised.

Measured onshore Stokes over in-band steps (Aug 2019, `surface_stokes`; only
57.5 % of in-band steps have `w_onshore > 0` at all):

| quantile | p10 | p25 | p50 | p75 | p90 | p99 | max |
|---|---|---|---|---|---|---|---|
| `w_onshore` (m/s) | 0.0037 | 0.0109 | 0.0256 | 0.0566 | 0.0884 | 0.1527 | 0.457 |

And a hard ceiling on the field itself: over 42 M samples of `|VSD|` spanning
2019, the **global maximum anywhere in the Baltic is 0.544 m/s**. Nothing
approaches 1 m/s, let alone 1.5.

Consequences of `w_half = 1.5`:

- `s` ranges 0.0025 → 0.234, never near 1. Since `w ≪ w_half` always,
  `s ≈ w/w_half` and the model degenerates to `τ ≈ τ0·w_half / w` — a pure
  inverse law in `w` with **one** effective parameter, the product
  `τ0·w_half`. `τ0` and `w_half` are not separately identifiable.
- The high-forcing tail is unbounded: τ falls to 4.3 d at the strongest
  forcing, with no floor. The saturating form exists precisely to impose that
  floor (beyond some forcing, stranding is near-certain within a step and the
  rate should cap).
- `w_half` cannot be described as a half-saturation constant in the paper
  without being wrong. It is a rate scale wearing the wrong name.

## Proposal

Keep `w_half` **inside** the measured distribution and move the rate scale to
`τ0`. Holding median τ ≈ 60 d (the current production behaviour):

| | `τ0` | `w_half` | `s` range | τ range (p10→max) |
|---|---|---|---|---|
| current | 24 h | 1.5 m/s | 0.003 – 0.23 | 406 → 4.3 d |
| proposed | ~490 h | 0.05 m/s | 0.07 – 0.90 | 296 → 22.6 d |

`w_half = 0.05` sits between the p50 and p75 of the measured forcing, so the
ramp spans its useful range: weak forcing is strongly suppressed, strong
forcing saturates against a floor of `τ0`. Both parameters then mean what
their names say, and they become separately identifiable — `τ0` sets the
floor, `w_half` sets where the transition happens.

`τ0 ≈ 490 h ≈ 20 d` reads as "a propagule pinned against a wave-battered
shore strands in about three weeks", which is a statement that can be argued
about on its merits. `τ0 = 24 h` currently reads as a 1-day floor that the
model never gets anywhere near.

## Work

1. Re-run `024d` + `024e` at `(τ0, w_half) = (490 h, 0.05)`; confirm the
   beached fraction lands near the current 30.8 %. Tune `τ0` if not — it is a
   pure scaling of `A`, so the beached fraction moves monotonically with it.
2. Re-render `029`/`030`.
3. Sweep `τ0` (not `w_half`) for the sensitivity range, since with `w_half`
   fixed inside the data range `τ0` is the meaningful axis. `031` takes the
   swept parameter from the store filename, so it needs a `τ0` variant of the
   same partitioning.
4. Update `beaching.md`: `w_half` becomes a fixed structural choice justified
   by the forcing distribution; `τ0` becomes the reported uncertain parameter.

## Why not just leave it

The numbers are not wrong — `(24 h, 1.5)` and `(490 h, 0.05)` give the same
median timescale and similar totals. What is wrong is the *description*: a
reviewer reading `w_half = 1.5 m/s` as a half-saturation constant for a field
that never exceeds 0.54 m/s will conclude the parameterisation was never
checked against its own forcing. Cheaper to fix than to defend.

Related: [../docs/beaching.md](../docs/beaching.md) (rate model and provenance),
[../docs/wam_extrapolation.md](../docs/wam_extrapolation.md) (how `w_onshore`
is obtained).
