# Subbasin borders on the model coastline

Authors 2026-11-17 (figure_layout_round §G): "For basins, the model coastline
is the reference." Borders must be continuous and end exactly on the BSH
HBMnoku model coastline the maps draw.

## Construction

Borders are the edges of the **nearest-HELCOM-polygon partition of model
water**:

1. Rasterise the 17 HELCOM level-2 polygons on a 200 m EPSG:3035 grid
   (holes = islands left unlabelled); drop label components < 25 km²
   (archipelago slivers that would leave enclave loops).
2. Fill every unlabelled pixel with the label of the nearest labelled pixel
   (`scipy.ndimage.distance_transform_edt`), up to 25 km reach; beyond it
   stays unlabelled (Skagerrak / North Sea, where no pair of subbasins
   meets and an unbounded fill would draw a Kattegat–Kiel Bay line across
   Jutland into the North Sea).
3. Pixel edges between two different labels, per label pair → `linemerge`,
   Douglas–Peucker at 1.5 pixels (removes the 200 m staircase).
4. To EPSG:4326, intersect with the model water polygon (the maps' fine-first
   merge of `bsh_hbmnoku_static/coastline.geojson`), so every line ends on
   the model coastline exactly and none runs over model land.

Inside the HELCOM polygons the partition is the HELCOM one (distance 0); in
the strips where model water and HELCOM disagree, it is the equidistant line,
the natural continuation of the border to the model coast. Prototype
(`tmp_borders/`): the same 24 subbasin pairs as the HELCOM shared edges,
0 dangling ends (every end on the model coast or a junction), 0 m on model
land.

Considered: snap / `set_precision` HELCOM topology (closes vertex gaps, not
the coast mismatch); extend open ends along their last direction (ambiguous
in archipelagos, can cross islands); 024a hex-edge borders (hex-stepped at
6 km).

## Where

~60 lines of raster + vector work fail the idiom test and would be copied
into 11 map notebooks, so it is precomputed once. New stage
`notebooks/006_BuildSubbasinBorders.py` (00x: static inputs only, after 004's
coastline), writing `output_root/SubbasinBorders/subbasin_borders.geojson`
(one row per subbasin pair: `subbasin_a`, `subbasin_b`, geometry). Considered
024a; rejected — it is per radius and rerunning it rewrites the key every
store hangs off. The map notebooks (025, 026, 026a, 026b, 027, 029, 030,
033, 034, 036, 037) replace their pairwise-intersection block by reading the
geojson (033 keeps the HELCOM polygons for its code labels).

## Steps

1. Plan (this file). 2. 006 notebook + job script, run. 3. Map notebooks.
4. Docs (visualisations.md backdrop, job_scripts.md, AGENTS.md stages).
5. Rerun map figures (ALL on the critical path, seasons in background),
   recopy paper PNGs.
