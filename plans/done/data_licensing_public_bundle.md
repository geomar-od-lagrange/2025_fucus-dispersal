# Data licensing for a public companion repo

> Implementation landed as the curated input bundle in the data twin
> (`<https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>`)
> plus the per-dataset attributions in `ATTRIBUTION.md` at the repo
> root. Original licensing audit retained as historical record.

## Goal

Bundle a minimal reproducibility set — HELCOM subbasin shapes, the Fucus
location shapefile, and a tiny BSH operational-model subset — into a
*public* git repo alongside the code. That requires confirming, per
dataset, that we have the right to redistribute and what attribution is
required.

This plan records what was checked, where the evidence lives, and what
still needs action. Written so another session can pick it up without
re-doing the research.

## Datasets in scope

1. `data/HELCOM_subbasins_2022_level2/` — HELCOM subbasin polygons
2. `data/Fucus_location_shp/` — `REDLIST_SIS_Macrophytes.shp`
3. `min_data/bsh_operationalmodel_data/` — BSH HBMnoku outputs
   (currents, temperature, salinity, static grid files) for a short
   demo window
4. `min_data/stokes/baltic_stokes_20200101.nc` — Stokes drift sample
   from CMEMS Baltic Sea Wave Hindcast
   (`cmems_mod_bal_wav_my_PT1H-i`, product
   `BALTICSEA_MULTIYEAR_WAV_003_015`), pulled by
   `min_data/download_stokes_sample.py`

Out of scope here: DWD ICON forcing (not shipped directly).

## Findings

### HELCOM subbasins — OK to bundle with attribution

Source of terms: `data/HELCOM_subbasins_2022_level2/metadata.xml`,
`gmd:resourceConstraints` block.

- Use constraint (verbatim):
  *"Data can be used freely given that the source (HELCOM) is cited."*
- Access constraint: INSPIRE `noLimitations`
  — *"Access constraints: No limitations on public access."*
- Catalog record:
  `https://metadata.helcom.fi/geonetwork/srv/eng/catalog.search#/metadata/d4b6296c-fd19-462c-94d2-4c81b9313d77`

**Action:** redistribute as-is; cite HELCOM + catalog UUID in
`ATTRIBUTION.md`.

### BSH HBMnoku — OK to bundle under dl-de/by-2-0

The NetCDF globals in `min_data/bsh_operationalmodel_data/*/*.nc` carry
`institution = Bundesamt fuer Seeschifffahrt und Hydrographie` and
`references = http://www.bsh.de` but **no** `license` attribute. The
licence is not in the parent Atom index
(`https://gdi.bsh.de/de/feed/Data-of-the-Operational-Models.xml`); it
is in the per-dataset child feeds.

Verified in the child feeds for both the 5 km series and the 900 m
2016 series. Identical `<rights>` element:

> *"Zugriffsbeschränkung: unclassified --- Nutzungsbedingungen: Es
> gelten keine Zugriffsbeschränkungen --- Lizenz: Dieser Datensatz
> kann gemäß der 'Datenlizenz Deutschland – Namensnennung – Version
> 2.0' (www.govdata.de/dl-de/by-2-0) genutzt werden."*

- Licence: **Datenlizenz Deutschland – Namensnennung 2.0**
  (`dl-de/by-2-0`) — <https://www.govdata.de/dl-de/by-2-0>
- Required source string: `Quelle: Bundesamt für Seeschifffahrt und Hydrographie`
- Access: no restrictions (`unclassified`)
- Contact for questions: `opmod@bsh.de`
- Fact sheet (no terms, just model docs):
  `https://gdi.bsh.de/de/data/..._FactSheet_HBMnoku_english.pdf`
- GDI-DE registry entry: <https://registry.gdi-de.org/id/de.bund.bsh>

Concretely checked feeds:

- `https://gdi.bsh.de/de/feed/Modelled-current-forecast-of-the-operational-circulation-model-of-BSH-in-the-North-and-Baltic-Sea-horizontal-resolution-5-km-series.xml`
- `https://gdi.bsh.de/de/feed/Modelled-currents-of-the-operational-circulation-model-of-BSH-in-the-German-Bight-and-the-western-Baltic-Sea-horizontal-resolution-ca-900-m-2016-series.xml`

`dl-de/by-2-0` permits copying, redistribution, adaptation, and
commercial use; the only requirement is attribution (name source +
link licence + note any modifications).

**How the bundle uses this data:** we ship individual `.nc` files
verbatim — no modification of file contents, no reprojection, no
variable stripping. The only "reduction" is that the bundle includes
only a small subset of the upstream years/files. File-selection of
unmodified originals is not a modification of the data itself, so no
"modified" notice is required beyond naming the subset extent for the
reader's benefit.

**Action:** redistribute the included `.nc` files as-is; include the
attribution line plus a plain-English note of which years/files are
present so users understand the scope.

### CMEMS Baltic Wave Hindcast / Stokes drift — CLEARED

Pulled by `min_data/download_stokes_sample.py` via the
`copernicusmarine` Python client:

```
dataset_id = "cmems_mod_bal_wav_my_PT1H-i"
variables  = ["VSDX", "VSDY"]
window     = 2020-01-01T00:00 .. 2020-01-02T00:00
output     = min_data/stokes/baltic_stokes_20200101.nc
```

NetCDF globals confirm:

- `cmems_product_id = BALTICSEA_MULTIYEAR_WAV_003_015`
- `institution = Baltic MFC, PU Finnish Meteorological Institute`
- `source = FMI-WAM_CMEMS`

**Product:** Baltic Sea Wave Hindcast — BALTICSEA_MULTIYEAR_WAV_003_015
**DOI:** <https://doi.org/10.48670/moi-00014>
**Licence:** Copernicus Marine Service Commitments and Licence —
<https://marine.copernicus.eu/user-corner/service-commitments-and-licence>

Permission for redistribution is explicit (Section 2.2(c)): licensees
may *"redistribute, disseminate any Copernicus Marine Service Product
in their original form via any media"*. No non-commercial clause,
no share-alike, no requirement that recipients register.

Required attribution (Section 2.4):

- For unchanged originals: `E.U. Copernicus Marine Service Information; <DOI link>`
- For derivatives: `Generated using E.U. Copernicus Marine Service Information; <DOI link>`
- For publications: `This study has been conducted using E.U. Copernicus Marine Service Information; <DOIs>`

Licensee obligations: *maintain records of use* and *propagate the
licence in all descending distributions* — i.e. the bundle's
`ATTRIBUTION.md` must carry the CMEMS licence forward (not just cite
the source).

**Note on form:** the sample was generated with
`copernicusmarine.subset`, which returns a subset (1 day, variables
`VSDX` and `VSDY` only). That means the shipped `.nc` is technically
a *derivative* of the upstream product, not an original copy. Use the
"Generated using …" attribution wording. CMEMS allows derivatives too,
so this is not a blocker — only an attribution-wording choice.

**Action:** redistribute the sample; use the "Generated using"
attribution form; cite DOI 10.48670/moi-00014; carry the CMEMS licence
forward in `ATTRIBUTION.md`.

### Fucus / REDLIST_SIS_Macrophytes — CLEARED

Resolution check done: the IUCN concern does not apply. IUCN is named
in the abstract only as the source of the *Red List assessment
criteria* (methodology), not as a provider of spatial data. The
`gmd:lineage` element in `data/Fucus_location_shp/metadata.xml`
enumerates the actual data sources, and IUCN is not among them.

Lineage statement (verbatim):

> *"The records of species compiled from the Danish national database
> for marine data (MADS), the German database for macrophyte
> occurrences (MARIDATA), the database of Swedish Species Information
> Centre, Botanical Museum Lund (LD), and Uppsala Museum of Evolution
> Herbarium (UPS). For the Swedish coastline the continuous
> distribution area is mainly based on expert view. (Downloaded from
> MADS and from Finnish Environment Institute on 18 June 2012)"*

pointOfContact is `HELCOM Secretariat`. The `SIS` in the filename is
not IUCN Species Information Service — the catalog record is served
via `https://maps.helcom.fi/website/MADS/...`, making it a MADS-side
naming convention.

Catalog record (Fucus vesiculosus, LC):
`https://metadata.helcom.fi/geonetwork/srv/eng/catalog.search#/metadata/5848c347-dd45-4135-bbb0-228be9ddeffb`

Use constraint (verbatim):
> *"Data can be freely used when cited the original data is from MADS
> and the Finnish Environment Institute."*

Access: `noLimitations` — *"No limitations on public access."*

**Note on scope:** the catalog UUID is for *Fucus vesiculosus*, but
the local shapefile is named `REDLIST_SIS_Macrophytes` and the
abstract says *"Distribution of the species can be found in
corresponding name column"*. Treat it as a multi-species macrophyte
bundle. Before publishing, inspect the `.dbf` columns and list the
species actually covered in the attribution block.

**Action:** redistribute with the attribution block below; inspect
`.dbf` columns to enumerate species.

## Implementation checklist

All three datasets are cleared for redistribution with attribution.
Remaining work:

- [ ] Create `ATTRIBUTION.md` at the repo root (or next to `LICENSE`)
      with one block per redistributed dataset. Each block: source
      name, link to catalog/licence, required attribution string,
      note about any modification.
- [ ] Add a sidecar `LICENSE-data` or equivalent mentioning that
      *code* is MIT (existing `LICENSE`) while *data* is under the
      per-dataset terms referenced in `ATTRIBUTION.md`.
- [ ] For BSH: files are shipped verbatim, no modification notice
      required; just list which years/files are present so users
      understand the bundle's scope.
- [ ] For CMEMS: use the "Generated using …" attribution wording
      (the Stokes sample is a subset, hence a derivative); include
      the DOI link and carry the CMEMS licence forward in
      `ATTRIBUTION.md`.
- [ ] If the public repo is separate from this one, decide whether
      HELCOM subbasins lives there or stays here (currently gitignored
      per `aaff65b gitignore: skip data/HELCOM_subbasins_2022_level2`).
- [ ] Consider removing `data/HELCOM_subbasins_2022_level2` from
      `.gitignore` if the plan is to track it *here*, or leave it
      ignored and track it *only* in the public bundle repo.

## Attribution block drafts (for copy-paste when implementing)

HELCOM subbasins:

```
HELCOM subbasins 2022 level 2
  Source: HELCOM (Baltic Marine Environment Protection Commission)
  Catalog: https://metadata.helcom.fi/geonetwork/srv/eng/catalog.search#/metadata/d4b6296c-fd19-462c-94d2-4c81b9313d77
  Terms: "Data can be used freely given that the source (HELCOM) is cited."
```

BSH HBMnoku subset:

```
BSH operational circulation model (HBMnoku) — minimal demo subset
  Source: Bundesamt für Seeschifffahrt und Hydrographie (BSH)
  Licence: Datenlizenz Deutschland – Namensnennung 2.0
           https://www.govdata.de/dl-de/by-2-0
  Attribution: Quelle: Bundesamt für Seeschifffahrt und Hydrographie
  Modifications: none (individual .nc files included verbatim;
                 a subset of upstream years/files is present)
  Scope included: <fill in — list of years/grids/file types>
  Contact: opmod@bsh.de
```

CMEMS Stokes drift sample:

```
CMEMS Baltic Sea Wave Hindcast — Stokes drift sample (VSDX, VSDY)
  Product: BALTICSEA_MULTIYEAR_WAV_003_015
  Dataset: cmems_mod_bal_wav_my_PT1H-i
  DOI:     https://doi.org/10.48670/moi-00014
  Producer: Baltic MFC, Finnish Meteorological Institute (FMI-WAM)
  Licence: Copernicus Marine Service Commitments and Licence
           https://marine.copernicus.eu/user-corner/service-commitments-and-licence
  Attribution: "Generated using E.U. Copernicus Marine Service
                Information; https://doi.org/10.48670/moi-00014"
  Scope included: 2020-01-01T00:00 .. 2020-01-02T00:00, variables
                  VSDX and VSDY only (subset produced via
                  copernicusmarine.subset)
  Notes: derivative (subset) of original product; CMEMS licence
         must be propagated to downstream users
```

Fucus / HELCOM Red List macrophytes:

```
HELCOM Red List — macrophyte distributions (incl. Fucus vesiculosus)
  Source: MADS and the Finnish Environment Institute (SYKE), via HELCOM
  Methodology: HELCOM Red List (2013), applying IUCN Red List criteria
  Compilation sources: MADS (Denmark), MARIDATA (Germany),
    Swedish Species Information Centre, Botanical Museum Lund (LD),
    Uppsala Museum of Evolution Herbarium (UPS); Swedish coastline
    distribution partially based on expert view
  Catalog: https://metadata.helcom.fi/geonetwork/srv/eng/catalog.search#/metadata/5848c347-dd45-4135-bbb0-228be9ddeffb
  Species covered: <fill in from .dbf columns>
  Terms: "Data can be freely used when cited the original data is from
          MADS and the Finnish Environment Institute."
```

## Open questions

- Where should the public bundle live — a new repo under
  `geomar-od-lagrange/` or a release asset of this one? (Licence
  obligations are the same either way; only discoverability differs.)
- Does BSH want notification when we publish the subset? The terms do
  not require it, but a courtesy email to `opmod@bsh.de` is cheap
  insurance and often yields a preferred citation form.
