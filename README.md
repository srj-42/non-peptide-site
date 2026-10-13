# Chemical LR Resource — static browse site

A zero-dependency static site over the chemical half of `data/release_v2/`. Same
construction as `../site/` (no framework, no CDN, no build step), and shares its
stylesheet so the two read as siblings.

**It is a separate site, not a tab in the peptide one.** The two halves are not
comparable — chemical rows have no confidence score and no verified evidence — and
putting them in one table would invite exactly the comparison PROJECT_REPORT.md §7
argues against.

## Run it

    xdg-open index.html

Data loads through `<script src>` rather than `fetch()`, so this works from `file://`
(browsers block `fetch` on local files, but not script tags). Or:

    ./serve.sh 8080

## Rebuild after changing the data

    cd ../curation
    python3 np6_export.py         # curation.db -> data/release_v2/chemical_*.csv
    python3 np7_build_site.py     # -> site_chemical/data/

`np7_build_site.py` also regenerates `styles.css` by appending the chemical-specific
rules to the peptide site's stylesheet, so a change there propagates here on rebuild.

## Layout

    index.html             five tabs: Browse pairs, Chemicals, Statistics,
                           Methods & caveats, Downloads
    app.js                 filtering, sorting, paging, detail drawer
    styles.css             generated: ../site/styles.css + chemical additions
    data/interactions.js   13,523 pairs, loaded up front (~900 KB)
    data/chemicals.js      689 chemical identities (~230 KB)
    data/enzymes.js        15,021 enzyme links grouped by chemical (~630 KB)
    data/prov/p0..15.js    provenance, sharded; a shard loads only when you open a
                           pair in it, so the initial load stays ~1.8 MB

`np7_build_site.py` and `app.js` must agree on the shard function
(`sum of char codes mod 16`) — that is the one coupling between them. It is checked
by routing every pair id through both sides after a build.

## What it shows

Clicking a row opens a drawer with three things:

1. **Chemical identity** — HMDB and PubChem links, formula, SMILES, both class
   vocabularies. ChEBI is shown as deliberately empty; no source supplies one.
2. **How the ligand becomes measurable** — the enzyme set grouped by role
   (precursor synthesis → synthesis → efflux / influx → degradation). A chemical with
   no enzymes is called out as *not inferable from transcriptomic data*, which is the
   thing that actually determines whether a row is usable. A `?` marks a transporter
   whose direction the source never stated.
3. **Provenance** — every upstream row verbatim, as the source wrote it, with PMIDs
   linked where the source supplied any. Sources that carry per-pair citations are
   marked; the rest are labelled "no citation".

The Browse tab leads with a warning giving the two numbers that matter most: the share
of pairs resting on a single source, and the share carrying any citation at all.

## What it deliberately does not do

- No score, anywhere. `n_sources` is provenance; cross-database consensus is rejected
  as a scoring signal (PLAN.md invariant 4).
- No claim of verification. Every citation is reproduced as published and none has
  been checked.
- No MEBOCOST. That source forbids redistribution, so its rows are withheld upstream
  of this site and it is stripped from the provenance lists shown here.
