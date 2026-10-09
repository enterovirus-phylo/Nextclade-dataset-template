# Nextclade dataset release checklist

Work through this before publishing a dataset. The testing section is the part
most easily skipped and most worth doing: a dataset that silently accepts
sequences from a different virus is worse than no dataset.

## Setup and reproducibility

- [ ] Workflow lives in its own repository, with the template's history or a note on its origin
- [ ] Ancestral root generated and `resources/inferred-root.fasta` committed (first release only)
- [ ] `README` records the build command, the parameters used, and anything a rerun would need
- [ ] The build runs from a clean checkout — try `snakemake --cores 1 clean` followed by a full run

## Nomenclature and configuration

- [ ] Clade nomenclature defined, and checked against the published literature or with people who work on this virus
- [ ] Clade-defining mutations in `resources/clades.tsv` verified against the tree
- [ ] QC thresholds in `dataset/pathogen.json` adjusted for your genome length and expected diversity
- [ ] `minSeedCover` tuned so that sequences from related viruses fail to align rather than aligning badly
- [ ] `defaultCds` and `cdsOrderPreference` match your annotation
- [ ] Alignment parameters checked on your most divergent real sequences

## Testing

Run each of these through the CLI against the built `dataset.zip` and read the
output, rather than only checking that the command exits cleanly:

```bash
nextclade3 run --input-dataset dataset.zip --output-all test_out/ SEQUENCES.fasta
```

- [ ] **Representative sequences** of every clade — each is assigned the clade you expect
- [ ] **Related viruses** from the same species — assigned sensibly, or rejected
- [ ] **A different virus entirely** — fails to align, rather than being placed somewhere arbitrary
- [ ] **Short sequences and partial genomes** — placed or rejected as intended, with coverage reported honestly
- [ ] **High-N sequences** — flagged by the missing-data QC rule
- [ ] **Known recombinants**, if your virus recombines — behaviour understood and documented
- [ ] **The example sequences** shipped in the dataset — all pass, since they are the first thing users run

## Documentation

- [ ] `dataset/README.md` states the scope, the clade system, what QC does, and the known limitations
- [ ] Example sequences chosen to span the diversity, with their expected output described
- [ ] `dataset/CHANGELOG.md` updated under an `## Unreleased` heading
- [ ] Reference and its coordinate system named explicitly
- [ ] Maintainers and a contact route given

## Release

- [ ] Tested in Nextclade Web against a locally served `out-dataset/` — tree renders, alignment looks right, performance acceptable
- [ ] Playground link created and checked
- [ ] Pull request opened on [`nextclade_data`](https://github.com/nextstrain/nextclade_data)
- [ ] Collaborators and stakeholders notified
- [ ] Merged and published (done by Nextclade Team)
- [ ] Watch for issue reports over the following weeks
