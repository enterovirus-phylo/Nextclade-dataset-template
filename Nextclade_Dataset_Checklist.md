# Nextclade Dataset Release Checklist

## Setup & Documentation
- [ ] Create repository in `enterovirus-phylo`
- [ ] Generate inferred root (first release only)
- [ ] Document build pipeline: commands, parameters, reproducibility notes

## Nomenclature & Configuration
- [ ] Define and validate clade nomenclature (consult experts if needed)
- [ ] Confirm clade-defining mutations
- [ ] Adjust QC thresholds for your virus
- [ ] Set `--min-seed-match` to exclude non-target Enteroviruses

## Testing
- [ ] Target virus sequences (all representatives)
- [ ] Other viruses from same species
- [ ] Different species (should fail gracefully)
- [ ] Edge cases: short sequences, high N-content, partial genomes, recombinants

## Documentation
- [ ] Write README: scope, clades, QC behavior, limitations
- [ ] Provide example sequences with expected outputs
- [ ] Update CHANGELOG with `## Unreleased` header

## Release
- [ ] Create PR on `nextclade_data`
- [ ] Test in Nextclade Web UI (tree, alignment, performance)
- [ ] Create and test playground link
- [ ] Notify stakeholders (ENPEN, collaborators)
- [ ] Merge and publish
- [ ] Monitor for issues post-release
