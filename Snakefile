# =============================================================================
# 1. Settings you need to change for your virus
# =============================================================================
# Everything in this block is a placeholder. The workflow refuses to start
# until the values marked "CHANGE ME" below have been filled in, so you get a
# readable error instead of a confusing crash halfway through a long run.

REFERENCE_ACCESSION = "<accession>"     # CHANGE ME: GenBank accession of your reference, e.g. "AY426531.1"

# CDS names that augur ancestral should translate. These must match the
# product/gene names in resources/reference.gbk *and* the Name attributes in
# dataset/genome_annotation.gff3, character for character.
# The default list is the picornavirus polyprotein; replace it with your own.
GENES = ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]

# =============================================================================
# 2. Filtering and subsampling
# =============================================================================

ALLOWED_DIVERGENCE = "3000"     # Sequences with more substitutions than this are treated as outliers
MIN_DATE           = "1960-01-01"   # Drop anything collected before this date ("" disables)
MIN_LENGTH         = "6000"     # Drop sequences shorter than this, in bp ("" disables)
MAX_SEQS           = "10000"    # Upper bound on tree size after subsampling
ROOTING            = "ancestral_sequence"  # "ancestral_sequence", "mid_point", "outgroup" or an accession
ID_FIELD           = "accession"  # Metadata column holding the sequence IDs: "accession" or "strain"

# =============================================================================
# 3. Which parts of the workflow to run
# =============================================================================

# Download sequences and metadata from NCBI via the ingest/ sub-workflow.
# Set to False if you want to drop your own data/sequences.fasta and
# data/metadata.tsv in place instead.
FETCH_SEQUENCES = True

# Root the reference tree on a precomputed ancestral sequence
# (resources/inferred-root.fasta) rather than on the reference or midpoint.
# See inferred-root/README.md for the background.
STATIC_ANCESTRAL_INFERRENCE = True

# Recompute resources/inferred-root.fasta from the current sequence set.
# This is slow and overwrites the existing file, so it is off by default.
# You must turn it on for the very first build, because the template does not
# ship an inferred root.
INFERRENCE_RERUN = False

# =============================================================================
# 4. Paths
# =============================================================================
# You rarely need to touch these. If you do, keep them in sync with the copies
# in inferred-root/snakefile.

SEQUENCES = "data/sequences.fasta"      # Input nucleotide sequences, FASTA
METADATA  = "data/metadata.tsv"         # Input metadata, TSV

# When the static ancestral root is in use, the ancestor is appended to the
# sequence and metadata files and the rest of the workflow reads these copies.
INFERRED_SEQ_PATH  = "results/sequences_with_ancestral.fasta" if STATIC_ANCESTRAL_INFERRENCE else SEQUENCES
INFERRED_META_PATH = "results/metadata_with_ancestral.tsv" if STATIC_ANCESTRAL_INFERRENCE else "results/metadata.tsv"

# Files that make up the Nextclade dataset. reference.fasta and
# genome_annotation.gff3 are produced by scripts/generate_from_genbank.py.
GFF_PATH       = "dataset/genome_annotation.gff3"   # Reference annotation, GFF3
PATHOGEN_JSON  = "dataset/pathogen.json"            # Nextclade QC and alignment config
README_PATH    = "dataset/README.md"                # Dataset description shown in Nextclade
CHANGELOG_PATH = "dataset/CHANGELOG.md"             # Release history
REFERENCE_PATH = "dataset/reference.fasta"          # Reference genome, FASTA

AUSPICE_CONFIG   = "resources/auspice_config.json"      # Auspice / Nextclade tree display config
EXCLUDE          = "resources/exclude.txt"              # Sequences to keep out of the tree
CLADES           = "resources/clades.tsv"               # Clade-defining mutations
ACCESSION_STRAIN = "resources/accession_strain.tsv"     # Optional accession -> corrected strain name map
INCLUDE_EXAMPLES = "resources/include_examples.txt"     # Sequences to force into the example set
COLORS           = "resources/colors.tsv"               # Fixed colour assignments
COLORS_SCHEMES   = "resources/color_schemes.tsv"        # Palette lookup used by scripts/assign-colors.py
GENBANK_PATH     = "resources/reference.gbk"            # Annotated reference, GenBank (used by augur ancestral)
INFERRED_ANCESTOR = "resources/inferred-root.fasta"     # Precomputed ancestral root sequence

configfile: PATHOGEN_JSON

# =============================================================================
# 5. Included helpers
# =============================================================================
# Checks that the configuration above is usable, and the messages printed when
# the workflow starts, succeeds or fails. Neither defines any rule.

include: "scripts/config_checks.snkm"
include: "scripts/workflow_messages.snkm"

# =============================================================================
# 6. Workflow rules
# =============================================================================
# Each rule is one stage of the dataset build. Snakemake works out the order
# from the inputs and outputs, so the listing below is roughly, but not
# strictly, the order things run in.

rule all:
    """
    Default target: the packaged Nextclade dataset plus the Auspice tree you
    can drag onto auspice.us to inspect the build.
    """
    input:
        auspice = "results/auspice.json",
        dataset_zip = "dataset.zip",
        pathogen = "out-dataset/pathogen.json",
        examples = "results/example_sequences.fasta",
        **({"root": INFERRED_ANCESTOR} if STATIC_ANCESTRAL_INFERRENCE else {})


if FETCH_SEQUENCES == True:
    rule fetch:
        """
        Download and curate sequences and metadata from NCBI by running the
        ingest/ sub-workflow. Configure the taxon in ingest/config/config.yaml.
        """
        message:
            "Fetching sequences and metadata from NCBI via the ingest sub-workflow"
        output:
            sequences = SEQUENCES,
            metadata = METADATA
        threads: workflow.cores
        shell:
            """
            set -euo pipefail
            (cd ingest && snakemake --cores {threads} all)
            """

rule curate:
    message:
        """
        Cleaning up metadata with augur merge & augur curate
        """
    input:
        meta = METADATA,  # Path to input metadata file
        strains = ACCESSION_STRAIN  # Strain - accession lookup table
    params:
        strain_id_field = ID_FIELD,
    output:
        metadata = "results/metadata.tsv",  # Final output file for publications metadata
    shell:
        """
        augur merge --metadata metadata={input.meta} strains={input.strains}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata results/metadata.merged.tmp
        augur curate normalize-strings \
            --metadata results/metadata.merged.tmp \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}

        rm results/metadata.merged.tmp
        """

rule add_reference_to_include:
    """
    Build the include list for augur filter: whatever the user listed in
    resources/include.txt, plus the reference and the inferred ancestor, which
    must never be dropped by subsampling.
    """
    input:
        include = "resources/include.txt",
    output:
        include = "results/include.txt",
    shell:
        """
        cat {input.include} > {output.include}
        echo "{REFERENCE_ACCESSION}" >> {output.include}
        echo ancestral_sequence >> {output.include}
        """

if STATIC_ANCESTRAL_INFERRENCE and INFERRENCE_RERUN:
    rule static_inferrence:
        message:
            """
            Running "inferred-root" snakefile for inference of the ancestral root. 
            This reference will be included in the Nextclade reference tree.
            WARNING: This will overwrite your {output.seq} & {output.meta} files!
            """
        input:
            dir = "inferred-root",
            dataset_path = "dataset",
            meta = rules.curate.output.metadata,
            seq = SEQUENCES,
            meta_ancestral = "resources/static_inferred_metadata.tsv",
            include = "results/include.txt"
        params:
            strain_id_field = ID_FIELD,
        output:
            inref = INFERRED_ANCESTOR,
            seq = INFERRED_SEQ_PATH,
            meta = INFERRED_META_PATH,
        threads: workflow.cores
        shell:
            r"""
            set -euo pipefail

            echo "Cleaning previous results..."
            rm -rf {input.dir}/results/* {input.dir}/resources/inferred-root.fasta

            echo "Running inferred-root workflow..."
            cd {input.dir}
            snakemake --cores {threads} all_sub
            cd - > /dev/null

            echo "Combining sequences with ancestral root..."
            cat {input.seq} {output.inref} > {output.seq}

            echo "Merging metadata..."
            augur merge \
                --metadata metadata={input.meta} ancestral={input.meta_ancestral} \
                --metadata-id-columns {params.strain_id_field} \
                --output-metadata {output.meta}

            echo "Static ancestral inference completed successfully!"
            """

if STATIC_ANCESTRAL_INFERRENCE and not INFERRENCE_RERUN:
    rule add_ancestral:
        input:
            meta = rules.curate.output.metadata,
            seq = SEQUENCES,
            meta_ancestral = "resources/static_inferred_metadata.tsv",
            inref = INFERRED_ANCESTOR,
        output:
            seq = INFERRED_SEQ_PATH,
            meta = INFERRED_META_PATH,
        params:
            strain_id_field = "accession",
        shell:
            """
            echo "Combining sequences with ancestral root..."
            cat {input.seq} {input.inref} > {output.seq}

            echo "Merging metadata..."
            augur merge \
                --metadata metadata={input.meta} ancestral={input.meta_ancestral} \
                --metadata-id-columns {params.strain_id_field} \
                --output-metadata {output.meta}

            echo "Static ancestral sequence imported successfully!"
            """

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = INFERRED_SEQ_PATH,
    output:
        sequence_index = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    """
    Cut the downloaded sequences down to a tree-sized set: drop anything older
    than {MIN_DATE} or shorter than {MIN_LENGTH}, then subsample to {MAX_SEQS}
    balanced across country and year, so one heavily sequenced outbreak does
    not swamp the tree. The reference and the ancestral sequence are always
    kept, via {input.include}.
    """
    input:
        sequences = INFERRED_SEQ_PATH,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        include = rules.add_reference_to_include.output,
    output:
        filtered_sequences = "results/filtered_sequences_raw.fasta",
        filtered_metadata = "results/filtered_metadata_raw.tsv",
    params: 
        min_date="" if MIN_DATE == "" else "--min-date " + MIN_DATE,
        min_length="" if MIN_LENGTH == "" else "--min-length " + MIN_LENGTH,
        max_seqs = MAX_SEQS,
        categories = "country year", #TODO: add subsampling per category?
        strain_id_field = ID_FIELD,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            {params.min_length} \
            {params.min_date} \
            --include {input.include} \
            --group-by {params.categories} \
            --subsample-max-sequences {params.max_seqs} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}
        """

rule align:
    """
    Align everything to the reference with Nextclade, using the alignment
    parameters from dataset/pathogen.json, so the tree is built with exactly
    the alignment the finished dataset will use.

    The {output.tsv} it writes also carries the per-sequence substitution
    counts that the get_outliers rule works from.
    """
    message:
        "Aligning sequences to the reference with Nextclade"
    input:
        sequences = rules.filter.output.filtered_sequences,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
    output:
        alignment = "results/aligned.fasta",
        tsv = "results/nextclade.tsv",
    params:
        translation_template = lambda w: "results/translations/cds_{cds}.translation.fasta",
        penalty_gap_extend = config["alignmentParams"]["penaltyGapExtend"],
        penalty_gap_open = config["alignmentParams"]["penaltyGapOpen"],
        penalty_gap_open_in_frame = config["alignmentParams"]["penaltyGapOpenInFrame"],
        penalty_gap_open_out_of_frame = config["alignmentParams"]["penaltyGapOpenOutOfFrame"],
        kmer_length = config["alignmentParams"]["kmerLength"],
        kmer_distance = config["alignmentParams"]["kmerDistance"],
        min_match_length = config["alignmentParams"]["minMatchLength"],
        allowed_mismatches = config["alignmentParams"]["allowedMismatches"],
        min_length = config["alignmentParams"]["minLength"],
        gap_alignment_side = config["alignmentParams"]["gapAlignmentSide"],  
        min_seed_cover = config["alignmentParams"]["minSeedCover"],
    threads: workflow.cores
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.annotation} \
        --alignment-preset high-diversity \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --gap-alignment-side {params.gap_alignment_side} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-seed-cover {params.min_seed_cover} \
        --min-length {params.min_length} \
        --max-alignment-attempts 5 \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations {params.translation_template} \
        --output-fasta {output.alignment} 
        """


rule get_outliers:
    """
    Flag sequences carrying more than {ALLOWED_DIVERGENCE} substitutions, or no
    usable substitution count at all. These are nearly always sequencing
    errors, misannotated records, or a different virus, and left in they show
    up as implausibly long branches.

    The list lands in {output.outliers} and feeds the exclude rule. Review
    it: if real diversity is being thrown away, raise {ALLOWED_DIVERGENCE}.
    """
    input:
        nextclade = rules.align.output.tsv,
    output:
        outliers = "results/outliers.txt",
        tmp = "tmp/outliers.txt",
    params:
        allowed_divergence = lambda w: ALLOWED_DIVERGENCE,
    shell:
        """
        tsv-filter -H -v --is-numeric totalSubstitutions {input.nextclade} \
        > {output.tmp}
        tsv-filter -H \
            --is-numeric totalSubstitutions \
            --gt totalSubstitutions:{params.allowed_divergence} \
            {input.nextclade} \
        | tail -n +2 >> {output.tmp}
        cat {output.tmp} \
        | tsv-select -H -f seqName \
        | tail -n +2 > {output.outliers}
        """


rule exclude:
    """
    Drop the manually excluded sequences ({input.exclude}), the
    automatically detected outliers, and the example sequences, which are
    shipped separately and should not also sit in the tree.

    This happens after alignment rather than during filtering on purpose.
    Removing sequences before subsampling would let the subsample pull in
    replacements, which could themselves be outliers, and round it would go.
    """
    input:
        sequences = rules.align.output.alignment,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        example = INCLUDE_EXAMPLES,
        include = rules.add_reference_to_include.output, ##TODO: check if you would like to keep these sequences
    params:
        strain_id_field = ID_FIELD,
    output:
        filtered_sequences = "results/filtered_aligned.fasta",
        filtered_metadata = "results/filtered_metadata.tsv",
        strains = "results/tree_strains.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} {input.outliers} {input.example} \
            --include {input.include} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains}
        """

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.exclude.output.filtered_sequences,
    output:
        tree = "results/tree_raw.nwk",
    threads: workflow.cores
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree} \
        """

rule refine:
    """
    Turn the raw tree into the reference tree: place the root where ROOTING
    says, and annotate branch lengths in mutations per site. Polytomies are
    kept rather than resolved arbitrarily, since a tree with no time signal
    should not pretend to have one.
    """
    message:
        "Rooting and annotating the tree (root: " + ROOTING + ")"
    input:
        tree = rules.tree.output.tree,
        alignment = rules.exclude.output.filtered_sequences,
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root {ROOTING} \
            --keep-polytomies \
            --divergence-unit mutations-per-site \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """


rule ancestral:
    """
    Reconstruct sequences at every internal node and translate the CDS listed
    in GENES. This is what lets Nextclade report mutations relative to a
    sequence's parent or to its clade founder, not just to the reference.

    This is the rule that fails if the CDS names in GENES, the GFF3 and the
    GenBank file do not agree. See step 1 of the README.
    """
    message:
        "Reconstructing ancestral sequences and translating CDS"
    input:
        tree=rules.refine.output.tree,
        alignment=rules.exclude.output.filtered_sequences,
        annotation=GENBANK_PATH,
    output:
        node_data="results/muts.json",
        ancestral_sequences="results/ancestral_sequences.fasta",
    params:
        translation_template=r"results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"results/translations/cds_%GENE.ancestral.fasta",
        genes=" ".join(GENES),
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --root-sequence {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template}\
            --output-sequences {output.ancestral_sequences}
        """

rule clades:
    """
    Assign each node a clade, by walking the tree and matching the
    clade-defining mutations in {input.clades} against the reconstructed
    mutations. Nodes matching nothing come out as "unassigned".
    """
    message:
        "Assigning clades from " + CLADES
    input:
        tree = rules.refine.output.tree,
        mutations = rules.ancestral.output.node_data,
        clades = CLADES
    output:
        json = "results/clades.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output-node-data {output.json}
        """

rule get_dates:
    """
    List the collection dates present in the data, in order, as input for
    colour assignment. Sorting them means the date colour scale runs oldest to
    newest rather than in whatever order the metadata happened to be in.
    """
    input:
        metadata = rules.exclude.output.filtered_metadata
    output:
        ordering = "results/color_ordering.tsv"
    run:
        import sys
        import pandas as pd

        column = "date"
        meta = pd.read_csv(input.metadata, delimiter='\t')

        if column not in meta.columns:
            print(f"The column '{column}' does not exist in the file.")
            sys.exit(1)

        deflist = meta[column].dropna().tolist()
        # Store unique values (ordered)
        deflist = sorted(set(deflist))
        if "XXXX-XX-XX" in deflist:
            deflist.remove("XXXX-XX-XX")

        result_df = pd.DataFrame({
            'column': ['date'] * len(deflist),
            'value': deflist
        })

        result_df.to_csv(output.ordering, sep='\t', index=False, header=False)
        

rule colors:
    """
    Build the colour table for Auspice: a graded scale for dates, plus the
    fixed assignments in {input.colors} for everything else. Sequences
    with an unknown date get a neutral grey.
    """
    input:
        ordering = rules.get_dates.output.ordering,
        color_schemes = COLORS_SCHEMES,
        colors = COLORS,
    output:
        colors = "results/colors_dates.tsv",
        final_colors = "results/final_colors.tsv"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors}

        echo -e '\ndate\tXXXX-XX-XX\t#a6acaf' >> {output.colors}

        cat {output.colors} {input.colors} > {output.final_colors}
        """

rule export:
    """
    Write the Auspice v2 JSON. This doubles as the Nextclade reference tree and
    as something you can open on auspice.us to inspect the build.
    """
    message:
        "Exporting the reference tree"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.exclude.output.filtered_metadata,
        mutations = rules.ancestral.output.node_data,
        branch_lengths = rules.refine.output.node_data,
        clades = rules.clades.output.json, # dummy_clades if not set yet
        auspice_config = AUSPICE_CONFIG,
        colors = rules.colors.output.final_colors,
    params:
        strain_id_field = ID_FIELD,
        fields = "region country date",
    output:
        auspice = "results/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --auspice-config {input.auspice_config} \
            --color-by-metadata {params.fields} \
            --colors {input.colors} \
            --include-root-sequence \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice}
        """


rule subsample_example_sequences:
    """
    Pick a small, recent, well-dated set of sequences to ship with the dataset.
    These are what Nextclade's "Load example" button loads, so they should
    align and pass QC cleanly: anything listed in exclude.txt or flagged as an
    outlier is kept out, and tree sequences are fine to reuse here.
    """
    message:
        "Selecting example sequences to ship with the dataset"
    input:
        all_sequences = INFERRED_SEQ_PATH,
        metadata = INFERRED_META_PATH,
        incl_examples = INCLUDE_EXAMPLES,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        tree_strains = "results/tree_strains.txt",  # strains in the tree
    output:
        example_sequences = "results/example_sequences.fasta",
    params:
        strain_id_field = ID_FIELD,
    shell:
        """
        augur filter \
            --sequences {input.all_sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --min-length 4000 \
            --include {input.incl_examples} \
            --exclude {input.exclude} {input.outliers} \
            --exclude-ambiguous-dates-by year \
            --min-date 2015 \
            --subsample-max-sequences 25  \
            --probabilistic-sampling \
            --output-sequences {output.example_sequences}
        """

rule assemble_dataset:
    """
    Collect the finished dataset files into out-dataset/ and pack them into
    dataset.zip, which is what you hand to Nextclade CLI or upload for review.
    """
    message:
        "Assembling the Nextclade dataset in out-dataset/ and zipping it"
    input:
        tree = rules.export.output.auspice,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
        sequences = rules.subsample_example_sequences.output.example_sequences,
        pathogen = PATHOGEN_JSON,
        readme = README_PATH,
        changelog = CHANGELOG_PATH,
    params:
        pathogen = "out-dataset/pathogen.json",
    output:
        tree = "out-dataset/tree.json",
        reference = "out-dataset/reference.fasta",
        annotation = "out-dataset/genome_annotation.gff3",
        sequences = "out-dataset/sequences.fasta",
        readme = "out-dataset/README.md",
        changelog = "out-dataset/CHANGELOG.md",
        dataset_zip = "dataset.zip",
    shell:
        """
        set -euo pipefail

        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {params.pathogen}
        cp {input.readme} {output.readme}
        cp {input.changelog} {output.changelog}
        zip -rj dataset.zip  out-dataset/*
        """

rule extract_clades_tsv:
    """
    Flatten the clade assignments into a two-column TSV, so the mutLabels rule
    can join them onto the Nextclade mutation calls.
    """
    input:
        json = rules.clades.output.json,
    output:
        tsv = "results/clades_metadata.tsv"
    run:
        import json
        import csv

        with open(input.json) as f:
            data = json.load(f)

        nodes = data.get("nodes", {})

        with open(output.tsv, "w", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerow(["accession", "clade"])

            for accession, values in nodes.items():
                clade = values.get("clade_membership", None)
                if clade:
                    writer.writerow([accession, clade])

rule mutLabels:
    """
    Work out which nucleotide substitutions are characteristic of each clade
    and merge them into pathogen.json as mutation labels. Nextclade uses these
    to tell "labelled" from "private" mutations when scoring QC.
    """
    message:
        "Deriving per-clade mutation labels and merging them into pathogen.json"
    input:
        table = rules.align.output.tsv,
        clade = rules.extract_clades_tsv.output.tsv,
        json = PATHOGEN_JSON,
    params:
        # A mutation is labelled for a clade if at least this fraction of the
        # clade's sequences carry it.
        min_proportion = 0.2,
        # Stricter threshold for large, well-sampled clades, where a low cutoff
        # would label a lot of noise.
        high_threshold_proportion = 0.60,
        # CHANGE ME: the clades that should use the stricter threshold above.
        clades_high_threshold = ["A", "B", "C"],
        # Clades to ignore entirely; "unassigned" collects everything that did
        # not fall into a defined clade.
        clades_to_drop = ["unassigned"],
    output:
        clade_meta = "results/clades_mut_metadata.tsv",
        properties = "results/virus_properties.json",
        json = "out-dataset/pathogen.json",
        newly_relevant = "results/newly_relevant_mutations.tsv",
    shell:
        """
        set -euo pipefail

        augur merge \
            --metadata meta={input.table} clade={input.clade} \
            --metadata-id-columns meta=seqName clade=accession \
            --output-metadata {output.clade_meta}

        python3 scripts/generate_virus_properties.py \
            --clade_meta {output.clade_meta} \
            --properties {output.properties} \
            --min-prop {params.min_proportion} \
            --high-min-prop {params.high_threshold_proportion} \
            --high-prop-clades "{params.clades_high_threshold}" \
            --exclude-clades "{params.clades_to_drop}" \
            --newly-relevant-output {output.newly_relevant}

        jq --slurpfile v {output.properties} \
           '.mutLabels.nucMutLabelMap = $v[0].nucMutLabelMap |
            .mutLabels.nucMutLabelMapReverse = $v[0].nucMutLabelMapReverse' \
           {input.json} > {output.json}

        zip -rj dataset.zip  out-dataset/*
        """

rule clean:
    """
    Remove everything the workflow generates, including the downloaded
    sequences. Leaves resources/inferred-root.fasta alone, since regenerating
    it is a deliberate, separate step.
    """
    shell:
        """
        rm -rf results out-dataset tmp dataset.zip
        rm -rf data ingest/data ingest/logs ingest/benchmarks
        """

rule viz:
    input: "results/auspice.json"
    shell: "auspice view --datasetDir results"

rule serve:
    input: "out-dataset/tree.json","results/virus_properties.json"
    params: "out-dataset"
    shell: "serve --cors {params} -l 3000"