rule download:
    message: "Downloading metadata and fasta files from S3"
    output:
        sequences = config["sequences"],
        metadata = config["metadata"]
    conda: config["conda_environment"]
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz - | gunzip -cq >{output.metadata:q}
        aws s3 cp s3://nextstrain-ncov-private/sequences.fasta.gz - | gunzip -cq > {output.sequences:q}
        """
from pathlib import Path
merged_data = Path(outdir).parent.joinpath("data")


rule merge_input:
    message: "Merging denmark and global data to prepare input for the workflow"
    input:
        denmark_fasta = config["denmark_fasta"],
        denmark_meta = config["denmark_meta"],
        gisaid_fasta = config["gisaid_fasta"],
        gisaid_meta = config["gisaid_meta"]
    output:
        metadata = f"{merged_data}/metadata_nextstrain_DKglobal.tsv",
        sequences = f"{merged_data}/merged_sequences.fasta"
    conda: config["conda_environment"]
    params:
        mergedir = merged_data
    shell:
        """
        Rscript --vanilla --no-environ /opt/workflows/merge_clean_metadata.R -l {input.denmark_meta} -g {input.gisaid_meta} -o {params.mergedir}
        awk 'NR > 1 {{print $1}}' {params.mergedir}/metadata_nextstrain_DKglobal.tsv |sort|uniq > {params.mergedir}/include.txt
        # Dedup fasta: https://www.biostars.org/p/143617/#466790
        cat {input.denmark_fasta} {input.gisaid_fasta}|awk '/^>/{{f=!d[$1];d[$1]=1}}f' | seqtk subseq - {params.mergedir}/include.txt > {output.sequences}
        """


rule copy_input:
    message: "Copy Denmark data to place to prepare input for the workflow"
    input:
        denmark_fasta = config["denmark_fasta"],
        denmark_meta = config["denmark_meta"],
        rootseqs = "my_profiles/denmark/rootseqs.fasta"
    output:
        metadata = f"{merged_data}/{Path(config['denmark_meta']).name}",
        sequences = f"{merged_data}/{Path(config['denmark_fasta']).name}"
    conda: config["conda_environment"]
    params:
        mergedir = merged_data
    shell:
        """
        cp {input.denmark_meta} {output.metadata}
        cat {input.denmark_fasta} {input.rootseqs} > {output.sequences}
        """

rule prepare_input:
    message: "Merging denmark and global data to prepare input for the workflow"
    input:
        fasta = rules.merge_input.output.sequences if "Global" in BUILD_NAMES else rules.copy_input.output.sequences, 
        meta = rules.merge_input.output.metadata if "Global" in BUILD_NAMES else rules.copy_input.output.metadata
    output:
        metadata = f"{merged_data}/input_metadata.tsv",
        sequences = f"{merged_data}/input_sequences.fasta"
    conda: config["conda_environment"]
    params:
        mergedir = merged_data
    shell:
        """
        ln -sfr {input.meta} {output.metadata}
        touch -h {output.metadata}
        ln -sfr {input.fasta} {output.sequences}
        touch -h {output.sequences}

        """
from datetime import date, timedelta
from treetime.utils import numeric_date

rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.prepare_input.output.sequences,
        metadata = rules.prepare_input.output.metadata,
        include = config["files"]["include"],
        exclude = config["files"]["exclude"]
    output:
        sequences = f"{outdir}/filtered.fasta"
    log:
        "logs/filtered.txt"
    params:
        min_length = config["filter"]["min_length"],
        exclude_where = config["filter"]["exclude_where"],
        min_date = config["filter"]["min_date"],
        date = date.today().strftime("%Y-%m-%d")
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --max-date {params.date} \
            --min-date {params.min_date} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where}\
            --min-length {params.min_length} \
            --output {output.sequences} 2>&1 | tee {log}
        """

rule excluded_sequences:
    message:
        """
        Generating fasta file of excluded sequences
        """
    input:
        sequences = rules.prepare_input.output.sequences,
        metadata = rules.prepare_input.output.metadata,
        include = config["files"]["exclude"]
    output:
        sequences = f"{outdir}/excluded.fasta"
    log:
        "logs/excluded.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
	    --min-length 50000 \
            --include {input.include} \
            --output {output.sequences} 2>&1 | tee {log}
        """

rule align_excluded:
    message:
        """
        Aligning excluded sequences to {input.reference}
          - gaps relative to reference are considered real
        """
    input:
        sequences = rules.excluded_sequences.output.sequences,
        reference = config["files"]["reference"]
    output:
        alignment = f"{outdir}/excluded_alignment.fasta"
    log:
        "logs/align_excluded.txt"
    threads: 2
    conda: config["conda_environment"]
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference 2>&1 | tee {log}
        """

rule diagnose_excluded:
    message: "Scanning excluded sequences {input.alignment} for problematic sequences"
    input:
        alignment = rules.align_excluded.output.alignment,
        metadata = rules.prepare_input.output.metadata,
        reference = config["files"]["reference"]
    output:
        diagnostics = f"{outdir}/excluded-sequence-diagnostics.tsv",
        flagged = f"{outdir}/excluded-flagged-sequences.tsv",
        to_exclude = f"{outdir}/check_exclusion.txt"
    log:
        "logs/diagnose-excluded.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/diagnostic.py \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --reference {input.reference} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --output-flagged {output.flagged} \
            --output-diagnostics {output.diagnostics} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """


checkpoint partition_sequences:
    input:
        sequences = rules.filter.output.sequences
    output:
        split_sequences = directory(f"{outdir}/split_sequences/")
    log:
        "logs/partition_sequences.txt"
    params:
        sequences_per_group = config["partition_sequences"]["sequences_per_group"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/partition-sequences.py \
            --sequences {input.sequences} \
            --sequences-per-group {params.sequences_per_group} \
            --output-dir {output.split_sequences} 2>&1 | tee {log}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - gaps relative to reference are considered real
        Cluster:  {wildcards.cluster}
        """
    input:
        sequences = outdir + "/split_sequences/"+"{cluster}.fasta",
        reference = config["files"]["reference"]
    output:
        alignment = outdir + "/split_alignments/"+"{cluster}.fasta"
    log:
        "logs/align_{cluster}.txt"
    benchmark:
        "benchmarks/align_{cluster}.txt"
    threads: 2
    conda: config["conda_environment"]
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference 2>&1 | tee {log}
        """

def _get_alignments(wildcards):
    checkpoint_output = checkpoints.partition_sequences.get(**wildcards).output[0]
    return expand(outdir + "/split_alignments/{i}.fasta",
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)

rule aggregate_alignments:
    message: "Collecting alignments"
    input:
        alignments = _get_alignments
    output:
        alignment = f"{outdir}/aligned.fasta"
    log:
        "logs/aggregate_alignments.txt"
    conda: config["conda_environment"]
    shell:
        """
        cat {input.alignments} > {output.alignment} 2> {log}
        """

rule diagnostic:
    message: "Scanning aligned sequences {input.alignment} for problematic sequences"
    input:
        alignment = rules.aggregate_alignments.output.alignment,
        metadata = rules.prepare_input.output.metadata,
        reference = config["files"]["reference"]
    output:
        diagnostics = f"{outdir}/sequence-diagnostics.tsv",
        flagged = f"{outdir}/flagged-sequences.tsv",
        to_exclude = f"{outdir}/to-exclude.txt"
    log:
        "logs/diagnostics.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/diagnostic.py \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --reference {input.reference} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --output-flagged {output.flagged} \
            --output-diagnostics {output.diagnostics} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """

rule refilter:
    message:
        """
        excluding sequences flagged in the diagnostic step in file {input.exclude}
        """
    input:
        sequences = rules.aggregate_alignments.output.alignment,
        metadata = rules.prepare_input.output.metadata,
        exclude = rules.diagnostic.output.to_exclude
    output:
        sequences = f"{outdir}/aligned-filtered.fasta"
    log:
        "logs/refiltered.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} 2>&1 | tee {log}
        """


rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.refilter.output.sequences
    output:
        alignment = f"{outdir}/masked.fasta"
    log:
        "logs/mask.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"],
        mask_sites = config["mask"]["mask_sites"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --mask-terminal-gaps \
            --output {output.alignment} 2>&1 | tee {log}
        """

def _get_subsampling_settings(wildcards):
    # Allow users to override default subsampling with their own settings keyed
    # by location type and name. For example, "region_europe" or
    # "country_iceland". Otherwise, default to settings for the location type.
    subsampling_scheme = _get_subsampling_scheme_by_build_name(wildcards.build_name)
    subsampling_settings = config["subsampling"][subsampling_scheme]

    if hasattr(wildcards, "subsample"):
        subsampling_settings = subsampling_settings[wildcards.subsample]

        # If users have supplied both `max_sequences` and `seq_per_group`, we
        # throw an error instead of assuming the user prefers one setting over
        # another by default.
        if subsampling_settings.get("max_sequences") and subsampling_settings.get("seq_per_group"):
            raise Exception(f"The subsampling scheme '{subsampling_scheme}' for build '{wildcards.build_name}' defines both `max_sequences` and `seq_per_group`, but these arguments are mutually exclusive. If you didn't define both of these settings, this conflict could be caused by using the same subsampling scheme name as a default scheme. In this case, rename your subsampling scheme, '{subsampling_scheme}', to a unique name (e.g., 'custom_{subsampling_scheme}') and run the workflow again.")

        # If users have supplied neither `max_sequences` nor `seq_per_group`, we
        # throw an error because the subsampling rule will still group by one or
        # more fields and the lack of limits on this grouping could produce
        # unexpected behavior.
        if not subsampling_settings.get("max_sequences") and not subsampling_settings.get("seq_per_group"):
            raise Exception(f"The subsampling scheme '{subsampling_scheme}' for build '{wildcards.build_name}' must define `max_sequences` or `seq_per_group`.")

    return subsampling_settings


def get_priorities(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return outdir + f"/{wildcards.build_name}/proximity_{subsampling_settings['priorities']['focus']}.tsv"
    else:
        # TODO: find a way to make the list of input files depend on config
        return config["files"]["include"]


def get_priority_argument(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return "--priority " + get_priorities(wildcards)
    else:
        return ""


def _get_specific_subsampling_setting(setting, optional=False):
    def _get_setting(wildcards):
        if optional:
            value = _get_subsampling_settings(wildcards).get(setting, "")
        else:
            value = _get_subsampling_settings(wildcards)[setting]

        if isinstance(value, str):
            # Load build attributes including geographic details about the
            # build's region, country, division, etc. as needed for subsampling.
            build = config["builds"][wildcards.build_name]
            if (setting == "min_date" or setting == "max_date") and value:
                value = eval(value)
                return f"--{setting.replace('_', '-')} {value}"
            value = value.format(**build)
        elif value is not None:
            # If is 'seq_per_group' or 'max_sequences' build subsampling setting,
            # need to return the 'argument' for augur
            if setting == 'seq_per_group':
                value = f"--sequences-per-group {value}"
            elif setting == 'max_sequences':
                value = f"--subsample-max-sequences {value}"

            return value
        else:
            value = ""

        # Check format strings that haven't been resolved.
        if re.search(r'\{.+\}', value):
            raise Exception(f"The parameters for the subsampling scheme '{wildcards.subsample}' of build '{wildcards.build_name}' reference build attributes that are not defined in the configuration file: '{value}'. Add these build attributes to the appropriate configuration file and try again.")

        return value

    return _get_setting

rule subsample:
    message:
        """
        Subsample all sequences by '{wildcards.subsample}' scheme for build '{wildcards.build_name}' with the following parameters:

         - group by: {params.group_by}
         - sequences per group: {params.sequences_per_group}
         - subsample max sequences: {params.subsample_max_sequences}
         - min-date: {params.min_date}
         - max-date: {params.max_date}
         - exclude: {params.exclude_argument}
         - include: {params.include_argument}
         - query: {params.query_argument}
         - priority: {params.priority_argument}
        """
    input:
        sequences = rules.mask.output.alignment,
        metadata = rules.prepare_input.output.metadata,
        include = config["files"]["include"],
        priorities = get_priorities
    output:
        sequences = outdir + "/{build_name}/sample-{subsample}.fasta"
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    params:
        group_by = _get_specific_subsampling_setting("group_by"),
        sequences_per_group = _get_specific_subsampling_setting("seq_per_group", optional=True),
        subsample_max_sequences = _get_specific_subsampling_setting("max_sequences", optional=True),
        exclude_argument = _get_specific_subsampling_setting("exclude", optional=True),
        include_argument = _get_specific_subsampling_setting("include", optional=True),
        query_argument = _get_specific_subsampling_setting("query", optional=True),
        min_date = _get_specific_subsampling_setting("min_date", optional=True),
        max_date = _get_specific_subsampling_setting("max_date", optional=True),
        priority_argument = get_priority_argument
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            {params.min_date} \
            {params.max_date} \
            {params.exclude_argument} \
            {params.include_argument} \
            {params.query_argument} \
            {params.priority_argument} \
            --group-by {params.group_by} \
            {params.sequences_per_group} \
            {params.subsample_max_sequences} \
            --output {output.sequences} 2>&1 | tee {log}
        """

rule proximity_score:
    message:
        """
        determine priority for inclusion in as phylogenetic context by
        genetic similiarity to sequences in focal set for build '{wildcards.build_name}'.
        """
    input:
        alignment = rules.mask.output.alignment,
        metadata = rules.prepare_input.output.metadata,
        reference = config["files"]["reference"],
        focal_alignment = outdir + "/{build_name}/sample-{focus}.fasta"
    output:
        priorities = outdir + "/{build_name}/proximity_{focus}.tsv"
    log:
        "logs/subsampling_priorities_{build_name}_{focus}.txt"
    resources:
        mem_mb = 4000
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/priorities.py --alignment {input.alignment} \
            --metadata {input.metadata} \
            --reference {input.reference} \
            --focal-alignment {input.focal_alignment} \
            --output {output.priorities} 2>&1 | tee {log}
        """

def _get_subsampled_files(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    return [
        f"{outdir}/{wildcards.build_name}/sample-{subsample}.fasta"
        for subsample in subsampling_settings
    ]

rule combine_samples:
    message:
        """
        Combine and deduplicate FASTAs
        """
    input:
        _get_subsampled_files
    output:
        alignment = outdir + "/{build_name}/subsampled_alignment.fasta"
    log:
        "logs/subsample_regions_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py \
            --input {input} \
            --output {output} 2>&1 | tee {log}
        """

# TODO: This will probably not work for build names like "country_usa" where we need to know the country is "USA".
rule adjust_metadata_regions:
    message:
        """
        Adjusting metadata for build '{wildcards.build_name}'
        """
    input:
        metadata = rules.prepare_input.output.metadata
    output:
        metadata = outdir + "/{build_name}/metadata_adjusted.tsv"
    params:
        region = lambda wildcards: config["builds"][wildcards.build_name]["region"]
    log:
        "logs/adjust_metadata_regions_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/adjust_regional_meta.py \
            --region {params.region:q} \
            --metadata {input.metadata} \
            --output {output.metadata} 2>&1 | tee {log}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.combine_samples.output.alignment
    output:
        tree = outdir + "/{build_name}/tree_raw.nwk"
    params:
        args =  config["tree"].get("tree-builder-args","") if "tree" in config else "",
        method = config['tree'].get("method", "iqtree") if "tree" in config else ""

    log:
        "logs/tree_{build_name}.txt"
    benchmark:
        "benchmarks/tree_{build_name}.txt"
    threads: 16
    resources:
        # Multiple sequence alignments can use up to 40 times their disk size in
        # memory, especially for larger alignments.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 40 * int(input.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            --nthreads {threads} 2>&1 \
            --method {params.method} | tee {log}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.combine_samples.output.alignment,
        metadata = _get_metadata_by_wildcards
    output:
        tree = outdir + "/{build_name}/tree.nwk",
        node_data = outdir + "/{build_name}/branch_lengths.json"
    log:
        "logs/refine_{build_name}.txt"
    benchmark:
        "benchmarks/refine_{build_name}.txt"
    threads: 1
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    params:
        root = config["refine"]["root"],
        clock_rate = config["refine"]["clock_rate"],
        clock_std_dev = config["refine"]["clock_std_dev"],
        coalescent = config["refine"]["coalescent"],
        date_inference = config["refine"]["date_inference"],
        divergence_unit = config["refine"]["divergence_unit"],
        clock_filter_iqd = config["refine"]["clock_filter_iqd"],
        keep_polytomies = "--keep-polytomies" if config["refine"].get("keep_polytomies", False) else "",
        timetree = "" if config["refine"].get("no_timetree", False) else "--timetree"
    conda: config["conda_environment"]
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            {params.timetree} \
            {params.keep_polytomies} \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --divergence-unit {params.divergence_unit} \
            --date-confidence \
            --no-covariance \
            --clock-filter-iqd {params.clock_filter_iqd} 2>&1 | tee {log}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.combine_samples.output.alignment
    output:
        node_data = outdir + "/{build_name}/nt_muts.json"
    log:
        "logs/ancestral_{build_name}.txt"
    params:
        inference = config["ancestral"]["inference"]
    conda: config["conda_environment"]
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous 2>&1 | tee {log}
        """

rule haplotype_status:
    message: "Annotating haplotype status relative to {params.reference_node_name}"
    input:
        nt_muts = rules.ancestral.output.node_data
    output:
        node_data = outdir + "/{build_name}/haplotype_status.json"
    log:
        "logs/haplotype_status_{build_name}.txt"
    params:
        reference_node_name = config["reference_node_name"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/annotate-haplotype-status.py \
            --ancestral-sequences {input.nt_muts} \
            --reference-node-name {params.reference_node_name:q} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = config["files"]["reference"]
    output:
        node_data = outdir + "/{build_name}/aa_muts.json"
    log:
        "logs/translate_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        node_data = outdir + "/{build_name}/traits.json"
    log:
        "logs/traits_{build_name}.txt"
    params:
        columns = _get_trait_columns_by_wildcards,
        sampling_bias_correction = _get_sampling_bias_correction_for_wildcards
    conda: config["conda_environment"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction} 2>&1 | tee {log}
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["files"]["clades"]
    output:
        ncov_clade = outdir + "/{build_name}/clades.json"
    log:
        "logs/clades_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.ncov_clade} 2>&1 | tee {log}
        """

rule pangolin_clade:
    message: "Running pangolin 2, local clade calculation"
    input:
        alignment = rules.combine_samples.output.alignment
    output:
        clades = outdir + "/{build_name}/lineage_report.csv"
    threads: 60
    params:
        tmpdir = outdir + "/{build_name}/pangtmp"
    shell:
        """
       bash -c "
        source activate pangolin
        pangolin {input.alignment} -t {threads} \
        --tempdir {params.tmpdir} \
        --outfile {output.clades}
        rm -rf {params.tmpdir}
        "
        """
rule pangolin:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        clades = rules.pangolin_clade.output.clades
    output:
        ncov_clade = outdir + "/{build_name}/pangolin.json"
    log:
        "logs/pangolin_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/add_pangolin_lineages.py \
            --clades {input.clades} \
            --tree {input.tree} \
            --output {output.ncov_clade}
        """


rule legacy_clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["files"]["legacy_clades"]
    output:
        ncov_clade = outdir + "/{build_name}/temp_legacy_clades.json"
    log:
        "logs/legacy_clades_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.ncov_clade} 2>&1 | tee {log}
        """

rule rename_legacy_clades:
    input:
        node_data = rules.legacy_clades.output.ncov_clade
    output:
        ncov_clade = outdir + "/{build_name}/legacy_clades.json"
    run:
        import json
        with open(input.node_data, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
            new_data = {}
            for k,v in d['nodes'].items():
                if "clade_membership" in v:
                    new_data[k] = {"legacy_clade_membership": v["clade_membership"]}
        with open(output.ncov_clade, "w") as fh:
            json.dump({"nodes":new_data}, fh)

rule resolve_clades:
    message: "Resolve subclades "
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["files"]["clades"],
        metadata = _get_metadata_by_wildcards,
        reference = config["files"]["reference"],
        alignment = rules.mask.output.alignment
    output:
        clades_dk = outdir + "/{build_name}/resolved_clades.json",
        tip_cluster = outdir + "/{build_name}/tip_cluster.json",
        new_clades = outdir + "/{build_name}/new_clades.tsv"
    log:
        "logs/resolve_clades_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur resolve-clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --alignment {input.alignment} \
            --reference {input.reference} \
            --metadata {input.metadata} \
            --new-clades {output.new_clades} \
            --max-depth 3 \
            --output-tip-cluster {output.tip_cluster} \
            --output-node-data {output.clades_dk} 2>&1 | tee {log}
        """

rule subclades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        subclades = config['files']['subclades'],
        clades = config["files"]["clades"]
    output:
        ncov_clade = outdir + "/{build_name}/temp_subclades.json"
    params:
        clade_file = outdir + "/{build_name}/temp_subclades.tsv"
    log:
        "logs/subclades_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        cat {input.clades} {input.subclades} > {params.clade_file} && \
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {params.clade_file} \
            --output-node-data {output.ncov_clade} 2>&1 | tee {log}
        """

rule rename_subclades_dk:
    input:
        clades_dk = rules.resolve_clades.output.clades_dk
    output:
        clades_dk = outdir + "/{build_name}/subclades_dk.json"
    run:
        import json
        with open(input.clades_dk, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
            new_data = {}
            for k,v in d['nodes'].items():
                if "clade_membership" in v:
                    new_data[k] = {"subclade_dk": v["clade_membership"]}
        with open(output.clades_dk, "w") as fh:
            json.dump({"nodes":new_data}, fh, indent=2)

rule rename_subclades:
    input:
        node_data = rules.subclades.output.ncov_clade
    output:
        ncov_clade = outdir + "/{build_name}/subclades.json"
    run:
        import json
        with open(input.node_data, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
            new_data = {}
            for k,v in d['nodes'].items():
                if "clade_membership" in v:
                    new_data[k] = {"subclade_membership": v["clade_membership"]}
        with open(output.ncov_clade, "w") as fh:
            json.dump({"nodes":new_data}, fh)


rule colors:
    message: "Constructing colors file"
    input:
        ordering = config["files"]["ordering"],
        color_schemes = config["files"]["color_schemes"],
        metadata = _get_metadata_by_wildcards
    output:
        colors = outdir + "/{build_name}/colors.tsv"
    log:
        "logs/colors_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata} 2>&1 | tee {log}
        """

rule recency:
    message: "Use metadata on submission date to construct submission recency field"
    input:
        metadata = _get_metadata_by_wildcards
    output:
        node_data = outdir + "/{build_name}/recency.json"
    log:
        "logs/recency_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --output {output} 2>&1 | tee {log}
        """

rule tip_frequencies:
    message: "Estimating censored KDE frequencies for tips"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        tip_frequencies_json = outdir + "/{build_name}/tip-frequencies.json"
    log:
        "logs/tip_frequencies_{build_name}.txt"
    params:
        min_date = config["frequencies"]["min_date"],
        pivot_interval = config["frequencies"]["pivot_interval"],
        narrow_bandwidth = config["frequencies"]["narrow_bandwidth"],
        proportion_wide = config["frequencies"]["proportion_wide"]
    conda: config["conda_environment"]
    shell:
        """
        augur frequencies \
            --method kde \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --min-date {params.min_date} \
            --pivot-interval {params.pivot_interval} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --output {output.tip_frequencies_json} 2>&1 | tee {log}
        """

rule nucleotide_mutation_frequencies:
    message: "Estimate nucleotide mutation frequencies"
    input:
        alignment = rules.combine_samples.output.alignment,
        metadata = _get_metadata_by_wildcards
    output:
        frequencies = outdir + "/{build_name}/nucleotide_mutation_frequencies.json"
    log:
        "logs/nucleotide_mutation_frequencies_{build_name}.txt"
    params:
        min_date = config["frequencies"]["min_date"],
        minimal_frequency = config["frequencies"]["minimal_frequency"],
        pivot_interval = config["frequencies"]["pivot_interval"],
        stiffness = config["frequencies"]["stiffness"],
        inertia = config["frequencies"]["inertia"]
    conda: config["conda_environment"]
    shell:
        """
        augur frequencies \
            --method diffusion \
            --alignments {input.alignment} \
            --gene-names nuc \
            --metadata {input.metadata} \
            --min-date {params.min_date} \
            --minimal-frequency {params.minimal_frequency} \
            --pivot-interval {params.pivot_interval} \
            --stiffness {params.stiffness} \
            --inertia {params.inertia} \
            --output {output.frequencies} 2>&1 | tee {log}
        """

def export_title(wildcards):
    # TODO: maybe we could replace this with a config entry for full/human-readable build name?
    location_name = wildcards.build_name

    # If specified in config file generally, or in a config file build
    if "title" in config["builds"][location_name]:
        return config["builds"][location_name]["title"]
    elif "title" in config:
        return config["title"]

    # Else return an auto-generated title
    if not location_name:
        return "Genomic epidemiology of novel coronavirus"
    elif location_name == "global":
        return "Genomic epidemiology of novel coronavirus - Global"
    else:
        location_title = location_name.replace("-", " ").title()
        return f"Genomic epidemiology of novel coronavirus - {location_title}"

def _get_tip_cluster_file(wildcards):
    # return rules.resolve_clades.output.new_clades if wildcards.build_name == 'Denmark' else "results/Denmark/new_clades.tsv"
    return checkpoints.resolve_clades.get(build_name = wildcards.build_name).output.tip_cluster

def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    # rules.pangolin.output.clade_data,
    # _get_tip_cluster_file(wildcards),
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.rename_legacy_clades.output.ncov_clade,
        rules.rename_subclades.output.ncov_clade,
        rules.clades.output.ncov_clade,
        rules.rename_subclades_dk.output.clades_dk,
        rules.recency.output.node_data,
        rules.traits.output.node_data
    ]

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs
rule update_description:
    message: "Prepare description file(s)"
    input:
        description = lambda w: config["builds"][w.build_name]["description"] if "description" in config["builds"][w.build_name] else config["files"]["description"]
    params:
        data_date = config["data_date"],
        gisaid_date = config["gisaid_date"]
    output:
        description = outdir + "/{build_name}/description.md"
    conda: config["conda_environment"]
    shell:
        """
        sed -e 's/__DenmarkDate__/{params.data_date}/;s/__GisaidDate__/{params.gisaid_date}/' {input.description} > {output.description}
        """

rule update_auspice_config:
    message: "Update auspice config to contain metadata columns"
    input:
        metadata = config['denmark_meta'],
        auspice_config = lambda w: config["builds"][w.build_name]["auspice_config"] if "auspice_config" in config["builds"][w.build_name] else config["files"]["auspice_config"],
    output:
        auspice_config = outdir + "/{build_name}/auspice_config_{build_name}.json"
    params:
        exclude_columns = config.get('exclude_columns', '')
    run:
        exclude_columns = params.exclude_columns.split(",")
        # exclude some dummy fields
        exclude_columns.extend(['strain', 'virus', 'host', 'country', 'region', 'date_submitted', 'region_exposure', 'division_exposure', 'length' ])
        columns  = list()
        with open(input.metadata, 'r', encoding='utf-8') as fh:
            columns = fh.readline().strip().split("\t")

        import json
        with open(input.auspice_config, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
            for c in columns:
                if c not in exclude_columns and c not in d["filters"]:
                    d["filters"].append(c)
        with open(output.auspice_config, "w") as fh:
            json.dump(d, fh, indent=2)

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        node_data = _get_node_data_by_wildcards,
        auspice_config = rules.update_auspice_config.output.auspice_config,
        colors = lambda w: config["builds"][w.build_name]["colors"] if "colors" in config["builds"][w.build_name] else ( config["files"]["colors"] if "colors" in config["files"] else rules.colors.output.colors.format(**w) ),
        lat_longs = config["files"]["lat_longs"],
        description = rules.update_description.output.description
    output:
        auspice_json = outdir + "/{build_name}/ncov_with_accessions.json"
        root_sequence_json = outdir + "/{build_name}/ncov_with_accessions_root-sequence.json"
    log:
        "logs/export_{build_name}.txt"
    params:
        title = export_title
    conda: config["conda_environment"]
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --title {params.title:q} \
            --description {input.description} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
rule incorporate_travel_history:
    message: "Adjusting main auspice JSON to take into account travel history"
    input:
        auspice_json = rules.export.output.auspice_json,
        colors = lambda w: config["builds"][w.build_name]["colors"] if "colors" in config["builds"][w.build_name] else ( config["files"]["colors"] if "colors" in config["files"] else rules.colors.output.colors.format(**w) ),
        lat_longs = config["files"]["lat_longs"]
    params:
        sampling = _get_sampling_trait_for_wildcards,
        exposure = _get_exposure_trait_for_wildcards
    output:
        auspice_json = outdir + "/{build_name}/ncov_with_accessions_and_travel_branches.json"
    log:
        "logs/incorporate_travel_history_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule finalize:
    message: "Remove extraneous colorings for main build and move frequencies"
    input:
        # Would probably need to return the original genome, not these masked and aligned ones
        genomes = rules.combine_samples.output.alignment,
        auspice_json = rules.incorporate_travel_history.output.auspice_json,
        frequencies = rules.tip_frequencies.output.tip_frequencies_json,
        root_sequence_json = rules.export.output.root_sequence_json
    output:
        genomes = out_auspice + "/ncov_{build_name}.fasta",
        auspice_json = out_auspice + "/ncov_{build_name}.json",
        tip_frequency_json = out_auspice + "/ncov_{build_name}_tip-frequencies.json",
        root_sequence_json = out_auspice + "/ncov_{build_name}_root-sequence.json"
    log:
        "logs/fix_colorings_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json} 2>&1 | tee {log} &&
        cp {input.frequencies} {output.tip_frequency_json} &&
        cp {input.root_sequence_json} {output.root_sequence_json}
        # Copy genomes to enable genome download feature by default
        cp {input.genomes} {output.genomes} 
        """

checkpoint ncov_clade_assignment:
    message: "Asign clades using nextstrain's asignment script"
    input:
        alignment = rules.combine_samples.output.alignment,
        subclades = rules.rename_subclades.output.ncov_clade,
        clades = config["files"]["clades"]
    output:
       clades = outdir + "/{build_name}/nextstrain_clade_assignment.tsv"
    params:
        clades = outdir + "/{build_name}/ncov_subclades.tsv"
    conda: config["conda_environment"]
    threads: 64
    shell:
        """
        cat {input.clades} {input.subclades} > {params.clades} && \
        python3 scripts/assign_clades.py --nthreads {threads}  \
        --alignment {input.alignment} \
        --clades {params.clades} \
        --chunk-size  48 \
        --output {output.clades}
        """

rule get_mutations:
    message: "Get direct mutations by comparing aligned strain sequence with reference"
    input:
        reference = config["files"]["reference"],
        alignment = rules.aggregate_alignments.output.alignment
    output:
        muts = outdir + "/{build_name}/mutations/tip_all_mutations.tsv",
    log:
        "logs/{build_name}/extract_SNPS.txt"
    params:
        refid = 'Wuhan/Hu-1/2019',
        out = outdir + "/{build_name}/mutations"
    conda: config["conda_environment"]
    shell:
        """
        python scripts/extract_SNPs.py --reference {input.reference} \
        --alignment {input.alignment} \
        -o {params.out} --refid {params.refid} 2>&1 | tee {log}
        """
