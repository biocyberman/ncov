import copy
from os import environ
from os import getcwd
from os import path
from pathlib import Path
from snakemake.logging import logger
from snakemake.utils import validate


# Get output paths and convert them to relative paths
outdir = config['outdir'] if "outdir" in config else "mutations"
config['workdir'] = outdir + '/ncov-aau'
outdir = path.relpath(outdir, getcwd())

merged_data = Path(outdir).joinpath("data")

configfile: "defaults/parameters.yaml"

validate(config, schema=config['workdir'] + "/workflow/schemas/config.schema.yaml")

# default build if none specified in config
if "builds" not in config:
    config["builds"] = {
        "global": {
            "subsampling_scheme": "region_global",
        }
    }

# This is customization for Denmark. The full build names below (i.e. Denmark_Light, Denmark_Full, Test, etc) must match configuration in my_profiles/denmark/builds.yml
short_build_names = {'light': 'Denmark_Light', 'full': 'Denmark_Full', 'global': 'Global', 'test': 'Test'}
BUILD_NAMES = [short_build_names[k] for k in config['custom_build'].split(',')] if "custom_build" in config else list(config["builds"].keys())

# Define patterns we expect for wildcards.
wildcard_constraints:
    # Allow build names to contain alpha characters, underscores, and hyphens
    # but not special strings used for Nextstrain builds.
    build_name = r'(?:[_a-zA-Z-](?!(tip-frequencies|gisaid|zh)))+',
    date = r"[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

# from pprint import PrettyPrinter as pp
# pp(indent=2).pprint(config['builds'])

# Create a standard ncov build for auspice, by default.
rule all:
    input:
        direct_mutaions  = outdir + "/all_mutations.tsv"

rule clean:
    message: "Removing directories: {params}"
    params:
        outdir
    shell:
        "rm -rfv {params}"



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
        date = numeric_date(date.today())
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


rule get_mutations:
    message: "Get direct mutations by comparing aligned strain sequence with reference"
    input:
        reference = config["files"]["reference"],
        alignment = config["alignment"] if (config.get("alignment") and path.exists(config.get("alignment")))  else rules.aggregate_alignments.output.alignment
    output:
        muts = outdir + "/all_mutations.tsv",
    log:
        "logs/extract_SNPS.txt"
    params:
        refid = config.get("reference_id", 'Wuhan/Hu-1/2019'),
        variant_list = f'--variant-list {config["variant_list"]}' if (config.get("variant_list") is not None)  else ""
    conda: config["conda_environment"]
    shell:
        """
        augur extract-snps --reference {input.reference} \
        --alignment {input.alignment} \
        {params.variant_list} \
        -o {outdir} --refid {params.refid} 2>&1 | tee {log}
        """

