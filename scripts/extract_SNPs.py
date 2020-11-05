"""
Translate gene regions from nucleotides to amino acids.
"""
import argparse
# from efficient_apriori import apriori
import os, sys
import numpy as np
from pandas.core.common import flatten
from Bio import SeqIO, SeqFeature, Seq, SeqRecord, Phylo
from augur.utils import read_node_data, load_features, write_json, write_VCF_translation, get_json_name
from sys import argv, exit


class MissingNodeError(Exception):
    pass


class MismatchNodeError(Exception):
    pass


def safe_translate(sequence, report_exceptions=False):
    """Returns an amino acid translation of the given nucleotide sequence accounting
    for gaps in the given sequence.

    Optionally, returns a tuple of the translated sequence and whether an
    exception was raised during initial translation.

    >>> safe_translate("ATG")
    'M'
    >>> safe_translate("ATGGT-")
    'MX'
    >>> safe_translate("ATG---")
    'M-'
    >>> safe_translate("ATGTAG")
    'M*'
    >>> safe_translate("")
    ''
    >>> safe_translate("ATGT")
    'MX'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Data.CodonTable import TranslationError
    from Bio.Seq import CodonTable
    translation_exception = False

    # sequences not mod 3 give messy BiopythonWarning, so avoid by padding.
    if len(sequence) % 3:
        sequence_padded = sequence + "N" * (3 - len(sequence) % 3)
    else:
        sequence_padded = sequence
    # try:
    #     # Attempt translation by extracting the sequence according to the
    #     # BioPhython SeqFeature in frame gaps of three will translate as '-'
    #     translated_sequence = str(Seq.Seq(sequence_padded).translate(gap='-'))
    # except TranslationError:
    #     translation_exception = True
    # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
    # one by one.
    codon_table = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
    str_seq = str(sequence_padded)
    codons = np.frombuffer(str_seq[:len(str_seq) - len(str_seq) % 3].encode(), dtype='S3').astype("U")
    assert len(codons) > 0
    aas = []

    for c in codons:
        # Parse result of single codon translation, add amino acids as
        # appropriate.
        try:
            aa = codon_table.get(c)
            if aa is None:
                if c == '---':
                    aas.append('-')
                else:
                    aas.append('X')
            else:
                aas.append(aa)
        except (TranslationError, ValueError):
            aas.append('X')

    translated_sequence = "".join(aas)
    assert (len(translated_sequence) == len(sequence_padded) / 3)
    if report_exceptions:
        return translated_sequence, translation_exception
    else:
        return translated_sequence


def translate_feature(aln, feature):
    '''
    Translates a subsequence of input nucleotide sequences.

    Parameters
    ----------
    aln : dict
        sequences indexed by node name

    feature : Bio.Seq.Seq
        BioPython sequence feature

    Returns
    -------
    dict :
        translated sequences indexed by node name

    '''
    translations = {}
    for sname, seq in aln.items():
        fseq = str(feature.extract(seq))
        aa_seq = safe_translate(fseq)
        translations[sname] = {}
        translations[sname]['aa'] = aa_seq
        translations[sname]['nt'] = fseq

    return translations


def construct_mut(start, pos, end):
    return str(start) + str(pos) + str(end)


def _nt_mut_to_codon_coord(annotation, feature, nt_pos):
        """
        :param annotation: A dictionary of features with starts, ends coordinates and strand
        :param feature: name of feature in the annotation dictionary. It's gene name for example.
        :param zero-based nt_pos: position of the mutation on feature/CDS/gene
        :return:
            an integer indicating coordinate of the mutation on genome
            an integer indicating position of the first base of the code that the muation stand on
        Example: mutations: C7G,A14-,G25N, G34A. Condon starts: C7G:7,
        Result ( different format than the real output): C7G:7,7, A14-:14,13, G25N:25,25, G34A:34,33
        refseq
        012 345 678 901 234 567 890 123 456 7 # 'ruler' for counting base by hand, use zero-base to keep one digit
        act TCT CCC CAG GAC AAC CAA ATG GCA A # feature starts at 'TCT' ends at 'A' - truncated
        0123 456 789 012 345 678 901 2
        aGGG GGT GGA ACA GGA GAG ACC C
        seq1
        012 345 678 901 234 567 890 123 456 7
        act TCT GCC CAG G-C AAC CAA ATG NCA A
        0123 456 789 012 345 678 901 2
        tGGG GAT GGA ACA GGA GAG ACC T
        """
        start = annotation[feature]['start']
        end = annotation[feature]['end']
        strand = annotation[feature]['strand']
        codon_pos = 3 * int(nt_pos / 3)
        if strand == '+':
            return start + nt_pos, int(start + codon_pos)
        else:
            return start - nt_pos, int(start - codon_pos)

def aa_to_nt_pos(annotation, feature, aa_pos):
    """
    :param annotation: A dictionary of features with starts, ends coordinates and strand.
    Coordinates in annotation are one-based count
    :param feature: name of feature in the annotation dictionary. It's gene name for example.
    :param zero-based aa_pos: position on amino acid sequence
    :return: an integer indicating position on genome

    Example: mutations: 2K.P2A, 2K.D4X, 2K.A8X, NS5.G2D. Features: 2K:4..28, NS5:30..51
            Genomic positions for mutations: 2K.P2A:7, 2K.D4X:13, 2K.A8X:25, NS5.G2D:33
                       012345678
    refseq '2K' gene: 'SPQDNQMAX'
    seq1   '2K' gene: 'SAQXNQMXX'
                        01234567
    refseq 'NS5' gene: 'GGGTGETX'
    seq1   'NS5' gene: 'GDGTGETX'

    """
    start = annotation[feature]['start']
    end = annotation[feature]['end']
    strand = annotation[feature]['strand']
    if strand == '+':
        return int(start + (aa_pos) * 3)
    else:
        return int(start - (aa_pos) * 3)


def _find_snps(rseq, nseq, annotations, fname, type = 'dna'):
    snps = {}
    for idx, (a, d) in enumerate(zip(rseq, nseq)):
        # mut = d if a != d else '-'  # Naive check for mutation
        if a != d :
            if type == 'prot' and d != 'X': # give coordinate to the first base of the codon
                pos = aa_to_nt_pos(annotations, fname, idx)
                p = idx + 1
                snps[pos] = {"g": fname, "p": p, "ref": a, "alt": d}
            if type == 'nuc' and a != 'N' and d != 'N': # give coordinate also to the first base of the codon
                p, pos = _nt_mut_to_codon_coord(annotations, fname, idx)
                snps[pos] = {"g": fname, "p": p, "ref": a, "alt": d}
    return snps


def generate_SNPs_table(seq_ids, translations, annotations):
    """
    :param seq_ids: list of sequences to be processed that doesn't contain 'refseq'
    :param refid: Sequence Id of the reference sequence in the alignment
    :param translations:
    :param annotations:
    :param nualn: nucleotide alignment for SNP calculation
    :return:
    """
    seq_ids = [k for k in seq_ids if k != 'refseq'] # sanitizing
    aa_muts = {}
    # fasta input shouldn't have mutations on root, so give empty entry
    aa_muts['refseq'] = {"aa_sequences": {}, "nt_sequences": {}}
    for fname, val in translations.items():
        prot = val['refseq']['aa']
        dna = val['refseq']['nt']
        aa_muts['refseq']["aa_sequences"][fname] = "".join(prot)
        aa_muts['refseq']["nt_sequences"][fname] = "".join(dna)
    for seqid in seq_ids:
        aa_muts[seqid] = {"aa_sequences": {}, "nt_sequences": {}, 'aa': {}, 'nt': {}}
        for fname, val in translations.items():
            if 'refseq' in val and seqid in val:
                rprot = val['refseq']['aa']
                nprot = val[seqid]['aa']
                rdna = val['refseq']['nt']
                ndna = val[seqid]['nt']
                aa_muts[seqid]["aa_sequences"][fname] = "".join(nprot)
                aa_muts[seqid]["nt_sequences"][fname] = "".join(ndna)
                aa_muts[seqid]['aa'].update(_find_snps(rprot, nprot, annotations, fname, 'prot'))
                aa_muts[seqid]['nt'].update(_find_snps(rdna, ndna, annotations, fname, 'nuc'))

            elif 'refseq' not in val[fname] and seqid not in val[fname]:
                print("\n*** Can't find 'refseq' OR %s in the alignment provided!" % (seqid))
                raise MismatchNodeError()
            else:
                print("no sequence pair for nodes refseq-%s" % (seqid))
    return aa_muts


def get_genes_from_file(fname):
    genes = []
    if os.path.isfile(fname):
        with open(fname) as ifile:
            for line in ifile:
                fields = line.strip().split('#')
                if fields[0].strip():
                    genes.append(fields[0].strip())
    else:
        print("File with genes not found. Looking for", fname)

    unique_genes = np.unique(np.array(genes))
    if len(unique_genes) != len(genes):
        print("You have duplicates in your genes file. They are being ignored.")
    print("Read in {} specified genes to translate.".format(len(unique_genes)))

    return unique_genes


# Biopython's trees don't store links to node parents, so we need to build
# a map of each node to its parent.
# Code from the Bio.Phylo cookbook: http://biopython.org/wiki/Phylo_cookbook
def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents


def annotate_phylotree_parents(tree):
    # Get all parent nodes by node.
    parents_by_node = all_parents(tree)

    # Next, annotate each node with its parent.
    for node in tree.find_clades():
        if node == tree.root:
            node.parent = None
        else:
            node.parent = parents_by_node[node]

    # Return the tree.
    return tree


def update_dict(dict, key, val):
    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1
    return dict


def run(args):
    # If genes is a file, read in the genes to translate
    if args.genes and len(args.genes) == 1 and os.path.isfile(args.genes[0]):
        genes = get_genes_from_file(args.genes[0])
    else:
        genes = args.genes
    outdir = os.getcwd()
    if args.outdir:
        outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    ## check file format and read in sequences
    from Bio import AlignIO
    alignment = AlignIO.read(open(args.alignment), "fasta")
    seq_dict = {}
    refid = SeqIO.read(args.reference, 'genbank').id
    if args.refid:
       refid = args.refid 

    for r in alignment:
        if r.id == refid:
            seq_dict['refseq'] = r.seq
        else:
            seq_dict[r.id] = r.seq

    ## load features; only requested features if genes given
    features = load_features(args.reference, genes)
    print("Read in {} features from reference sequence file".format(len(features)))
    if features is None:
        print("ERROR: could not read features of reference sequence file")
        return 1

    ### translate every feature - but not 'nuc'!
    translations = {}
    for fname, feat in features.items():
        if feat.type != 'source':
            translations[fname] = translate_feature(seq_dict, feat)

    ## glob the annotations for later auspice export
    #
    # Note that BioPython FeatureLocations use
    # "Pythonic" coordinates: [zero-origin, half-open)
    # Starting with augur v6 we use GFF coordinates: [one-origin, inclusive]
    annotations = {}
    for fname, feat in features.items():
        annotations[fname] = {'seqid': args.reference,
                              'type': feat.type,
                              'start': int(feat.location.start) + 1,
                              'end': int(feat.location.end),
                              'strand': '+' if feat.location.strand else '-'}

    # for f in annotations:
    #     print(f"{f}\t{annotations[f]['start']}\t{annotations[f]['end']}\t{annotations[f]['strand']}")

    ## determine amino acid mutations for each sequence
    seq_ids = [k for k in seq_dict.keys() if k != 'refseq']
    aa_muts = generate_SNPs_table(seq_ids, translations, annotations)

    print("Writing SNP data ...")
    fwa = open(outdir + "/tip_aa_snps.tsv", "w")
    fwn = open(outdir + "/tip_all_mutations.tsv", "w")
    fwa.write(f'Node\tMutations\tPath\n')
    fwn.write(f'Strain\tPosition\tNuRef\tNuAlt\tNuMut\tGene\tAaPos\tAaRef\tAaAlt\tAaMut\n')

    fwae = open(outdir + "/tip_aa_snps_ex.tsv", "w")
    # fwne = open(outdir + "/tip_nt_snps_ex.tsv", "w")
    fwae.write(f'Node\tMutations\n')
    # fwne.write(f'Node\tMutations\n')

    for seqid in seq_ids:
        amuts = aa_muts[seqid]['aa']
        tmp = [f'{amuts[pos]["g"]}.{amuts[pos]["ref"]}{amuts[pos]["p"]}{amuts[pos]["alt"]}' for pos in amuts]
        fwa.write(f'{seqid}\t{",".join(tmp)}\n')
        nmuts = aa_muts[seqid]['nt'] # mutations within gene features
        gmuts = {}
        for pos, (a, d) in enumerate(zip(seq_dict['refseq'], seq_dict[seqid])):
            if a != d and a != 'N' and d != 'N' and (pos + 1) not in nmuts: # capture mutations outside gene features
                gmuts[pos + 1] = {'g': 'no_gene', 'ref': a, 'alt': d, 'p': pos + 1}
        for pos, mut in nmuts.items():
            gene = nmuts[pos]['g']
            n_ref = mut['ref']
            n_alt = mut['alt']
            n_mut =f'{n_ref}{mut["p"]}{n_alt}'
            # fwne.write(f'{seqid}\t{n_mut}\n')
            if pos in amuts:
                gene = amuts[pos]['g']
                a_ref = amuts[pos]['ref']
                a_alt = amuts[pos]['alt']
                a_pos = amuts[pos]['p']
                a_mut = f'{gene}.{a_ref}{amuts[pos]["p"]}{a_alt}'
                fwae.write(f'{seqid}\t{a_mut}\n')
            else:
               a_ref = a_alt = a_pos = a_mut = '='
            fwn.write(f'{seqid}\t{mut["p"]}\t{n_ref}\t{n_alt}\t{n_mut}\t{gene}\t{a_pos}\t{a_ref}\t{a_alt}\t{a_mut}\n')
        for pos, mut in gmuts.items():
            gene = gmuts[pos]['g']
            a_ref = a_alt = a_pos = a_mut = 'unk'
            n_ref = mut['ref']
            n_alt = mut['alt']
            n_mut =f'{mut["ref"]}{mut["p"]}{mut["alt"]}'
            # fwne.write(f'{seqid}\t{n_mut}\n')
            fwn.write(f'{seqid}\t{mut["p"]}\t{n_ref}\t{n_alt}\t{n_mut}\t{gene}\t{a_pos}\t{a_ref}\t{a_alt}\t{a_mut}\n')
    fwa.close()
    fwn.close()
    fwae.close()
    # fwne.close()
    print("Done!")
    print(f"Output are: {outdir}/tip_all_mutations.tsv")


def main():
    parser = argparse.ArgumentParser(
        prog="extract_SNPs.py",
        description="Report direct mutations based on alignments")

    parser.add_argument('--refid', help="Reference sequence ID in the alignment file in case there is a mismatch between reference sequence and alignment")
    parser.add_argument('-r', '--reference', required=True,
                        help='GenBank file containing the annotation')
    parser.add_argument('-a', '--alignment', help="alignment in fasta format")
    parser.add_argument('-o', '--outdir', help="Output directory for saving files", type=str)
    parser.add_argument('--genes', nargs='+', help="genes to translate (list or file containing list)")

    return run(parser.parse_args())


# Run when called as `python -m augur`, here for good measure.
if __name__ == "__main__":
    exit(main())
