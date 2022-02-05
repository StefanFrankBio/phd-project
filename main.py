import itertools
from Bio import SeqIO, Align

def read_reference(filepath):
    with open(filepath) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record.seq


def read_mutated_sequences(filepath):
    mutated_sequences = []
    with open(filepath) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            mutated_sequences.append(record.seq)
    return mutated_sequences


def split_codons(sequence):
    overhang = len(sequence) % 3
    if overhang != 0:
        sequence = sequence[:-overhang]
    return [sequence[i:i+3] for i in range(0, len(sequence), 3)]


def single_substitutions(codon):
    substitutions = []
    for i, nt in enumerate(codon):
        possible_subs = ["A", "C", "G", "T"]
        possible_subs.remove(nt)
        for sub in possible_subs:
            substitutions.append(codon[:i] + sub + codon[i+1:])
    return substitutions


def build_trans_table():
    nucs = "TCAG"
    aminos = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = (itertools.product(nucs, nucs, nucs))
    codons = ["".join(tpl) for tpl in codons]
    return dict(zip(codons, aminos))


def count_synonymous_sites(bool_list):
    return [sum(bool_list[x*3:(x+1)*3])/3 for x in range(3)]


def align_pair(reference, mutated):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -1
    alignment = aligner.align(reference, mutated)
    alignment = str(alignment[0])[:-1].split("\n")
    return alignment


def main():
    trans_table = build_trans_table()
    synonymous_sites = []
    reference = read_reference("simulated_reference.fasta")
    codons = split_codons(reference)
    for codon in codons:
        sns = single_substitutions(codon)
        translations = [trans_table[s] for s in sns]
        is_syn = [amino == trans_table[codon] for amino in translations]
        synonymous_sites += count_synonymous_sites(is_syn)

    mutated_sequences = read_mutated_sequences("simulated_mutations.fasta")
    for seq in mutated_sequences[:1]:
        alignment = align_pair(reference, seq)


if __name__ == "__main__":
    main()
