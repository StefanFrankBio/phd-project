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
    import itertools
    nucs = "TCAG"
    aminos = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = (itertools.product(nucs, nucs, nucs))
    codons = ["".join(tpl) for tpl in codons]
    return dict(zip(codons, aminos))


def count_synonymous_sites(bool_list):
    return [sum(bool_list[x*3:(x+1)*3])/3 for x in range(3)]

def main():
    trans_table = build_trans_table()
    synonymous_sites = []
    with open("geneS.fasta") as f:
        next(f)
        sequence = f.read().replace("\n", "")
    codons = split_codons(sequence)
    for codon in codons:
        sns = single_substitutions(codon)
        translations = [trans_table[s] for s in sns]
        is_syn = [amino == trans_table[codon] for amino in translations]
        synonymous_sites += count_synonymous_sites(is_syn)

if __name__ == "__main__":
    main()
