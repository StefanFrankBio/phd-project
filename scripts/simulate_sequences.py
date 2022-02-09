import itertools
import random

def generate_reference():
    nucleotides = "TCAG"
    codons = itertools.product(nucleotides, nucleotides, nucleotides)
    codons = ["".join(tpl) for tpl in codons]
    random.shuffle(codons)
    reference = "".join(codons)
    return reference


def mutate_sequence(sequence):
    nucleotides = "TCAG"
    idx = random.randint(0, len(sequence)-1)
    nucleotides = nucleotides.replace(sequence[idx], "")
    mutation = random.choice(nucleotides)
    return sequence[:idx] + mutation + sequence[idx+1:]


def main():
    reference = generate_reference()
    mutation = reference
    with open("input/tests/simulated_reference.fasta", "w") as output:
        print(">simulated_reference", file=output)
        print(reference, file=output)
    with open("input/tests/simulated_mutations.fasta", "w") as output:
        for i in range(100):
            mutation = mutate_sequence(mutation)
            print(f">simulated_mutation_{i}", file=output)
            print(mutation, file=output)


if __name__ == "__main__":
    main()