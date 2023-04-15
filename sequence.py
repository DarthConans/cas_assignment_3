import random
import hashlib
from utils import AMINO_CHARS_LETTERS, translate_codon_sequence_to_aas
from antigenic_calculators import bloom_antigenic_calculator

class Sequence:

    def __init__(self, sequence, base_sequence, original_amino_acids):
        self.__sequence__ = sequence
        self.__mutated_codons__ = set()
        self.__mutated_site_indexes__ = set()
        for i in range(len(base_sequence)):
            new_seq_char = sequence[i]
            base_seq_char = base_sequence[i]
            if new_seq_char != base_seq_char:
                self.__mutated_codons__.add(i)
                self.__mutated_site_indexes__.add(i // 3)
        self.__hash_val__ = hash(sequence)
        self.__base_sequence__ = base_sequence
        self.__original_amino_acids__ = original_amino_acids
        self.__amino_acids__ = translate_codon_sequence_to_aas(self.__sequence__)
        self.__non_neutral_mutations__ = []
        for i in range(len(self.__amino_acids__)):
            new_amino = self.__amino_acids__[i]
            old_amino = original_amino_acids[i]
            if new_amino != old_amino:
                self.__non_neutral_mutations__.append({
                    "site": i,
                    "old_amino": old_amino,
                    "new_amino": new_amino
                })
        self.__bloom_antigenic_fitness__ = None

    def get_bloom_antigenic_fitness(self):
        if self.__bloom_antigenic_fitness__ is None:
            self.__bloom_antigenic_fitness__ = bloom_antigenic_calculator.calculate_fitness_of_sequence(self)
        return self.__bloom_antigenic_fitness__

    def get_non_neutral(self):
        return self.__non_neutral_mutations__

    def get_hash_val(self):
        return self.__hash_val__

    def get_sequence(self):
        return self.__sequence__

    def get_mutated_site_indexes(self):
        return self.__mutated_site_indexes__

    def __str__(self):
        return self.__repr__()

    def __repr__(self) -> str:
        return f"SEQUENCE {self.__sequence__}, CODON INDEXES {self.__mutated_codons__}, SITE INDEXES " \
               f"{self.__mutated_site_indexes__}, ORIGINAL SEQUENCE {self.__base_sequence__}, AAS {self.__amino_acids__}" \
               f", ORIGINAL AAS {self.__original_amino_acids__}, NON NEUTRAL MUTATIONS {self.__non_neutral_mutations__}"

    def __hash__(self) -> int:
        return self.__hash_val__

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Sequence):
            return False
        if self.__hash_val__ != o.__hash_val__:
            return False
        return o.__sequence__ == self.__sequence__

    @classmethod
    def generate_unique_randoms(cls, maximum, number_to_generate, minimum=0):
        randoms = set()

        while len(randoms) < number_to_generate:
            randoms.add(random.randint(minimum, maximum))
        return randoms

    def __mutate_sequence__(self, num_mutations, force_mutations=True):
        amino_indexes = self.generate_unique_randoms(len(self.__sequence__) - 1, num_mutations)
        new_seq = list(self.__sequence__).copy()
        for amino_index in amino_indexes:
            new_seq_char = new_seq[amino_index]
            if force_mutations:
                potential_chars = list(AMINO_CHARS_LETTERS - {new_seq_char})
            else:
                potential_chars = list(AMINO_CHARS_LETTERS)
            new_char = random.choice(potential_chars)
            new_seq[amino_index] = new_char
        return Sequence(''.join(new_seq), self.__base_sequence__, self.__original_amino_acids__)

    def __gen_unique_mutations__(self, num_to_generate, num_mutations, force_mutations=True):
        ret = set()
        while len(ret) < num_to_generate:
            ret.add(self.__mutate_sequence__(num_mutations, force_mutations))
        return ret

    def __gen_mutations__(self, num_to_generate, num_mutations, force_mutations=True):
        ret = []
        for i in range(num_to_generate):
            ret.append(self.__mutate_sequence__(num_mutations, force_mutations))

    def generate_mutations(self, num_to_generate=100, num_mutations=1, unique=True, force_mutations=True):
        if unique:
            return self.__gen_unique_mutations__(num_to_generate, num_mutations, force_mutations)
        else:
            return self.__gen_mutations__(num_to_generate, num_mutations, force_mutations)

    def probabilistic_combine(self, to_combine_with, probability):
        if isinstance(to_combine_with, Sequence):
            other_chars = list(to_combine_with.get_sequence())
        else:
            other_chars = list(to_combine_with)
        new_chars = list(self.get_sequence())
        for i in range(len(new_chars)):
            prob = random.random()
            should_combine = prob > probability
            if should_combine:
                new_chars[i] = other_chars[i]
        return Sequence(''.join(new_chars), self.__base_sequence__, self.__original_amino_acids__)

    def force_combine(self, to_combine_with, number_of_sites_to_swap):
        if isinstance(to_combine_with, Sequence):
            other_chars = list(to_combine_with.get_sequence())
        else:
            other_chars = list(to_combine_with)
        new_chars = list(self.get_sequence())
        sites_to_swap = self.generate_unique_randoms(len(self.__sequence__) - 1, number_of_sites_to_swap)
        for site_to_swap in sites_to_swap:
            new_chars[site_to_swap] = other_chars[site_to_swap]
        return Sequence(''.join(new_chars), self.__base_sequence__, self.__original_amino_acids__)


if __name__ == '__main__':
    orig_seq = Sequence("AAAAAA", "AAAAAA", translate_codon_sequence_to_aas("AAAAAA"))
    t = {orig_seq}
    t2 = orig_seq.generate_mutations(3, 2, True, True)
    for seq in t2:
        print(seq)
    t3 = orig_seq.probabilistic_combine("GGGGGG", .5)
    t4 = orig_seq.force_combine("CCCCCC", 6)
    t5 = t4.get_bloom_antigenic_fitness()
    print("whatever")
