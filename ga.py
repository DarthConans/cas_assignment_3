from sequence import Sequence, get_antigenic_one_hop, get_one_hop
from utils import load_sequence, AMINO_CHARS_LETTERS
from multiprocessing import Pool
import random
import pickle as pkl

class Genetic_Algorithm:

    valid_preserve_lowest_strategies = [
        None, "rand", "lowest"
    ]

    def __init__(self, out_file, initial_sequences, number_of_generations=100, top_to_preserve=100, number_of_children=100,
                 number_of_mutations=3, antigen_weight=.5, fitness_weight=None, force_mutations=True, unique=True,
                 interbreed_random_prob=None, interbreed_specific_sequence_prob=None, interbreed_specific_sequence=None,
                 preserve_lowest_strategy=None, interbreed_top_prob=None, strains_to_check_for=None):
        self.initial_sequences = initial_sequences
        self.number_of_generations = number_of_generations
        self.top_to_preserve = top_to_preserve
        self.number_of_children = number_of_children
        self.number_of_mutations = number_of_mutations
        self.antigen_weight = antigen_weight
        self.fitness_weight = fitness_weight
        self.force_mutations = force_mutations
        self.unique = unique
        self.fitnesses = [self.generate_fitness_tuple(initial_sequences[0])]
        self.interbreed_random_prob = interbreed_random_prob
        if interbreed_specific_sequence_prob is not None and interbreed_specific_sequence is None:
            raise ValueError("SPECIFIC SEQUENCE PROBABILITY SPECIFIED, BUT NO SEQUENCE SUPPLIED")
        if interbreed_specific_sequence is not None and interbreed_specific_sequence_prob is None:
            raise ValueError("SPECIFIC SEQUENCE SPECIFIED, BUT NO PROBABILITY SUPPLIED")
        self.interbreed_specific_sequence_prob = interbreed_specific_sequence_prob
        self.interbreed_specific_sequence = interbreed_specific_sequence
        if preserve_lowest_strategy not in self.valid_preserve_lowest_strategies:
            raise ValueError(f"{preserve_lowest_strategy} IS NOT A VALID STRATEGY")
        self.preserve_lowest_strategy = preserve_lowest_strategy
        self.interbreed_top_prob = interbreed_top_prob
        self.strains_to_check_for = None if strains_to_check_for is None else {strain_to_check_for: [] for strain_to_check_for in strains_to_check_for}
        self.out_file = out_file


    @classmethod
    def generate_random_base_pairs(cls, length_of_sequence):
        ret = []
        candidates = list(AMINO_CHARS_LETTERS)
        for i in range(length_of_sequence):
            ret.append(random.choice(candidates))
        return "".join(ret)

    def generate_fitness_tuple(self, sequence):
        return sequence.get_bloom_antigenic_fitness(), sequence.get_bloom_fitness()

    def generate_generation(self, parents):
        return generate_generation(parents, self.number_of_children, self.number_of_mutations,
                                   self.force_mutations, self.unique, self.antigen_weight, self.fitness_weight,
                                   self.top_to_preserve, self.interbreed_random_prob,
                                   self.interbreed_specific_sequence_prob, self.interbreed_specific_sequence,
                                   self.interbreed_top_prob, self.preserve_lowest_strategy, self.strains_to_check_for)

    def run_ga(self):
        generation, strains_found_in_generation = self.generate_generation(self.initial_sequences)
        print("DONE WITH ROUND 1")
        for strain_found in strains_found_in_generation:
            self.strains_to_check_for[strain_found].append(1)
        for i in range(2, self.number_of_generations + 1):
            print(f"DONE WITH ROUND {i}")
            generation, strains_found_in_generation = self.generate_generation(generation)
            self.fitnesses.append(self.generate_fitness_tuple(generation[0]))
            for strain_found in strains_found_in_generation:
                self.strains_to_check_for[strain_found].append(i)
        print(
            f"AFTER {self.number_of_generations} GENERATIONS, I ACHIEVED AN ANTIGENTIC FITNESS OF {generation[0].get_bloom_antigenic_fitness()} AND A REGULAR FITNESS OF {generation[0].get_bloom_fitness()}")
        print("krewl")

    def __run_generation__(self):
        pass

    def save_results(self):
        with open(f"{self.out_file}.pkl", "wb") as f:
            pkl.dump(self, f)


def generate_generation(parents, number_of_children=100, number_of_mutations=3, force_mutations=True,
                        unique=True, antigen_weight=.5,
                        fitness_weight=None, top_to_preserve=100, interbreed_random_prob=None,
                        interbreed_specific_prob=None, interbreed_specific_sequence=None, interbreed_top_prob=None,
                        preserve_lowest_strategy=None, strains_to_check_for=None):
    candidates = [parents[0]]
    for sequence in parents:
        candidates.extend(sequence.generate_mutations(number_of_children,
                                                      num_mutations=number_of_mutations,
                                                      force_mutations=force_mutations, unique=unique))
    if unique:
        candidates = set(candidates)
        while len(candidates) < top_to_preserve:
            for sequence in parents:
                candidates.update(sequence.generate_mutations(number_of_children,
                                                              num_mutations=number_of_mutations,
                                                              force_mutations=force_mutations, unique=unique))
                if len(candidates) >= top_to_preserve:
                    break
        candidates = list(candidates)
    if interbreed_random_prob:
        should_breed = random.random() < interbreed_random_prob
        if should_breed:
            random_sequence = Genetic_Algorithm.generate_random_base_pairs(len(candidates[0].get_sequence()))
            to_breed_with = candidates[0]
            child = to_breed_with.force_combine(random_sequence, number_of_mutations)
            candidates.append(child)
    if interbreed_specific_prob:
        should_breed = random.random() < interbreed_random_prob
        if should_breed:
            to_breed_with = candidates[0]
            child = to_breed_with.force_combine(interbreed_specific_sequence, number_of_mutations)
            candidates.append(child)

    if fitness_weight is None:
        cand_fitnesses = sorted(candidates, key=lambda x: x.get_bloom_fitness(), reverse=True)
        max_fitness = cand_fitnesses[0].get_bloom_fitness()
        fitness_weight = 1 / max_fitness
    candidates = sorted(candidates, key=lambda x: x.generate_fitness_score_non_neutral(antigen_weight,
                                                                                       fitness_weight),
                        reverse=True)
    if interbreed_top_prob:
        should_breed = random.random() < interbreed_top_prob
        if should_breed:
            parent_1 = candidates[0]
            parent_2 = candidates[1]
            child = parent_1.force_combine(parent_2, number_of_mutations)
            candidates.append(child)
        candidates = sorted(candidates, key=lambda x: x.generate_fitness_score_non_neutral(antigen_weight,
                                                                                           fitness_weight),
                            reverse=True)
    number_to_preserve = top_to_preserve if not preserve_lowest_strategy else top_to_preserve - 1
    ret = candidates[0:number_to_preserve]
    strains_found = []
    if strains_to_check_for:
        for strain_to_check_for in strains_to_check_for.keys():
            if strain_to_check_for in candidates:
                strains_found.append(strain_to_check_for)
    if preserve_lowest_strategy:
        to_preserve = None
        if preserve_lowest_strategy == "rand":
            to_preserve = random.choice(candidates[number_to_preserve:])
        elif preserve_lowest_strategy == "lowest":
            to_preserve = candidates[-1]
        if to_preserve is None:
            raise ValueError(f"{preserve_lowest_strategy} IS NOT A VALID STRATEGY")
        ret.append(to_preserve)
    return ret, strains_found


if __name__ == '__main__':
    antigen_neutral = get_antigenic_one_hop()
    neutral = get_one_hop()
    original_aas = load_sequence()
    original_sequence = Sequence(original_aas)
    specific_sequence = Sequence("AAUAUUACAAACUUGUGCCCUUUUGGUGAAGUUUUUAACGCCACCAGAUUUGCAUCUGUUUAUGCUUGGAACAGGAAGAGAAUCAGCAACUGUGUUGCUGAUUAUUCUGUCCUAUAUAAUUCCGCAUCAUUUUCCACUUUUAAGUGUUAUGGAGUGUCUCCUACUAAAUUAAAUGAUCUCUGCUUUACUAAUGUCUAUGCAGAUUCAUUUGUAAUUAGAGGUGAUGAAGUCAGACAAAUCGCUCCAGGCCAAACUGGAAAGAUUGCUGAUUAUAAUUAUAAAUUACCAGAUGAUUUUACAGGCUGCGUUAUAGCUUGGAAUUCUAACAAUCUUGAUUCUAAGGUUGGUGGUAAUUAUAAUUACCUGUAUAGAUUGUUUAGGAAGUCUAAUCUCAAACCUUUUGAGAGAGAUAUUUCAACUGAAAUCUAUCAGGCCGGUAGCACACCUUGUAAUGGUGUUCAAGGUUUUAAUUGUUACUUUCCCUUACAAUCAUAUGGUUUCCAACCCACUAAUGGUGUUGGUUACCAACCAUACAGAGUAGUAGUACUUUCUUUUGAACUUCUACAUGCACCAGCAACU")
    genetic_algorithm = Genetic_Algorithm("test", neutral, number_of_generations=3, number_of_children=3,
                                          interbreed_random_prob=.7, interbreed_specific_sequence_prob=.7,
                                          fitness_weight=None, interbreed_specific_sequence=specific_sequence,
                                          interbreed_top_prob=.7, preserve_lowest_strategy="lowest",
                                          strains_to_check_for=[original_sequence])
    genetic_algorithm.run_ga()
    genetic_algorithm.save_results()
    print("krewl")
