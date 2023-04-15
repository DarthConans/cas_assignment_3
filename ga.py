from sequence import Sequence
from utils import load_sequence
from multiprocessing import Pool

class Genetic_Algorithm:

    def __init__(self, initial_sequences, number_of_generations=100, top_to_preserve=100, number_of_children=100,
                 number_of_mutations=3, antigen_weight=.5, fitness_weight=.5, force_mutations=True, unique=True):
        self.initial_sequences = initial_sequences
        self.number_of_generations = number_of_generations
        self.top_to_preserve = top_to_preserve
        self.number_of_children = number_of_children
        self.number_of_mutations = number_of_mutations
        self.antigen_weight = antigen_weight
        self.fitness_weight = fitness_weight
        self.force_mutations = force_mutations
        self.unique = unique
        self.generations = [initial_sequences]

    def generate_generation(self, parents):
        return generate_generation(parents, self.number_of_children, self.number_of_mutations,
                                   self.force_mutations, self.unique, self.antigen_weight, self.fitness_weight,
                                   self.top_to_preserve)

    def run_ga(self):
        generation = self.generate_generation(self.initial_sequences)
        self.generations.append(generation)
        for i in range(1, self.number_of_generations):
            generation = self.generate_generation(generation)
            self.generations.append(generation)
        print("krewl")
        print(f"AFTER {self.generations} GENERATIONS, I ACHIEVED A FITNESS OF {generation[0].generate_fitness_score_non_neutral(self.antigen_weight, self.fitness_weight)}")

    def __run_generation__(self):
        pass


def generate_generation(parents, number_of_children=100, number_of_mutations=3, force_mutations=True,
                        unique=True, antigen_weight=.5,
                        fitness_weight=.5, top_to_preserve=100):
    candidates = [parents[0]]
    for sequence in parents:
        candidates.extend(sequence.generate_mutations(number_of_children,
                                                      num_mutations=number_of_mutations,
                                                      force_mutations=force_mutations, unique=unique))
    candidates = sorted(candidates, key=lambda x: x.generate_fitness_score_non_neutral(antigen_weight,
                                                                                fitness_weight),
                 reverse=True)
    ret = candidates[0:top_to_preserve]
    return ret


if __name__ == '__main__':
    original_aas = load_sequence()
    original_sequence = Sequence(original_aas, original_aas)
    first_generation = generate_generation([original_sequence], number_of_children=200)
    genetic_algorithm = Genetic_Algorithm(first_generation, number_of_children=3)
    genetic_algorithm.run_ga()
    print("krewl")
