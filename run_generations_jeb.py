import functools
import os
from multiprocessing import Pool
from functools import partial

from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_fitness_set

def run_ga_parallel(run, start_name, start_set, known_strains, num_to_gen=100):

    print(f"Running round: {run}")
    num_mutations = 3
    if start_name=="neutral_1_mutation_no_stop":
        num_mutations = 1
    genetic_algorithm = Genetic_Algorithm(f"results/Jeb/{start_name}/{start_name}_{run}", start_set, number_of_generations=1000,
                                          number_of_children=10, top_to_preserve=num_to_gen,number_of_mutations=num_mutations,
                                          interbreed_random_prob=None, interbreed_specific_sequence_prob=None,
                                          fitness_weight=None, antigen_weight=1, interbreed_specific_sequence=None,
                                          interbreed_top_prob=None, preserve_lowest_strategy=None,
                                          strains_to_check_for=known_strains, report=True)
    genetic_algorithm.run_ga()
    genetic_algorithm.save_results()

if __name__ == '__main__':
    antigen_neutral = get_antigenic_one_hop()
    neutral = get_one_hop()
    original_aas = load_sequence()
    original_sequence = Sequence(original_aas)
    neut_1_antigenic_1 = get_fitness_set()

    num_to_gen = 100
    known_strains = [Sequence(x) for x  in generate_known_strains()]

    #starts = [antigen_neutral, neutral, neut_1_antigenic_1, neutral]
    #start_names = ["antigen_neutral", "neutral", "neutral_1_hops_antigenic", "neutral_1_mutation"]
    starts = [ neutral]
    start_names = ["neutral_1_mutation_no_stop"]
    for name in start_names:
        if not os.path.exists(f"results/Jeb/{name}/"):
            os.makedirs(f"results/Jeb/{name}/")
    for i, s in enumerate(starts):
        print("Running start: ", start_names[i])
        f = partial(run_ga_parallel, start_name=start_names[i], start_set=s, known_strains=known_strains)
        f(0)
    print("krewl")
