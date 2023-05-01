import functools

from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_fitness_set
import os
from numpy import arange
from multiprocessing import Pool
from functools import partial

known_strains = generate_known_strains()
alpha = Sequence(known_strains[0])
beta = Sequence(known_strains[1])
delta = Sequence(known_strains[2])
omicron = Sequence(known_strains[3])
specific_sequences=[alpha, beta, delta, omicron]


def generate_mutations(self, num_to_generate=100, num_mutations=1, unique=True, force_mutations=True):
    sequence = load_sequence()
    known_strains = generate_known_strains()
    alpha = Sequence(known_strains[0], sequence)
    beta = Sequence(known_strains[1], sequence)
    delta = Sequence(known_strains[2], sequence)
    omicron = Sequence(known_strains[3], sequence)
    specific_sequences = [alpha, beta, delta, omicron]
    return specific_sequences

path = f"results/Lawrence/test/test"
if not os.path.exists(f"{path}.pkl"):

    sequence = Sequence(load_sequence())
    sequence.generate_mutations = generate_mutations
    genetic_algorithm = Genetic_Algorithm(f"{path}", [sequence], top_to_preserve=3,
                                          number_of_generations=3,
                                          number_of_children=3, number_of_mutations=3,
                                          interbreed_random_prob=None, report=False, strains_to_check_for=[alpha],
                                          number_of_interbreed_random=10)
    genetic_algorithm.run_ga()
    genetic_algorithm.save_results()
