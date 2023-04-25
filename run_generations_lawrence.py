import functools

from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_fitness_set
import os
from numpy import arange
from multiprocessing import Pool
from functools import partial


def mutation_run(run, num_mutations):
    if not os.path.exists(f"results/Lawrence/mutations/{num_mutations}/"):
        os.makedirs(f"results/Lawrence/mutations/{num_mutations}/")
    if not os.path.exists(f"results/Lawrence/mutations/{num_mutations}/{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Lawrence/mutations/{num_mutations}/{run}", neutral,
                                              number_of_generations=100, number_of_children=3,
                                              number_of_mutations=num_mutations)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()
    print(f"FINISHED {run} RUN {num_mutations} NUMBER OF MUTATIONS")


if __name__ == '__main__':
    #antigen_neutral = get_antigenic_one_hop()
    #neutral = get_one_hop()
    #original_aas = load_sequence()

    neutral = get_one_hop()
    for num_mutations in [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
        args = [x for x in range(100)]
        f = partial(mutation_run, num_mutations=num_mutations)

        with Pool(8) as p:
            p.map(f, args)

    for top_prob in [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]:
        for run in range(100):
            if not os.path.exists(f"results/Lawrence/interbreed_top/{top_prob}/"):
                os.makedirs(f"results/Lawrence/interbreed_top/{top_prob}/")
            if not os.path.exists(f"results/Lawrence/interbreed_top/{top_prob}/{run}.pkl"):
                genetic_algorithm = Genetic_Algorithm(f"results/Lawrence/interbreed_top/{top_prob}/{run}", neutral,
                                                      number_of_generations=100,
                                                      number_of_children=3, number_of_mutations=3,
                                                      interbreed_top_prob=top_prob)
                genetic_algorithm.run_ga()
                genetic_algorithm.save_results()
            print(f"FINISHED {run} RUN {top_prob} TOP PROB")
    for random_prob in [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]:
        for run in range(100):
            if not os.path.exists(f"results/Lawrence/interbreed_random/{random_prob}/"):
                os.makedirs(f"results/Lawrence/interbreed_random/{random_prob}/")
            if not os.path.exists(f"results/Lawrence/interbreed_random/{random_prob}/{run}.pkl"):
                genetic_algorithm = Genetic_Algorithm(f"results/Lawrence/interbreed_random/{random_prob}/{run}", neutral,
                                                      number_of_generations=100,
                                                      number_of_children=3, number_of_mutations=3,
                                                      interbreed_random_prob=random_prob)
                genetic_algorithm.run_ga()
                genetic_algorithm.save_results()
            print(f"FINISHED {run} RUN {random_prob} RANDOM PROB")
    print("krewl")