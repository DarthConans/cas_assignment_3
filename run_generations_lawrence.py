import functools

from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_fitness_set
import os
from numpy import arange
from multiprocessing import Pool
from functools import partial

def antigen_weight_run(run):
    if not os.path.exists(f"results/Lawrence/antigen/{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Lawrence/antigen/{run}", neutral,
                                              number_of_generations=100, number_of_children=3, antigen_weight=1,
                                              fitness_weight=0,
                                              number_of_mutations=1, report=False)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()
    print(f"FINISHED ANTIGEN RUN {run}")
def mutation_run(run, num_mutations):
    if not os.path.exists(f"results/Lawrence/mutations/{num_mutations}/{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Lawrence/mutations/{num_mutations}/{run}", neutral,
                                              number_of_generations=100, number_of_children=3, antigen_weight=1,
                                              number_of_mutations=num_mutations, report=False)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()
    print(f"FINISHED {run} RUN {num_mutations} NUMBER OF MUTATIONS")


def top_prob_run(run, top_prob):
    path = f"results/Lawrence/interbreed_top/{str(top_prob).replace('.','_')}/{run}"
    if not os.path.exists(f"{path}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"{path}", neutral,
                                              number_of_generations=100, antigen_weight=1,
                                              number_of_children=3, number_of_mutations=1,
                                              interbreed_top_prob=top_prob, report=False)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()
    print(f"FINISHED {run} RUN {top_prob} TOP PROB")


def random_prob_run(run, random_prob):
    path = f"results/Lawrence/interbreed_random/{str(random_prob).replace('.','_')}/{run}"
    if not os.path.exists(f"{path}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"{path}", neutral,
                                              number_of_generations=100, antigen_weight=1,
                                              number_of_children=3, number_of_mutations=1,
                                              interbreed_random_prob=random_prob, report=False,
                                              number_of_interbreed_random=10)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()
    print(f"FINISHED {run} RUN {random_prob} RANDOM PROB")


if __name__ == '__main__':
    #antigen_neutral = get_antigenic_one_hop()
    #neutral = get_one_hop()
    #original_aas = load_sequence()

    neutral = get_one_hop()
    args = [x for x in range(100)]
    with Pool(19) as p:
        p.map(antigen_weight_run, args)
    num_mutations_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    
    for num_mutations in num_mutations_list:
        if not os.path.exists(f"results/Lawrence/mutations/{num_mutations}/"):
            os.makedirs(f"results/Lawrence/mutations/{num_mutations}/")
    for num_mutations in num_mutations_list:
        args = [x for x in range(100)]
        f = partial(mutation_run, num_mutations=num_mutations)
        with Pool(19) as p:
            p.map(f, args)
    top_probs = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
    for top_prob in top_probs:
        if not os.path.exists(f"results/Lawrence/interbreed_top/{str(top_prob).replace('.','_')}/"):
            os.makedirs(f"results/Lawrence/interbreed_top/{str(top_prob).replace('.','_')}/")
    print(f"UP TO TOP PROB")
    for top_prob in top_probs:
        args = [x for x in range(100)]
        f = partial(top_prob_run, top_prob=top_prob)
        with Pool(19) as p:
            p.map(f, args)
    random_probs = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
    print(f"UP TO RANDOM PROB")
    for random_prob in random_probs:
        if not os.path.exists(f"results/Lawrence/interbreed_random/{str(random_prob).replace('.', '_')}/"):
            os.makedirs(f"results/Lawrence/interbreed_random/{str(random_prob).replace('.', '_')}/")
    for random_prob in random_probs:
        args = [x for x in range(100)]
        f = partial(random_prob_run, random_prob=random_prob)
        with Pool(19) as p:
            p.map(f, args)
    print("krewl")