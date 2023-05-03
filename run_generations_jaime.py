import os
from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_one_hop, get_fitness_set, get_antigenic_one_hop
from multiprocessing import Pool
from functools import partial

neutral = get_one_hop()
antigenic_neutral = get_antigenic_one_hop()
neut_1_antigenic_1 = get_fitness_set()

known_strains = generate_known_strains()
alpha = Sequence(known_strains[0])
beta = Sequence(known_strains[1])
delta = Sequence(known_strains[2])
omicron = Sequence(known_strains[3])
specific_sequences=[alpha, beta, delta, omicron]

def alpha_run(arg_tuple):
    prob = arg_tuple[0]
    run = arg_tuple[1]
    if not os.path.exists(f"results/Jaime/prob_{prob}/alpha_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_{prob}/alpha_{run}", neut_1_antigenic_1,
                                            number_of_children=10,
                                            antigen_weight=1.0,
                                            interbreed_specific_sequence=specific_sequences[0],
                                            interbreed_specific_sequence_prob=prob,
                                            strains_to_check_for=specific_sequences)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()

def beta_run(arg_tuple):
    prob = arg_tuple[0]
    run = arg_tuple[1]
    if not os.path.exists(f"results/Jaime/prob_{prob}/beta_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_{prob}/beta_{run}", neut_1_antigenic_1,
                                            number_of_children=10,
                                            antigen_weight=1.0,
                                            interbreed_specific_sequence=specific_sequences[1],
                                            interbreed_specific_sequence_prob=prob,
                                            strains_to_check_for=specific_sequences)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()

def delta_run(arg_tuple):
    prob = arg_tuple[0]
    run = arg_tuple[1]
    if not os.path.exists(f"results/Jaime/prob_{prob}/delta_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_{prob}/delta_{run}", neut_1_antigenic_1,
                                            number_of_children=10,
                                            antigen_weight=1.0,
                                            interbreed_specific_sequence=specific_sequences[2],
                                            interbreed_specific_sequence_prob=prob,
                                            strains_to_check_for=specific_sequences)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()


def omicron_run(arg_tuple):
    prob = arg_tuple[0]
    run = arg_tuple[1]
    if not os.path.exists(f"results/Jaime/prob_{prob}/omicron_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_{prob}/omicron_{run}",
                                                neut_1_antigenic_1,
                                                number_of_children=10,
                                                antigen_weight=1.0,
                                                interbreed_specific_sequence=specific_sequences[3],
                                                interbreed_specific_sequence_prob=prob,
                                                strains_to_check_for=specific_sequences)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()

def no_mix_neutral(run):
    if not os.path.exists(f"results/Jaime/prob_0/neutral_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_0/neutral_{run}",
                                                neutral,
                                                number_of_children=10,
                                                antigen_weight=1.0,
                                                strains_to_check_for=specific_sequences)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()

def no_mix_antigenic_neutral(run):
    if not os.path.exists(f"results/Jaime/prob_0/antigenic_neutral_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_0/antigenic_neutral_{run}",
                                                antigenic_neutral,
                                                number_of_children=10,
                                                antigen_weight=1.0,
                                                strains_to_check_for=specific_sequences)
        genetic_algorithm.run_ga()
        genetic_algorithm.save_results()

if __name__ == '__main__':
    probabilities = [0.01, 0.1, 0.5]

    runs = [x for x in range(0, 10)]
    # alpha
    args0 = []
    for p in probabilities:
        for run in range(0, 10):
            args0.append((p, run))
    with Pool(7) as p:
        p.map(alpha_run, args0)
    # beta
    args1 = []
    for p in probabilities:
        for run in range(0, 10):
            args1.append((p, run))
    with Pool(7) as p:
        p.map(beta_run, args1)
    # delta
    args2 = []
    for p in probabilities:
        for run in range(0, 10):
            args2.append((p, run))
    with Pool(5) as p:
        p.map(delta_run, args2)
    # omicron
    args3 = []
    for p in probabilities:
        for run in range(0, 10):
            args3.append((p, run))
    with Pool(4) as p:
        p.map(omicron_run, args3)
    # no mix neutral
    with Pool(6) as p:
        p.map(no_mix_neutral, runs)
    # no mix antigenically neutral
    with Pool(6) as p:
        p.map(no_mix_antigenic_neutral, runs)

    print("krewl")