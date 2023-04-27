import os
from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_fitness_set
from multiprocessing import Pool
from functools import partial

neut_1_antigenic_1 = get_fitness_set()

[alpha, beta, delta, omicron] = generate_known_strains()
specific_sequences=[Sequence(delta), Sequence(omicron)]

def delta_fifty(run):
    if not os.path.exists(f"results/Jaime/prob_0.5/delta_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_0.5/delta_{run}", neut_1_antigenic_1,
                                          number_of_children=10, interbreed_specific_sequence_prob=0.5,
                                          antigen_weight=1.0,
                                          interbreed_specific_sequence=specific_sequences[0],
                                          strains_to_check_for=[alpha, beta, delta, omicron])
    genetic_algorithm.run_ga()
    genetic_algorithm.save_results()


def omicron_run(arg_tuple):
    prob = arg_tuple[0]
    run = arg_tuple[1]
    if not os.path.exists(f"results/Jaime/prob_{p}/omicron_{run}.pkl"):
        genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_{p}/omicron_{run}", neut_1_antigenic_1,
                                          number_of_children=10, interbreed_specific_sequence_prob=prob,
                                          antigen_weight=1.0,
                                          interbreed_specific_sequence=specific_sequences[1],
                                          strains_to_check_for=[alpha, beta, delta, omicron])
    genetic_algorithm.run_ga()
    genetic_algorithm.save_results()


if __name__ == '__main__':
    neut_1_antigenic_1 = get_fitness_set()

    #known_strain_names=["alpha", "beta", "delta", "omicron"]
    #[alpha, beta, delta, omicron] = generate_known_strains()
    #specific_sequences=[Sequence(alpha), Sequence(beta), Sequence(delta), Sequence(omicron)]
    #specific_sequences=[Sequence(omicron)]

    probabilities = [0.01, 0.1, 0.5]

    runs = [x for x in range(10)]
    with Pool(7) as p:
        p.map(delta_fifty, runs)
    args = []
    for p in probabilities:
        for run in range(0, 10):
            args.append((p, run))
    with Pool(7) as p:
        p.map(omicron_run, args)
    print("krewl")