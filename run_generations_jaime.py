from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_fitness_set

if __name__ == '__main__':
    neut_1_antigenic_1 = get_fitness_set()

    known_strain_names=["alpha", "beta", "delta", "omicron"]
    [alpha, beta, delta, omicron] = generate_known_strains()
    specific_sequences=[Sequence(alpha), Sequence(beta), Sequence(delta), Sequence(omicron)]

    probabilities = [0.01, 0.1, 0.5]

    for i in range(0,len(specific_sequences)):
        for p in probabilities:
            for run in range(0, 10):
                genetic_algorithm = Genetic_Algorithm(f"results/Jaime/prob_{p}/{known_strain_names[i]}_{run}", neut_1_antigenic_1,
                                                        number_of_children=10, interbreed_specific_sequence_prob=p, antigen_weight=1.0,
                                                        interbreed_specific_sequence=specific_sequences[i], strains_to_check_for=[alpha, beta, delta, omicron])
                genetic_algorithm.run_ga()
                genetic_algorithm.save_results()
    print("krewl")