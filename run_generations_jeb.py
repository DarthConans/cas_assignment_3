from ga import *
from utils import generate_known_strains
from sequence import Sequence, get_fitness_set

if __name__ == '__main__':
    antigen_neutral = get_antigenic_one_hop()
    neutral = get_one_hop()
    original_aas = load_sequence()
    original_sequence = Sequence(original_aas)
    neut_1_antigenic_1 = get_fitness_set()

    num_to_gen = 100
    random_orig_seq_start = list(original_sequence.generate_mutations(num_to_gen,5,unique=True,force_mutations=True))
    #print(random_orig_seq_start)
    new_start_list_just_seq = []
    for x in random_orig_seq_start:
        new_start_list_just_seq.append(x.__sequence__)
    known_strains = generate_known_strains()


    #get fitness set for project 2 
    starts = [antigen_neutral, neutral, new_start_list_just_seq, neut_1_antigenic_1]
    start_names = ["antigen_neutral","neutral","random", "neutral_1_hops_antigenic"]
    for i,s in enumerate(starts):
        print("Running start: ",start_names[i])
        for run in range(0, 10):
            print("Running round: ", run)
            genetic_algorithm = Genetic_Algorithm(f"results/Jeb/{start_names[i]}/{start_names[i]}_{run}", s,number_of_generations=1000, number_of_children=10,top_to_preserve=num_to_gen,
                                      interbreed_random_prob=None, interbreed_specific_sequence_prob=None,
                                      fitness_weight=None,antigen_weight=1, interbreed_specific_sequence=None,
                                      interbreed_top_prob=None, preserve_lowest_strategy=None,
                                      strains_to_check_for=known_strains)
            genetic_algorithm.run_ga()
            genetic_algorithm.save_results()
    print("krewl")