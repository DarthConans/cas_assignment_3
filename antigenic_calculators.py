from bindingcalculator import *
from os.path import isfile
import pickle as pkl


class bloom_antigenic_calculator():
    if isfile("caches/bloom.pkl"):
        with open("caches/bloom.pkl", "rb") as f:
            __sites_and_fitnesses__ = pkl.load(f)
    else:
        __sites_and_fitnesses__ = {

        }
    __binding__ = BindingCalculator()

    def __init__(self) -> None:
        super().__init__()

    @classmethod
    def calculate_fitness(cls, sites, offset=331):
        sorted_sites = [site + offset for site in sites]
        sorted(sorted_sites)
        if tuple(sorted_sites) in cls.__sites_and_fitnesses__:
            return cls.__sites_and_fitnesses__[tuple(sorted_sites)]
        else:
            to_calculate_with = set(cls.__binding__.sites).intersection(set(sorted_sites))
            if to_calculate_with:
                antigenic_difference = 1 - cls.__binding__.binding_retained(sorted_sites)
            else:
                antigenic_difference = 0
            cls.__sites_and_fitnesses__[tuple(sorted_sites)] = antigenic_difference
            with open("caches/bloom.pkl", "wb") as f:
                pkl.dump(cls.__sites_and_fitnesses__, f)
            return antigenic_difference

    @classmethod
    def calculate_fitness_of_sequence(cls, sequence, offset=301):
        sites = sequence.get_mutated_site_indexes()
        return cls.calculate_fitness(sites)

    @classmethod
    def calculate_fitness_of_only_different_proteins(cls, sequence, offset=331):
        non_neutral_mutations = sequence.get_non_neutral()
        sites = [non_neutral_mutation["site"] for non_neutral_mutation in non_neutral_mutations]
        return cls.calculate_fitness(sites, offset)
