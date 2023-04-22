import numpy as np
import pandas as pd
from os.path import exists
import pickle as pkl


def load_fitness_frame():
    frame = pd.read_csv("data/aamut_fitness_all.csv")
    frame = frame[frame["gene"] == "S"]
    frame = frame[["gene", "clade_founder_aa", "mutant_aa", "aa_site", "delta_fitness"]]
    return frame


class bloom_fitness:
    __fitness_frame__ = load_fitness_frame()
    __cache_path__ = "caches/bloom_fitness.pkl"
    if exists(__cache_path__):
        with open(__cache_path__, "rb") as f:
            __cache__ = pkl.load(f)
    else:
        __cache__ = {

        }

    @classmethod
    def get_fitness_change(cls, mutation_data, offset=331):
        total_change = 0
        for mutation_datum in mutation_data[0]:
            site = mutation_datum[0] + offset
            original_amino_acid = mutation_datum[1]
            new_amino_acid = mutation_datum[2]
            fitness_row = cls.__fitness_frame__[(cls.__fitness_frame__["aa_site"] == site) & (
                    cls.__fitness_frame__["clade_founder_aa"] == original_amino_acid) & (
                                                        cls.__fitness_frame__["mutant_aa"] == new_amino_acid)]
            if not fitness_row.empty:
                total_change += fitness_row["delta_fitness"].iloc[0]
            else:
                return np.NAN
        return total_change

    @classmethod
    def __fitness__(cls, site, old_amino, new_amino):
        fitness = 0
        dict_tuple = (site, old_amino, new_amino)
        if dict_tuple in cls.__cache__.keys():
            fitness = cls.__cache__[dict_tuple]
        else:
            fitness_row = cls.__fitness_frame__[(cls.__fitness_frame__["aa_site"] == site) &
                                                (cls.__fitness_frame__["clade_founder_aa"] == old_amino) &
                                                (cls.__fitness_frame__["mutant_aa"] == new_amino)
                                                ]
            if not fitness_row.empty:
                fitness = fitness_row["delta_fitness"].iloc[0]
            cls.__cache__[dict_tuple] = fitness
            with open(cls.__cache_path__, "wb") as f:
                pkl.dump(cls.__cache__, f)
        return fitness

    @classmethod
    def get_sequence_fitness_change(cls, sequence, offset=331):
        non_neutrals = sequence.get_non_neutral()
        total_change = 0
        """
        
                    "site": i,
                    "old_amino": old_amino,
                    "new_amino": new_amino
                    """
        for non_neutral in non_neutrals:
            site = non_neutral["site"] + offset
            old_amino = non_neutral["old_amino"]
            new_amino = non_neutral["new_amino"]
            fitness = cls.__fitness__(site, old_amino, new_amino)
            total_change += fitness
        for neutral in sequence.get_neutral():
            site = neutral["site"]
            old_amino = new_amino = neutral["amino"]
            fitness = cls.__fitness__(site, old_amino, new_amino)
            total_change += fitness
        return total_change

if __name__ == '__main__':
    frame = load_fitness_frame()
    print("krewl")