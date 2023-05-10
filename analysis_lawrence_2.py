import os
import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
from sequence import *
from ga import *


def process_directory(path):
    dirs = os.listdir(path)
    overall_random_fit = 0
    top_fit_avg = 0
    top_antigen_avg = 0
    top_fit_path = ""
    for dir in dirs:
        inner_path = f"{path}/{dir}"
        antigenic_vals = []
        fitness_vals = []
        for pickle in os.listdir(inner_path):
            pickle_path = f"{inner_path}/{pickle}"
            with open(pickle_path, "rb") as f:
                gene_a = pkl.load(f)
            esc = gene_a.fitnesses[-1][0]
            fit = gene_a.fitnesses[-1][1]
            antigenic_vals.append(esc)
            fitness_vals.append(fit)
        print("krewl")
        fitness_avg = sum(fitness_vals) / len(fitness_vals)
        antgen_avg = sum(antigenic_vals) / len(antigenic_vals)
        overall = fitness_avg + antgen_avg
        if overall > overall_random_fit:
            overall_random_fit = overall
            top_fit_avg = fitness_avg
            top_antigen_avg = antgen_avg
            top_fit_path = inner_path
    return top_fit_avg, top_antigen_avg, top_fit_path

if __name__ == '__main__':
    path = "results/Lawrence/mutations"
    fit_avg, antigen_avg, fit_path = process_directory(path)
    print("krewl")