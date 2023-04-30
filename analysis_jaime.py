import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style

def get_avg_low_high(name, prob):
    escapes = []
    fitnesses = []
    for i in range(10):
        ga = pd.read_pickle(f"results/Jaime_new/prob_{prob}/{name}_{i}.pkl")
        escapes.append([x[0] for x in ga.fitnesses])
        fitnesses.append([x[1] for x in ga.fitnesses])

    fit_low = [min(idx) for idx in zip(*fitnesses)]
    fit_high = [max(idx) for idx in zip(*fitnesses)]
    fit_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*fitnesses)]

    esc_low = [min(idx) for idx in zip(*escapes)]
    esc_high = [max(idx) for idx in zip(*escapes)]
    esc_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*escapes)]

    return fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg

def get_mutations(name, prob):
    mutations = []
    for i in range(10):
        ga = pd.read_pickle(f"results/Jaime_new/prob_{prob}/{name}_{i}.pkl")
        mutations.append([])


def plot_fitness_escape_over_time(name, prob, fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg):
    style.use('fivethirtyeight')

    plt.rc('axes', titlesize=14)
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('legend', fontsize=12)

    fig = plt.figure(figsize=(8, 6))
    fit_plot = fig.add_subplot(2, 1, 1)
    esc_plot = fig.add_subplot(2, 1, 2)
    fig.suptitle(f"Interbreed with {name} at {int(float(prob)*100)}%")
    fit_plot.set_title("Spike Protein Fitness Over Time")
    esc_plot.set_title("Antigenic Escape Over Time")

    xs = range(len(fit_avg))

    fit_plot.set_xlabel("generation")
    fit_plot.set_xticks(range(0, int(max(xs))+10, 10))

    esc_plot.set_xlabel("generation")
    esc_plot.set_xticks(range(0, int(max(xs))+10, 10))

    fit_plot.set_ylabel("fitness score")
    #fit_plot.set_yscale("log")

    esc_plot.set_ylabel("escape score")
    #esc_plot.set_yscale("log")

    plt.tight_layout()

    fit_plot.fill_between(xs, fit_low, fit_high, alpha=0.2, color="green")
    fit_plot.plot(xs, fit_avg, linewidth=1, color="green")

    esc_plot.fill_between(xs, esc_low, esc_high, alpha=0.2, color="red")
    esc_plot.plot(xs, esc_avg, linewidth=1, color="red")

    #fit_plot.legend(loc="upper left")

    plt.savefig(f'results/figures/Figure_{name}_{prob}.png')

for name in ['alpha', 'beta', 'delta', 'omicron']:
    for prob in ['0.01', '0.1', '0.5']:
        fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg = get_avg_low_high(name, prob)
        plot_fitness_escape_over_time(name, prob, fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg)
        #plot_mutation_per_site(name, prob)


print("whatever")