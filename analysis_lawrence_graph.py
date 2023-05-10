import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style

def get_avg_low_high():
    escapes = []
    fitnesses = []
    for i in range(100):
        ga = pd.read_pickle(f"results/Lawrence/interbreed_top/1/{i}.pkl")
        escapes.append([x[0] for x in ga.fitnesses])
        fitnesses.append([x[1] for x in ga.fitnesses])

    fit_low = [min(idx) for idx in zip(*fitnesses)]
    fit_high = [max(idx) for idx in zip(*fitnesses)]
    fit_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*fitnesses)]

    esc_low = [min(idx) for idx in zip(*escapes)]
    esc_high = [max(idx) for idx in zip(*escapes)]
    esc_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*escapes)]

    return fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg

def plot_fitness_escape_over_time(fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg):
    style.use('fivethirtyeight')

    plt.rc('axes', titlesize=14)
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('legend', fontsize=12)

    fig = plt.figure(figsize=(8, 6))
    fit_plot = fig.add_subplot(2, 1, 1)
    esc_plot = fig.add_subplot(2, 1, 2)
    fig.suptitle("Interbreed Top Two With 100% Probability")
    fit_plot.set_title("Spike Protein Fitness Over Time")
    esc_plot.set_title("Antigenic Escape Over Time")

    xs = range(len(fit_avg))

    fit_plot.set_xlabel("generation")
    fit_plot.set_xticks(range(0, int(max(xs))+10, 10))

    esc_plot.set_xlabel("generation")
    esc_plot.set_xticks(range(0, int(max(xs))+10, 10))

    fit_plot.set_ylabel("fitness score")
    esc_plot.set_ylabel("escape score")

    plt.tight_layout()

    fit_plot.fill_between(xs, fit_low, fit_high, alpha=0.2, color="green")
    fit_plot.plot(xs, fit_avg, linewidth=1, color="green")

    esc_plot.fill_between(xs, esc_low, esc_high, alpha=0.2, color="red")
    esc_plot.plot(xs, esc_avg, linewidth=1, color="red")

    #fit_plot.legend(loc="upper left")

    plt.savefig(f'results/figures/Figure_interbreed_top_100.png')


fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg = get_avg_low_high()
plot_fitness_escape_over_time(fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg)
