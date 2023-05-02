import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style

def get_avg_low_high_top_fitness(muts):
    esc_avg = []
    esc_low = []
    esc_high = []
    fit_avg = []
    fit_low = []
    fit_high = []

    for m in muts:
        e_avg = 0
        e_low = 1000
        e_high = 0
        f_avg = 0
        f_low = 1000
        f_high = 0
        for i in range(100):
            ga = pd.read_pickle(f"results/Lawrence/mutations/{m}/{i}.pkl")
            esc = ga.fitnesses[-1][0]
            fit = ga.fitnesses[-1][1]
            e_avg = e_avg + esc
            f_avg = f_avg + fit
            if esc < e_low:
                e_low = esc
            if esc > e_high:
                e_high = esc
            if fit < f_low:
                f_low = fit
            if fit > f_high:
                f_high = fit
        esc_avg.append(e_avg/100)
        esc_low.append(e_low)
        esc_high.append(e_high)
        fit_avg.append(f_avg/100)
        fit_low.append(f_low)
        fit_high.append(f_high)

    return fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg

def plot_fitness_escape_over_mutations(xs, fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg):
    style.use('fivethirtyeight')

    plt.rc('axes', titlesize=14)
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('legend', fontsize=12)

    fig = plt.figure(figsize=(8, 6))
    fit_plot = fig.add_subplot(2, 1, 1)
    esc_plot = fig.add_subplot(2, 1, 2)
    fig.suptitle(f"Average Fitness per Mutations per Generation for 100 runs")
    fit_plot.set_title("Avg. Spike Protein Fitness after 100 Generations")
    esc_plot.set_title("Avg. Antigenic Escape after 100 Generations")

    fit_plot.set_xlabel("mutations per generation")
    fit_plot.set_xticks(range(len(xs)))
    fit_plot.set_xticklabels(xs)

    esc_plot.set_xlabel("mutations per generation")
    esc_plot.set_xticks(range(len(xs)))
    esc_plot.set_xticklabels(xs)

    fit_plot.set_ylabel("fitness score")
    esc_plot.set_ylabel("escape score")

    plt.tight_layout()

    fit_plot.bar(range(len(xs)), fit_avg, color="teal", alpha=0.9)
    fit_plot.errorbar(range(len(xs)), fit_avg, yerr=[abs(fit_high[i]-fit_low[i]) for i in range(len(fit_high))], capsize=5, elinewidth=0.5, ecolor="black", ls="none")
    esc_plot.bar(range(len(xs)), esc_avg, color="red", alpha=0.8)
    esc_plot.errorbar(range(len(xs)), esc_avg, yerr=[abs(esc_high[i]-esc_low[i]) for i in range(len(esc_high))], capsize=5, elinewidth=0.5, ecolor="black", ls="none")

    plt.savefig(f'results/figures/Figure_Mutations.png')


# main
muts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg = get_avg_low_high_top_fitness(muts)
plot_fitness_escape_over_mutations(muts, fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg)

print("whatever")