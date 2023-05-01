import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style

def get_avg_low_high(name, prob):
    escapes = []
    fitnesses = []
    for i in range(10):
        ga = pd.read_pickle(f"results/Jaime_lawrence/prob_{prob}/{name}_{i}.pkl")
        escapes.append([x[0] for x in ga.fitnesses])
        fitnesses.append([x[1] for x in ga.fitnesses])

    fit_low = [min(idx) for idx in zip(*fitnesses)]
    fit_high = [max(idx) for idx in zip(*fitnesses)]
    fit_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*fitnesses)]

    esc_low = [min(idx) for idx in zip(*escapes)]
    esc_high = [max(idx) for idx in zip(*escapes)]
    esc_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*escapes)]

    return fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg

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
    fig.suptitle(f"Interbreed with {name.capitalize()} at {int(float(prob)*100)}% for 10 runs")
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

    plt.savefig(f'results/figures/Figure_{name}_{prob}.png')

def get_mutations(name, prob):
    sites = [0]*200
    for i in range(10):
        ga = pd.read_pickle(f"results/Jaime_lawrence/prob_{prob}/{name}_{i}.pkl")
        for j in ga.top_strains[-1].get_non_neutral():
            sites[j["site"]] = sites[j["site"]] + 1
    return [x / 10 for x in sites]

def plot_mutation_per_site(name, prob, mut_avg, orig):
    style.use('fivethirtyeight')

    plt.rc('axes', titlesize=14)
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('legend', fontsize=12)

    fig = plt.figure(figsize=(12, 6))
    fig.suptitle(f"Interbreed with {name.capitalize()} at {int(float(prob)*100)}%")
    mut_plot = fig.add_subplot(1, 1, 1)
    mut_plot.set_title("Probability of Mutations per Site in RBD for 10 runs")

    xs = range(331, 531)

    mut_plot.set_xlabel("site")
    mut_plot.set_xticks(range(xs[0], xs[-1]+10, 50))

    mut_plot.set_ylabel("Probality of Mutation")
    mut_plot.bar(xs, mut_avg, color="teal", alpha=0.8, label="average ga top result")
    mut_plot.bar(xs, orig, color="blue", alpha=0.4, label=name.capitalize())

    fig.legend(loc="upper left")

    plt.savefig(f'results/figures/Mutations_per_Site_{name}_{prob}.png')

blank = [0]*200
# N501Y
alpha = blank[0:]
alpha[501-331] = 1
# K417N E484K N501Y
beta = blank[0:]
beta[417-331] = 1
beta[484-331] = 1
beta[501-331] = 1
# L452R T478K
delta = blank[0:]
delta[452-331] = 1
delta[478-331] = 1
# G339D S371L S373P S375F K417N N440K G446S S477N T478K E484A Q493R G496S Q498R N501Y Y505H
omicron = blank[0:]
omicron[339-331] = 1
omicron[371-331] = 1
omicron[373-331] = 1
omicron[375-331] = 1
omicron[417-331] = 1 
omicron[440-331] = 1 
omicron[446-331] = 1 
omicron[477-331] = 1 
omicron[478-331] = 1 
omicron[484-331] = 1 
omicron[493-331] = 1 
omicron[496-331] = 1 
omicron[498-331] = 1 
omicron[501-331] = 1 
omicron[505-331] = 1
strains = {'alpha':alpha, 'beta':beta, 'delta':delta, 'omicron':omicron}

for name in ['alpha', 'beta', 'delta', 'omicron']:
    for prob in ['0.01', '0.1', '0.5']:
        fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg = get_avg_low_high(name, prob)
        plot_fitness_escape_over_time(name, prob, fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg)
        plot_mutation_per_site(name, prob, get_mutations(name, prob), strains[name])


print("whatever")