import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
from utils import load_sequence, translate_codon_sequence_to_aas
from matplotlib import style


def get_avg_low_high(name):
    escapes = []
    fitnesses = []
    for i in range(10):
        ga = pd.read_pickle(f"results/Jeb/{name}/{name}_{i}.pkl")
        escapes.append([x[0] for x in ga.fitnesses])
        fitnesses.append([x[1] for x in ga.fitnesses])

    fit_low = [min(idx) for idx in zip(*fitnesses)]
    fit_high = [max(idx) for idx in zip(*fitnesses)]
    fit_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*fitnesses)]

    esc_low = [min(idx) for idx in zip(*escapes)]
    esc_high = [max(idx) for idx in zip(*escapes)]
    esc_avg = [sum(sub_list) / len(sub_list) for sub_list in zip(*escapes)]

    return fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg

def plot_fitness_escape_over_time(name, fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg):
    style.use('fivethirtyeight')

    plt.rc('axes', titlesize=14)
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('legend', fontsize=12)

    fig = plt.figure(figsize=(8, 6))
    fit_plot = fig.add_subplot(2, 1, 1)
    esc_plot = fig.add_subplot(2, 1, 2)
    fig.suptitle(f"Starting Population {name.capitalize()} (10 runs)")
    fit_plot.set_title("Avg. Spike Protein Fitness Over Time")
    esc_plot.set_title("Avg. Antigenic Escape Over Time")

    xs = range(len(fit_avg))

    fit_plot.set_xlabel("generation")
    fit_plot.set_xticks(range(0, int(max(xs))+10, 100))

    esc_plot.set_xlabel("generation")
    esc_plot.set_xticks(range(0, int(max(xs))+10, 100))

    fit_plot.set_ylabel("fitness score")
    esc_plot.set_ylabel("escape score")

    plt.tight_layout()

    fit_plot.fill_between(xs, fit_low, fit_high, alpha=0.2, color="green")
    fit_plot.plot(xs, fit_avg, linewidth=1, color="green")

    esc_plot.fill_between(xs, esc_low, esc_high, alpha=0.2, color="red")
    esc_plot.plot(xs, esc_avg, linewidth=1, color="red")

    #fit_plot.legend(loc="upper left")

    plt.savefig(f'results/figures/Figure_{name}.png')

def get_mutations(name):
    sites = [0]*200
    for i in range(10):
        ga = pd.read_pickle(f"results/Jeb/{name}/{name}_{i}.pkl")
        for j in ga.top_strains[-1].get_non_neutral():
            sites[j["site"]] = sites[j["site"]] + 1
    return [x / 10 for x in sites]

def plot_mutation_per_site(name, mut_avg):
    style.use('fivethirtyeight')

    plt.rc('axes', titlesize=14)
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('legend', fontsize=12)

    fig = plt.figure(figsize=(12, 6))
    fig.suptitle(f"Starting Population {name.capitalize()}")
    mut_plot = fig.add_subplot(1, 1, 1)
    mut_plot.set_title("Probability of Non-Neutral Mutations per Site in RBD for 10 runs")

    xs = range(331, 531)

    mut_plot.set_xlabel("site")
    mut_plot.set_xticks(range(xs[0], xs[-1]+10, 50))

    mut_plot.set_ylabel("Probability of Non-Neutral Mutation")
    mut_plot.bar(xs, mut_avg, color="teal", alpha=0.8, label="average ga top result")

    fig.legend(loc="upper left")

    plt.savefig(f'results/figures/Mutations_per_Site_start_{name}.png')

def get_best(names):
    best_score = 0
    for name in names:
        for i in range(10):
            ga = pd.read_pickle(f"results/Jeb/{name}/{name}_{i}.pkl")
            score = ga.top_strains[-1].generate_fitness_score_non_neutral(1, 0.5)
            if score > best_score:
                best_score = score
                strain = ga.top_strains[-1]
    return strain


names = ['antigen_neutral', 'neutral', 'neutral_1_hops_antigenic', 'neutral_1_mutation']

#for name in names:
    #fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg = get_avg_low_high(name)
    #plot_fitness_escape_over_time(name, fit_low, fit_high, fit_avg, esc_low, esc_high, esc_avg)
    #plot_mutation_per_site(name, get_mutations(name))

rbd = translate_codon_sequence_to_aas(get_best(names).get_sequence())
print(rbd)

#MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPTXLRVVRDNAVFNCXTRPAVXQAKRTPDVFPHEDILTIAKANAFAADSXFFKRALDLNEGNGRTKACDCFVLLLPQPSGDQIIKNITEAAQNVLYPHGXXLPCVAINSXKRDATAVKTDDMRFQSGDQNKTQTFYLDVGNQLLEXSGKRAKRLQRVXXHVFXLPISFSRSYQXIKIPHPKSDLIFQVVLSRPAVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT

print("whatever")