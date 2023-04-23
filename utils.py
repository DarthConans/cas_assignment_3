AMINO_CHARS_NUMS = ["1", "2", "3", "4"]

AMINO_CHARS_LETTERS = {"G", "U", "C", "A"}

translate_dict = {
    "G": "1",
    "U": "2",
    "C": "3",
    "A": "4",
}


def char_to_number(char):
    return translate_dict[char]


translate_back_dict = {
    "1": "G",
    "2": "U",
    "3": "C",
    "4": "A",
}


def number_to_char(char):
    return translate_back_dict[char]


def translate_string_to_numbers(to_translate):
    ret = ""
    for char in to_translate:
        ret += char_to_number(char)
    return ret


def translate_sequence_to_numbers(to_translate):
    return set([translate_string_to_numbers(trans) for trans in to_translate])


def load_sequence(TRANSLATE_FLAG=False, file_name="data/spike_protein_base_pairs.txt"):
    with open(file_name, "r") as f:
        untranslated = f.read().replace(" ", "").replace("\n", "").strip().upper()[330*3:523*3].replace("T", "U")
    if TRANSLATE_FLAG:
        return tranlate_aas_to_base_pairs(untranslated)
    else:
        return untranslated

def generate_known_strains():
    rbd = load_sequence()
    # N501Y
    alpha = replace_site_codon(rbd, ["N501Y"])
    # K417N E484K N501Y
    beta = replace_site_codon(rbd, ["K417N", "E484K", "N501Y"])
    # L452R T478K
    delta = replace_site_codon(rbd, ["L452R", "T478K"])
    # G339D S371L S373P S375F K417N N440K G446S S477N T478K E484A Q493R G496S Q498R N501Y Y505H
    # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9376347/
    omicron = replace_site_codon(rbd, ["G339D", "S371L", "S373P", "S375F", "K417N",
                                       "N440K", "G446S", "S477N", "T478K", "E484A",
                                       "Q493R", "G496S", "Q498R", "N501Y", "Y505H"])

    return [alpha, beta, delta, omicron]

def replace_site_codon(sequence, mutations):
    offset = 331
    new_sequence = sequence[0:]
    for i in range(0, len(mutations)):
        site = int(mutations[i][1:-1])
        mutation = mutations[i][-1]
        new_sequence = (new_sequence[0:3*(site-offset)]
                        +tranlate_aas_to_base_pairs([mutation])
                        +new_sequence[3*(site-offset+1):])
    return new_sequence

EQUIVALENT_SEQUENCES = [
    ({"GUU", "GUC", "GUA", "GUG"}, "V"),
    ({"GCU", "GCC", "GCA", "GCG"}, "A"),
    ({"GAU", "GAC"}, "D"),
    ({"GAA", "GAG"}, "E"),
    ({"GGU", "GGC", "GGA", "GGG"}, "G"),
    ({"UUU", "UUC"}, "F"),
    ({"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"}, "L"),
    ({"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"}, "S"),
    ({"UAU", "UAC"}, "Y"),
    ({"UAA", "UAG", "UGA"}, "X"),
    ({"UGU", "UGC"}, "C"),
    ({"UGG"}, "W"),
    ({"CCU", "CCC", "CCA", "CCG"}, "P"),
    ({"CAU", "CAC"}, "H"),
    ({"CAA", "CAG"}, "Q"),
    ({"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"}, "R"),
    ({"AUU", "AUC", "AUA"}, "I"),
    ({"AUG"}, "M"),
    ({"ACU", "ACC", "ACA", "ACG"}, "T"),
    ({"AAU", "AAC"}, "N"),
    ({"AAA", "AAG"}, "K"),
]


def tranlate_aas_to_base_pairs(aas, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    ret = []
    for aa in aas:
        found = 0
        for base_pair, prot in translation_list:
            if prot == aa:
                ret.append(list(base_pair)[0])
                found = 1
        if not found:
            raise ValueError(f"AMINO ACID {aa} DIDN'T MATCH")
    return "".join(ret)


def translate_base_pairs_to_aas(base_pairs_to_translate, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    for base_pairs, aa in translation_list:
        if base_pairs_to_translate in base_pairs:
            return aa
    raise ValueError(f"INVALID BASE PAIR SEQUENCE {base_pairs_to_translate} DIDN'T MATCH")

def translate_numbers_to_string(to_translate):
    ret = ""
    for char in str(to_translate):
        ret += number_to_char(char)
    return ret


def translate_codon_sequence_to_aas(codon_sequence_to_translate, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    ret = []
    codon_string = str(codon_sequence_to_translate)
    for i in range(0, len(codon_string), 3):
        local_codon = codon_string[i:i + 3]
        ret.append(translate_base_pairs_to_aas(local_codon, translation_list))
    return "".join(ret)

