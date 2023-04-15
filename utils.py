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

def translate_codon_sequence_to_aas(codon_sequence_to_translate, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    ret = []
    codon_string = str(codon_sequence_to_translate)
    for i in range(0, len(codon_string), 3):
        local_codon = codon_string[i:i+3]
        ret.append(translate_base_pairs_to_aas(local_codon, translation_list))
    return "".join(ret)