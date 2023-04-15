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