import random
import sys

def check_dot_bracket(ss):
    """
    Checks if the given secondary structure string is in the correct dot-bracket notation.
    It ensures that every opening bracket has a matching closing bracket, and vice versa.

    Parameters:
    ss (str): The secondary structure string in dot-bracket notation.

    Returns:
    list: A list of pairs, where each pair is a list of two indices representing a base pair.

    Raises:
    SystemExit: If the secondary structure string is not in the correct dot-bracket notation.
    """

    db_list = [['(', ')'], ['[', ']'], ['<', '>'], ['{', '}'], ['A', 'a'], ['B', 'b'], ['C', 'c'], ['D', 'd'], ['E', 'e']]
    allowed_characters = '()[]<>{}AaBbCcDdEe.&'

    for i in range(0, len(ss)):
        if ss[i] not in allowed_characters:
            print("Not allowed characters in structures. Check input file.")
            sys.exit()

    stack_list = []
    pairs_list = []

    # stack-pop for all versions of brackets form the db_list
    for i in range(0, len(db_list)):
        for c, s in enumerate(ss):
            if s == db_list[i][0]:
                stack_list.append(c)
            elif s == db_list[i][1]:
                if len(stack_list) == 0:
                    print("There is no opening bracket for nt position " + str(c + 1) + '-' + ss[c])
                    sys.exit()
                elif s == db_list[i][1]:
                    pairs_list.append([stack_list.pop(), c])
        if len(stack_list) > 0:
            err = stack_list.pop()
            print("There is no closing bracket for nt position " + str(err) + '-' + ss[err])
            sys.exit()

    return pairs_list


def check_seq_restr(restr):
    """
    Checks if the given sequence restraints string only contains allowed characters.

    Args:
        restr (str): The sequence restraints string.

    Raises:
        SystemExit: If the sequence restraints string contains not allowed characters.
    """

    allowed_characters = 'ACGUWSMKRYBDHVN-&'

    for i in range(0, len(restr)):
        if restr[i] not in allowed_characters:
            print("Not allowed characters in sequence restraints. Check input file.")
            sys.exit()


def check_length(ss, restr):
    """
    Checks if the secondary structure string and the sequence restraints string are of the same length.

    Args:
        ss (str): The secondary structure string.
        restr (str): The sequence restraints string.

    Raises:
        SystemExit: If the secondary structure string and the sequence restraints string are of different lengths.
    """

    if len(ss) != len(restr):
        print("Secondary structure and sequence restraints are of different length. Check input file.")
        sys.exit()


def get_nt_list(input_file):
    """
    Constructs a list of Nucleotide objects from the given InputFile object.

    Args:
        input_file (InputFile): The InputFile object containing the data from the input file.

    Returns:
        list_of_nt (list): A list of Nucleotide objects.
    """

    pair_list = input_file.pairs
    restr_seq = input_file.seq_restr

    list_of_nt = []

    for i in range(0, len(restr_seq)):
        obj = Nucleotide(number=i)
        obj.add_letter(nt_dictionary(restr_seq[i]))
        list_of_nt.append(obj)

    for i in range(0, len(pair_list)):
        open_br = pair_list[i][0]
        close_br = pair_list[i][1]
        list_of_nt[open_br].add_pair(close_br)
        list_of_nt[close_br].add_pair(open_br)
        nt_open = list_of_nt[open_br]
        nt_close = list_of_nt[close_br]

        for k in range(0, len(nt_open.letters)):
            pairing_let = can_pair(nt_open.letters[k])
            nt_open.add_pairing_l(pairing_let)

        for l in range(0, len(nt_close.letters)):
            pairing_let = can_pair(nt_close.letters[l])
            nt_close.add_pairing_l(pairing_let)

        nt_open.add_allowed_l(list(set(nt_close.pair_letters).intersection(nt_open.letters)))
        nt_close.add_allowed_l(list(set(nt_open.pair_letters).intersection(nt_close.letters)))

    for i in range(0, len(list_of_nt)):
        if list_of_nt[i].letters_allowed == None:
            list_of_nt[i].letters_allowed = list_of_nt[i].letters

    return list_of_nt


def nt_dictionary(nt):
    """
    Maps a nucleotide character to a list of nucleotide characters based on the IUPAC convention.

    Parameters:
    nt (str): The nucleotide character.

    Returns:
    list: A list of nucleotide characters.
    """

    constraints_dict = {
        'N': ['A', 'C', 'G', 'U'],
        'W': ['A', 'U'],
        'S': ['C', 'G'],
        'M': ['A', 'C'],
        'K': ['G', 'U'],
        'R': ['A', 'G'],
        'Y': ['C', 'U'],
        'B': ['C', 'G', 'U'],
        'D': ['A', 'G', 'U'],
        'H': ['A', 'C', 'U'],
        'V': ['A', 'C', 'G'],
        'C': ['C'],
        'A': ['A'],
        'G': ['G'],
        'U': ['U'],
        '&': ['&']
    }

    nt_all = constraints_dict[nt]

    return nt_all


def check_input_logic(nt_list):
    """
    Checks if the logic of the input is correct, specifically if the pairing restrictions make sense.

    Parameters:
    nt_list (list of Nucleotide): A list of Nucleotide objects.

    Raises:
    SystemExit: If the logic of the input is incorrect.
    """

    for i in range(0, len(nt_list)):
        if nt_list[i].pairs_with != None:
            if_pair(nt_list[i], nt_list[nt_list[i].pairs_with])


def if_pair(nt1, nt2):
    """
    Checks if two nucleotides can form a base pair according to the pairing rules.

    Parameters:
    nt1 (Nucleotide): The first nucleotide object.
    nt2 (Nucleotide): The second nucleotide object.

    Raises:
    SystemExit: If the nucleotides cannot form a base pair.
    """

    pair_dict = {'A': ['U'], 'U': ['G', 'A'], 'G': ['U', 'C'], 'C': ['G']}

    store = False
    for i in range(0, len(nt1.letters)):
        nt1l = nt1.letters[i]
        for j in range(0, len(nt2.letters)):
            nt2l = nt2.letters[j]
            if nt2l in pair_dict[nt1l]:
                store = True

    if store == True:
        pass
    else:
        print("Wrong restraints in the input file. Nucleotide " + str(nt1.number + 1) + " " + str(nt1.letters) + ", cannot pair with nucleotide " + str(nt2.number + 1) + " " + str(nt2.letters))
        sys.exit()


def wc_pair(nt1):
    """
    Returns the Watson-Crick pair of a given nucleotide character.

    Args:
    nt1 (str): The nucleotide character.

    Returns:
    str: The Watson-Crick pair of the nucleotide character.
    """

    pair_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}

    nt2 = pair_dict[nt1]

    return nt2


def can_pair(nt):
    """
    Returns a list of nucleotide characters that can form a base pair with a given nucleotide character.

    Args:
    nt (str): The nucleotide character.

    Returns:
    list: A list of nucleotide characters that can pair with the given nucleotide.
    """

    pair_dict = {'A': ['U'], 'U': ['G', 'A'], 'G': ['U', 'C'], 'C': ['G']}

    pairing_l = pair_dict[nt]

    return pairing_l



def random_sequence_generator(nt_list, input_file, sim_options):
    """
    Generates a random sequence that satisfies the sequence restraints.

    Args:
    nt_list (list): A list of Nucleotide objects representing sequence restraints.
    input_file (InputFile): The InputFile object containing the data from the input file.

    Returns:
    str: The randomly generated sequence that complies with the specified restraints.
    """

    seq_l = list(input_file.seq_restr).copy()
    pair_list = sorted(input_file.pairs.copy())

    for i in range(0, len(nt_list)):
        if nt_list[i].pairs_with == None:
            random.choice(nt_dictionary(seq_l[i]))

    for i in range(0, len(pair_list)):
        nt1 = nt_list[pair_list[i][0]].letters_allowed
        nt2 = nt_list[pair_list[i][1]].letters_allowed
        nt1.sort()
        nt2.sort()
        nt1num = nt_list[pair_list[i][0]].number
        nt2num = nt_list[pair_list[i][1]].number

        allowed_choices = allowed_choice(nt1, sim_options.nt_percentages)

        seq_l[nt1num] = random.choices(nt1, weights=allowed_choices)[0]
        seq_l[nt2num] = wc_pair(seq_l[nt1num])

    for i in range(0, len(seq_l)):
        if seq_l[i] not in ["A", "C", "G", "U"]:
            seq_l[i] = random.choice(nt_dictionary(seq_l[i]))

    result_sequence = ''.join(seq_l)

    return result_sequence


def initial_sequence_generator(nt_list, input_file, sim_options):
    """
    Generates the initial RNA sequence. The function iterates over the secondary structure
    of the RNA (input.sec_struct) and generates a corresponding sequence based on the type
    of structure at each position (e.g., pair, unpair). The generated sequence is returned
    as a string.

    Args:
        input_file (object): An InputFile object that contains the secondary structure of the RNA.

    Returns:
        str: The generated initial RNA sequence.
    """

    seq_l = list(input_file.seq_restr).copy()
    pair_list = sorted(input_file.pairs.copy())

    # assign A to all nonbonded nucleotides
    for i in range(0, len(nt_list)):
        if nt_list[i].pairs_with == None and "A" in nt_list[i].letters:
            seq_l[i] = "A"

    # Loop boosting - G (or if not possible = U) as the first of nonbonded nts, unless its (.( - single bulge, then stay A
    for i in range(1, len(nt_list) - 1):
        if (nt_list[i].pairs_with == None and nt_list[i - 1].pairs_with != None) and (nt_list[i].pairs_with == None and nt_list[i + 1].pairs_with == None):
            if "G" in nt_list[i].letters:
                seq_l[i] = "G"
            elif "U" in nt_list[i].letters:
                seq_l[i] = "U"

    if sim_options.acgu_percentages == 'off':

        # paired nucleotides handling, if possible GC, if not AU
        for i in range(0, len(pair_list)):
            nt1 = nt_list[pair_list[i][0]].letters_allowed
            nt2 = nt_list[pair_list[i][1]].letters_allowed
            nt1num = nt_list[pair_list[i][0]].number
            nt2num = nt_list[pair_list[i][1]].number
            if ("C" in nt1 and "G" in nt1) and ("C" in nt2 and "G" in nt2):
                seq_l[nt1num] = random.choice(["C", "G"])
                seq_l[nt2num] = wc_pair(seq_l[nt1num])
            elif "C" in nt1 and "G" in nt2:
                seq_l[nt1num] = "C"
                seq_l[nt2num] = "G"
            elif "G" in nt1 and "C" in nt2:
                seq_l[nt1num] = "G"
                seq_l[nt2num] = "C"
            elif ("A" in nt1 and "U" in nt1) and ("A" in nt2 and "U" in nt2):
                seq_l[nt1num] = random.choice(["A", "U"])
                seq_l[nt2num] = wc_pair(seq_l[nt1num])
            elif "A" in nt1 and "U" in nt2:
                seq_l[nt1num] = "A"
                seq_l[nt2num] = "U"
            elif "U" in nt1 and "A" in nt2:
                seq_l[nt1num] = "U"
                seq_l[nt2num] = "A"

    elif sim_options.acgu_percentages == 'on':
        for i in range(0, len(pair_list)):
            nt1 = nt_list[pair_list[i][0]].letters_allowed
            nt2 = nt_list[pair_list[i][1]].letters_allowed
            nt1.sort()
            nt2.sort()
            nt1num = nt_list[pair_list[i][0]].number
            nt2num = nt_list[pair_list[i][1]].number

            allowed_choices = allowed_choice(nt1, sim_options.nt_percentages)

            seq_l[nt1num] = random.choices(nt1, weights=allowed_choices)[0]

            seq_l[nt2num] = wc_pair(seq_l[nt1num])

    # random choice (from allowed letters) of whatever left
    for i in range(0, len(seq_l)):
        if seq_l[i] not in ["A", "C", "G", "U"]:
            seq_l[i] = random.choice(nt_dictionary(seq_l[i]))

    result_sequence = ''.join(seq_l)

    if input_file.seed_seq:
        result_sequence = input_file.seed_seq

    return result_sequence


def allowed_choice(allowed, percs):
    """
    Determines the allowed mutations based on nucleotide percentages.

    Args:
    allowed (list): List of allowed mutations.
    percs (dict): Dictionary with the percentage of each nucleotide.

    Returns:
    list: List of weights for the allowed mutations, corresponding to the nucleotide percentages.
    """

    return [percs[nt] for nt in allowed]











class Nucleotide:
    """
    Represents a nucleotide in an RNA sequence.

    Attributes:
    number (int): The position of the nucleotide in the sequence.
    letters (list): The possible letters that this nucleotide could be (A, C, G, or U).
    pairs_with (int, optional): The position of the nucleotide that this one pairs with, if any.
    pair_letters (list, optional): The possible letters that the paired nucleotide could be.
    letters_allowed (list): The letters that are allowed for this nucleotide, considering constraints.

    The class provides functionalities to manage the nucleotide's properties and constraints.
    """

    def __init__(self, number):
        """
        Initialize a Nucleotide instance.

        Parameters:
        number (int): The position of the nucleotide in the sequence.
        """
        self.number = number
        self.letters = []
        self.pairs_with = None
        self.pair_letters = []
        self.letters_allowed = None

    def add_letter(self, letter):
        """
        Adds a nucleotide letter to the current list of possible letters for this nucleotide.

        Parameters:
        letter (str): The nucleotide letter to be added.
        """
        self.letters = letter

    def add_pair(self, pair):
        """
        Sets the pairing nucleotide's position for this nucleotide.

        Parameters:
        pair (int): The position of the nucleotide that this one pairs with.
        """
        self.pairs_with = pair

    def add_pairing_l(self, pairing_l):
        """
        Adds to the list of possible letters for the nucleotide's pairing partner.

        Parameters:
        pairing_l (list): A list of possible letters for the paired nucleotide.
        """
        self.pair_letters += pairing_l
        self.pair_letters = list(set(self.pair_letters))

    def add_allowed_l(self, list):
        """
        Sets the list of allowed letters for this nucleotide, considering any constraints.

        Parameters:
        list (list): The list of letters that are allowed for this nucleotide.
        """
        self.letters_allowed = list



