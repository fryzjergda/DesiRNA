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



