"""
This module provides utility functions for working with RNA secondary structures and sequences.

Functions:
- check_dot_bracket(ss): Checks if a secondary structure string is in the correct dot-bracket notation.
- check_seq_restr(restr): Checks if a sequence restraints string contains allowed characters.
- check_length(ss, restr): Checks if the secondary structure and sequence restraints strings have the same length.
- get_nt_list(input_file): Constructs a list of Nucleotide objects from an InputFile object.
- nt_dictionary(nt): Maps a nucleotide character to a list of nucleotide characters based on the IUPAC convention.
- check_input_logic(nt_list): Checks if the logic of the input is correct.
- if_pair(nt1, nt2): Checks if two nucleotides can form a base pair.
- wc_pair(nt1): Returns the Watson-Crick pair of a nucleotide character.
- can_pair(nt): Returns a list of nucleotide characters that can form a base pair with a given nucleotide character.
- random_sequence_generator(nt_list, input_file, sim_options): Generates a random sequence that satisfies sequence restraints.
- initial_sequence_generator(nt_list, input_file, sim_options): Generates the initial RNA sequence based on secondary structure.
- allowed_choice(allowed, percs): Determines allowed mutations based on nucleotide percentages.
- get_rep_temps(sim_options): Generates temperature shelves for each replica in the simulation.
- generate_initial_list(nt_list, input_file, sim_options): Generates an initial list of RNA sequences.
- generate_initial_list_random(nt_list, input_file, sim_options): Generates an initial list of random RNA sequences.
- get_mutation_position(seq_obj, available_positions, sim_options, input_file): Determines a position to mutate.
- expand_cases(cases, max_value, range_expansion=3): Expands a set of cases within a specified range.
- mutate_sequence(sequence_obj, nt_list, sim_options, input_file): Mutates a sequence and returns the mutated sequence.
- get_alt_pairs(input_file): Processes alternative secondary structures and updates the input_file object with a flattened list of base pairs.
- get_pairs_for_graphs(input_file): Merges base pairs from the primary and alternative secondary structures for graph representation.
- find_distinct_graphs(pairs): Identifies distinct graphs (connected components) from a list of base pairs.
- is_bipartite_and_color(graph): Determines if a graph is bipartite and assigns colors to its vertices accordingly.
- generate_acgu_states(coloring): Generates possible ACGU mappings for a graph based on its coloring.
- generate_graphs(pairs): Generates graphs from RNA base pairs and determines their bipartiteness and ACGU states.
- filter_states_by_constraints(graph_numbers, possible_states, sequence_restraint): Filters potential nucleotide states by sequence constraints.
- update_graphs(input_file): Updates the graphs in an InputFile object based on sequence constraints.
- get_allsnakes(input_file): Extracts all 'snake' sequences from the graphs of an InputFile object.
- get_alt_mcc(simulation_data, input_file): Calculates MCC scores for alternative secondary structures against a sequence.
- find_diff_position(str1, str2): Finds the first position where two strings differ.

Classes:
- Nucleotide: Represents a nucleotide with properties like letters, pairs_with, etc.
- InputFile: Represents an input file with secondary structure, sequence restraints, etc.
- ScoreSeq: Represents an RNA sequence with scoring information.
"""

import random
import sys

from collections import deque

import RNA

import numpy as np

from utils import energy_scores as es
from utils.sim_score import SimScore


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


def get_alt_pairs(input_file):
    """
    Processes alternative secondary structures and updates the input_file object with a flattened list of base pairs.

    This function iterates over all alternative secondary structures defined in the input_file object,
    checks their dot-bracket notation for correctness, and compiles a comprehensive list of all base pairs
    encountered across these structures. It flattens this list and assigns it back to the input_file object.

    Args:
        input_file (InputFile): The input file object containing alternative secondary structures.

    Returns:
        None: The function directly modifies the input_file object, adding a flattened list of alternative base pairs.
    """

    alt_pairs = []
    for i in range(0, len(input_file.alt_sec_structs)):
        alt_pairs.append(check_dot_bracket(input_file.alt_sec_structs[i]))

    input_file.alt_pairs = alt_pairs
    flattened_list = [pair for sublist in alt_pairs for pair in sublist]
    input_file.alt_pairs = flattened_list


def get_pairs_for_graphs(input_file):
    """
    Merges base pairs from the primary and alternative secondary structures for graph representation.

    This function calls get_alt_pairs to process alternative structures and then combines these pairs with the
    primary structure's pairs. It ensures that all pairs are unique and filters out any redundant pairs,
    particularly focusing on pairs that are exclusive to the alternative structures. These exclusive pairs
    are used to update the input_file object for further processing.

    Args:
        input_file (InputFile): The input file object with primary and alternative secondary structures.

    Returns:
        list: A list of filtered and corrected base pairs for use in graph representations.
    """

    list1 = input_file.pairs
    get_alt_pairs(input_file)

    list2 = input_file.alt_pairs
    pairs = []
    for pair1 in list1:
        for pair2 in list2:
            if (pair1[0] in pair2 or pair1[1] in pair2) and pair1 != pair2:
                if pair1 not in pairs:
                    pairs.append(pair1)
                if pair2 not in pairs and pair2 != pair1:
                    pairs.append(pair2)

    merged_list = list(set([tuple(pair) for pair in list1] + [tuple(pair) for pair in list2]))

    flattened_elements = [element for pair in merged_list for element in pair]

    excluded_pairs = []
    filtered_pairs_corrected = []

    for pair in merged_list:
        if flattened_elements.count(pair[0]) > 1 or flattened_elements.count(pair[1]) > 1:
            filtered_pairs_corrected.append(list(pair))
        else:
            excluded_pairs.append(list(pair))

    input_file.excluded_alt_pairs = excluded_pairs

    return filtered_pairs_corrected


def find_distinct_graphs(pairs):
    """
    Identifies distinct graphs (connected components) from a list of base pairs.

    Utilizing a depth-first search algorithm, this function explores the connectivity among bases defined by
    the provided pairs, thereby identifying distinct connected components within the RNA structure. Each distinct
    graph represents a set of bases that are interconnected through base pairing, enabling further analysis of
    structural components.

    Args:
        pairs (list of tuple): A list of base pairs, where each pair is a tuple of two integers representing connected bases.

    Returns:
        list of set: A list of sets, with each set containing the bases belonging to a distinct connected component (graph).
    """

    graph = {}
    for pair in pairs:
        a, b = pair
        if a not in graph:
            graph[a] = []
        if b not in graph:
            graph[b] = []
        graph[a].append(b)
        graph[b].append(a)

    visited = set()
    distinct_graphs = []

    for start_node in graph.keys():
        if start_node not in visited:
            current_graph = set()
            stack = [start_node]

            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    current_graph.add(node)
                    for neighbour in graph[node]:
                        if neighbour not in visited:
                            stack.append(neighbour)

            distinct_graphs.append(current_graph)

    return distinct_graphs


def is_bipartite_and_color(graph):
    """
    Determines if a graph is bipartite and assigns colors to vertices accordingly.

    Args:
        graph (dict): A dictionary representing the graph where keys are vertex identifiers and values are lists of adjacent vertices.

    Returns:
        tuple: A pair where the first element is a boolean indicating whether the graph is bipartite, and the second element is a dictionary of vertex color assignments if the graph is bipartite, or an empty dictionary otherwise.
    """

    color = {}
    for vertex in graph:
        if vertex not in color:
            queue = deque([vertex])
            color[vertex] = 0
            while queue:
                v = queue.popleft()
                for u in graph[v]:
                    if u not in color:
                        color[u] = 1 - color[v]
                        queue.append(u)
                    elif color[u] == color[v]:
                        return False, {}

    return True, color


def generate_acgu_states(coloring):
    """
    Generates possible ACGU states for each bipartite graph based on its coloring.

    Args:
        coloring (dict): A dictionary where keys are graph identifiers and values are dictionaries representing the coloring of vertices in the graph.

    Returns:
        dict: A dictionary with keys as graph identifiers and values as information about possible ACGU mappings, including the sequence of nucleotides, their colors, and the corresponding ACGU states.
    """

    possible_mappings = [
        {0: 'A', 1: 'U'},
        {0: 'U', 1: 'A'},
        {0: 'C', 1: 'G'},
        {0: 'G', 1: 'C'}
    ]
    states_per_graph = {}
    for state, info in coloring.items():
        if info != "Not bipartite":
            sorted_nodes = sorted(info.keys())
            sorted_colors = [info[node] for node in sorted_nodes]
            states = []
            for mapping in possible_mappings:
                states.append([mapping[color] for color in sorted_colors])
            states_per_graph[state] = {
                "numbers": sorted_nodes,
                "colors": sorted_colors,
                "states": states
            }
        else:
            states_per_graph[state] = info

    return states_per_graph


def generate_graphs(pairs):
    """
    Processes base pairs to generate graphs, checks for bipartiteness, assigns colors, and maps possible ACGU states.

    Args:
        pairs (list): A list of tuples representing base pairs.

    Returns:
        list: A sorted list of dictionaries, each containing the numbers and possible states for each graph, with the graphs sorted by the first number in each graph.
    """

    distinct_graphs = find_distinct_graphs(pairs)

    results = {}
    for index, graph_set in enumerate(distinct_graphs, start=1):
        graph_dict = {node: [] for node in graph_set}
        for a, b in pairs:
            if a in graph_set and b in graph_set:
                graph_dict[a].append(b)
                graph_dict[b].append(a)
        is_bipartite, coloring = is_bipartite_and_color(graph_dict)
        if is_bipartite == False:
            print("The structural restraints in the alternative structures cannot be solved. Check your input.")
            sys.exit()
        simple_coloring = {node: coloring[node] for node in sorted(coloring)}
        results[f"Graph{index}"] = {
            "is_bipartite": is_bipartite,
            "coloring": simple_coloring if is_bipartite else "Not bipartite"
        }

    acgu_states = generate_acgu_states({k: v['coloring'] for k, v in results.items() if v['is_bipartite']})

    for graph_name, states in acgu_states.items():
        results[graph_name].update({"ACGU_states": states})

    formatted_results = []

    for graph_name, graph_info in results.items():
        acgu_states = graph_info['ACGU_states']

        formatted_result = {"numbers": acgu_states['numbers'], "states": acgu_states['states']}

        formatted_results.append(formatted_result)

        formatted_graphs = sorted(formatted_results, key=lambda x: x["numbers"][0])

    return formatted_graphs


def filter_states_by_constraints(graph_numbers, possible_states, sequence_restraint):
    """
    Filters possible ACGU states for a graph based on sequence constraints.

    Args:
        graph_numbers (list): A list of integers representing vertices in the graph.
        possible_states (list): A list of lists, each representing a possible state mapping for the graph vertices.
        sequence_restraint (str): A string representing sequence restraints.

    Returns:
        list: A list of valid states that conform to the sequence restraints.
    """

    restraint_allowed = [constraints_dict[nt] for nt in sequence_restraint]

    valid_states = []
    for state in possible_states:
        is_valid = True
        for i, nt in enumerate(graph_numbers):
            if state[i] not in restraint_allowed[nt]:
                is_valid = False
                break
        if is_valid:
            valid_states.append(state)

    return valid_states


def update_graphs(input_file):
    """
    Updates the possible ACGU states for each graph in the input_file based on sequence constraints.

    Args:
        input_file (InputFile): The input file object containing graph and sequence information.

    Returns:
        None: The function directly modifies the input_file object, updating the states for each graph.
    """

    for i in range(len(input_file.graphs)):
        valid_states = filter_states_by_constraints(input_file.graphs[i]["numbers"], input_file.graphs[i]["states"], input_file.seq_restr)
        if valid_states == []:
            print("The structural restraints are contradictive to sequence restraints. Check your input.")
            sys.exit()
        input_file.graphs[i]["states"] = valid_states


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


def get_allsnakes(input_file):
    """
    Extracts all graph vertices from the input_file and stores them as a flat list in the 'allsnakes' attribute.

    Args:
        input_file (InputFile): The input file object containing graphs.

    Returns:
        None: The function directly modifies the input_file object, adding an 'allsnakes' attribute containing all vertices.
    """

    allsnakes = []

    for i in range(len(input_file.graphs)):
        allsnakes.append(input_file.graphs[i]["numbers"])

    allsnakes_flat = [item for sublist in allsnakes for item in sublist]
    input_file.allsnakes = allsnakes_flat


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

    if input_file.alt_sec_structs != None:
        get_allsnakes(input_file)

    if input_file.excluded_alt_pairs != None:
        pair_list = pair_list + input_file.excluded_alt_pairs
        seen = set()
        index = 0
        while index < len(pair_list):
            element = tuple(pair_list[index])
            if element in seen:
                pair_list.pop(index)
            else:
                seen.add(element)
                index += 1
        input_file.pairs = pair_list

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

    if input_file.graphs != None:
        for i in range(len(input_file.graphs)):
            graph = input_file.graphs[i]
            nts = graph["numbers"]
            states = graph["states"]
            for k in range(len(nts)):
                list_of_nt[nts[k]].snake = True
                list_of_nt[nts[k]].snake_number = i
                list_of_nt[nts[k]].snake_nts = nts
                list_of_nt[nts[k]].snake_states = states

    return list_of_nt


def nt_dictionary(nt):
    """
    Maps a nucleotide character to a list of nucleotide characters based on the IUPAC convention.

    Parameters:
    nt (str): The nucleotide character.

    Returns:
    list: A list of nucleotide characters.
    """

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
    Generate a random RNA sequence that complies with specified sequence restraints.

    Args:
    nt_list (list): A list of Nucleotide objects representing sequence restraints.
    input_file (InputFile): An InputFile object containing input data.
    sim_options (DesignOptions): A SimulationOptions object with simulation settings.

    Returns:
    str: The randomly generated RNA sequence that complies with the specified restraints.
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
    Generate the initial RNA sequence based on the secondary structure information.

    The function iterates over the secondary structure of the RNA (input.sec_struct) and generates
    a corresponding sequence based on the type of structure at each position (e.g., pair, unpair).

    Args:
        nt_list (list): A list of Nucleotide objects representing sequence restraints.
        input_file (InputFile): An InputFile object containing the secondary structure of the RNA.
        sim_options (DesignOptions): A DesignOptions object with simulation settings.

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

    if input_file.seed_seq:
        result_sequence = input_file.seed_seq

    if input_file.graphs != None:
        first_nt_of_snakes = []
        for i in range(len(input_file.graphs)):
            first_nt_of_snakes.append(input_file.graphs[i]["numbers"][0])

        for i in range(len(first_nt_of_snakes)):
            nt_in_snake = nt_list[first_nt_of_snakes[i]].snake_nts
            state_in_snake = nt_list[first_nt_of_snakes[i]].snake_states[0]

            for k in range(len(nt_in_snake)):
                seq_l[nt_in_snake[k]] = state_in_snake[k]
    result_sequence = ''.join(seq_l)

    return result_sequence


def get_alt_mcc(simulation_data, input_file):
    """
    Calculates the MCC for alternative structures against suboptimal structures generated from sequences.

    Args:
        simulation_data (list): A list of dictionaries, each representing a sequence and its simulation data.
        input_file (InputFile): The input file object containing alternative secondary structures.

    Returns:
        list: The updated simulation_data list with added MCC values for each alternative structure.
    """

    for i in (range(len(simulation_data))):
        sequence = simulation_data[i]["sequence"]
        fc = RNA.fold_compound(sequence)
        for k in range(len(input_file.alt_sec_structs)):

            subopt_struct = es.get_first_suboptimal_structure_and_energy(sequence, fc, k + 1)[0]
            target_struct = input_file.alt_sec_structs[k]
            ssc = SimScore(target_struct.replace("&", "Ee"), subopt_struct.replace("&", "Ee"))
            ssc.find_basepairs()
            ssc.cofusion_matrix()

            mcc = ssc.mcc()
            simulation_data[i]["mcc_" + str(k + 1)] = 1 - mcc
            simulation_data[i]["alt_struct_" + str(k + 1)] = subopt_struct

    return simulation_data


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


def get_rep_temps(sim_options):
    """
    Generate temperature shelves for each replica in the simulation.

    The function generates a list of temperatures using a geometric progression,
    with the first temperature being T_min and the last being T_max. The generated
    list of temperatures is returned.

    Args:
        sim_options (DesignOptions): A DesignOptions object with simulation settings.

    Returns:
        list: A list of temperatures for each replica.
    """

    if sim_options.replicas != 1:
        delta = (sim_options.T_max - sim_options.T_min) / (sim_options.replicas - 1)
    else:
        delta = (sim_options.T_max - sim_options.T_min) / 2
    rep_temps = []

    T_curr = sim_options.T_min

    for i in range(sim_options.replicas):
        if sim_options.replicas != 1:
            rep_temps.append(round(T_curr, 3))
            T_curr += delta
        else:
            T_curr += delta
            rep_temps.append(round(T_curr, 3))

    if sim_options.replicas == 1:
        rep_temps = [sim_options.T_max]

    # this piece of code is experimental, it may be used to generate non equally distributed temp shelves
    pot = 1.5
    rep_temps_mod = []

    for i in range(0, len(rep_temps)):
        if i == 0:
            rep_temps_mod.append(rep_temps[i])
        else:
            rep_temps_mod.append((1 - (1 - (rep_temps[i] / sim_options.T_max)**pot)**(1 / pot)) * sim_options.T_max)

    una = False
    if una == True:
        rep_temps = rep_temps_mod

    return rep_temps


def generate_initial_list(nt_list, input_file, sim_options):
    """
    Generate an initial list of RNA sequences based on the input nucleotide list and input file.

    Args:
        nt_list (list): A list of Nucleotide objects.
        input_file (InputFile): An InputFile object containing data from the input file.
        sim_options (DesignOptions): A DesignOptions object with simulation settings.

    Returns:
        list: A list of initialized sequence objects.
    """

    sequence = initial_sequence_generator(nt_list, input_file, sim_options)

    seq_list = []

    for i in range(0, sim_options.replicas):
        if sim_options.diff_start_replicas == "different":
            sequence = initial_sequence_generator(nt_list, input_file, sim_options)
        sequence_object = es.score_sequence(sequence, input_file, sim_options)
        sequence_object.get_replica_num(i + 1)
        sequence_object.get_temp_shelf(sim_options.rep_temps_shelfs[i])
        sequence_object.get_sim_step(0)
        seq_list.append(sequence_object)

    return seq_list


def generate_initial_list_random(nt_list, input_file, sim_options):
    """
    Generate an initial list of random RNA sequences based on the input nucleotide list and input file.

    Args:
        nt_list (list): A list of Nucleotide objects.
        input_file (InputFile): An InputFile object containing data from the input file.
        sim_options (DesignOptions): A DesignOptions object with simulation settings.

    Returns:
        list: A list of initialized random sequence objects.

    This function generates random RNA sequences based on the provided nucleotide list and input file.
    It creates a list of random sequences for the specified number of replicas and calculates scores
    for each sequence using the provided simulation options. The sequences and their details are
    saved to a CSV file with a name based on the simulation options.
    """

    seq_list = []
    for i in range(0, sim_options.replicas):
        sequence = random_sequence_generator(nt_list, input_file, sim_options)
        sequence_object = es.score_sequence(sequence, input_file, sim_options)
        sequence_object.get_replica_num(i + 1)
        sequence_object.get_temp_shelf(sim_options.rep_temps_shelfs[i])
        sequence_object.get_sim_step(0)
        seq_list.append(sequence_object)

    rand_sequences_txt = ""
    for i in range(0, len(seq_list)):
        rand_sequences_txt += str(vars(seq_list[i])) + "\n"

    with open(sim_options.outname + '_random.csv', 'w', newline='\n', encoding='utf-8') as myfile:
        myfile.write(rand_sequences_txt)


def get_mutation_position(seq_obj, available_positions, sim_options, input_file):
    """
    Determines a position in the sequence for mutation.

    Args:
        seq_obj (ScoreSeq): A ScoreSeq object representing an RNA sequence.
        available_positions (list): A list of positions eligible for mutation.
        sim_options (DesignOptions): A DesignOptions object with simulation settings.
        input_file (InputFile): An InputFile object containing data from the input file.

    Returns:
        int: The position in the sequence selected for mutation.

    This function determines a position in the RNA sequence for mutation based on the provided
    parameters and simulation settings. It considers various factors such as point mutations,
    temperature shelves, and structural information to make the selection.
    """

    if sim_options.point_mutations == "off":
        mutation_position = random.choice(available_positions)
    elif sim_options.point_mutations == "on":

        max_perc_prob = sim_options.tm_max
        min_perc_prob = sim_options.tm_min

        pair_list_mfe = check_dot_bracket(seq_obj.mfe_ss)

        query_structure = {tuple(pair) for pair in pair_list_mfe}

        false_negatives = input_file.target_pairs_tupl - query_structure
        false_negatives = [item for tup in false_negatives for item in tup]
        false_negatives = [item for item in false_negatives if item in available_positions]

        false_positives = query_structure - input_file.target_pairs_tupl
        false_positives = [item for tup in false_positives for item in tup]
        false_positives = [item for item in false_positives if item in available_positions]

        false_cases = false_negatives + false_positives
        false_cases = list(set(false_cases))

        prob_shelfs = [round(i, 2) for i in np.linspace(max_perc_prob, min_perc_prob, num=len(sim_options.rep_temps_shelfs))]

        shelf = sim_options.rep_temps_shelfs.index(seq_obj.temp_shelf)

        mutat_point_prob = prob_shelfs[shelf]

        if len(false_cases) == 0:
            range_pos = available_positions
        else:
            false_cases = expand_cases(false_cases, len(seq_obj.sequence) - 1)
            range_pos = random.choices([false_cases, available_positions], weights=[mutat_point_prob, 1 - mutat_point_prob])[0]

        mutation_position = random.choice(range_pos)

    return mutation_position


def expand_cases(cases, max_value, range_expansion=3):
    """
    Expands a set of cases within a specified range without exceeding the maximum value.

    Parameters:
    cases (set): A set of initial cases.
    max_value (int): The maximum allowable value in the range.
    range_expansion (int): The range expansion value around each case.

    Returns:
    list: A sorted list of expanded cases.
    """

    expanded_cases = set()  # Using a set to avoid duplicates

    for case in cases:
        # Adding the cases in the range, taking care not to exceed the maximum value
        for i in range(-range_expansion, range_expansion + 1):
            expanded_value = case + i
            if 0 < expanded_value <= max_value:  # Checking the boundaries
                expanded_cases.add(expanded_value)

    return sorted(list(expanded_cases))


def mutate_sequence(sequence_obj, nt_list, sim_options, input_file):
    """
    Mutate the given RNA sequence and return the mutated sequence as a ScoreSeq object.

    Args:
        sequence_obj (ScoreSeq): A ScoreSeq object representing an RNA sequence.
        nt_list (list): List of Nucleotide objects.
        sim_options (DesignOptions): A DesignOptions object with simulation settings.
        input_file (InputFile): An InputFile object containing data from the input file.

    Returns:
        ScoreSeq: A ScoreSeq object representing the mutated RNA sequence.

    This function performs mutations on the provided RNA sequence object based on the given
    nucleotide list and simulation settings. It returns a ScoreSeq object representing the
    mutated RNA sequence.
    """

    sequence = sequence_obj.sequence
    sequence_list = list(sequence)
    range_pos = []
    for i in range(0, len(sequence)):
        if len(nt_list[i].letters_allowed) != 1:
            range_pos.append(i)

    nt_pos = get_mutation_position(sequence_obj, range_pos, sim_options, input_file)
    nt2_pos = None

    if nt_list[nt_pos].snake == False:
        if (nt_list[nt_pos].pairs_with == None) and (len(nt_list[nt_pos].letters_allowed) != 1):
            available_mutations = nt_list[nt_pos].letters_allowed.copy()
            if len(available_mutations) > 1:
                if sequence_list[nt_pos] in available_mutations:
                    available_mutations.remove(sequence_list[nt_pos])
                mutated_nt = random.choice(available_mutations)
                sequence_list[nt_pos] = mutated_nt
            else:
                print("cannot mutate")
        elif (nt_list[nt_pos].pairs_with == None) and (len(nt_list[nt_pos].letters_allowed) == 1):
            # No mutation needed if there's only one letter allowed and it's not paired
            pass
        elif nt_list[nt_pos].pairs_with != None:
            nt1 = nt_list[nt_pos]

            available_mutations_nt1 = nt1.letters_allowed.copy()

            if (sequence_list[nt_pos] in available_mutations_nt1) and len(available_mutations_nt1) != 1:

                available_mutations_nt1.remove(sequence_list[nt_pos])
            available_mutations_nt1.sort()

            nt2 = nt_list[nt1.pairs_with]
            allowed_mutations_nt2 = nt2.letters_allowed.copy()
            nt2_pos = nt2.number
            if sim_options.acgu_percentages == "on":
                allowed_choices_nt1 = allowed_choice(available_mutations_nt1, sim_options.nt_percentages)
                mutated_nt1 = random.choices(available_mutations_nt1, weights=allowed_choices_nt1)[0]
            elif sim_options.acgu_percentages == "off":
                mutated_nt1 = random.choice(available_mutations_nt1)

            allowed_pairings_nt2 = can_pair(mutated_nt1)
            available_mutations_nt2 = list(set(allowed_mutations_nt2).intersection(allowed_pairings_nt2))

            allowed_mutations_nt2.sort()

            if sim_options.acgu_percentages == "on":
                allowed_choices_nt2 = allowed_choice(available_mutations_nt2, sim_options.nt_percentages)
                mutated_nt2 = random.choices(available_mutations_nt2, weights=allowed_choices_nt2)[0]
            elif sim_options.acgu_percentages == "off":
                mutated_nt2 = random.choice(available_mutations_nt2)

            sequence_list[nt_pos] = mutated_nt1
            sequence_list[nt2_pos] = mutated_nt2
    elif nt_list[nt_pos].snake == True:
        graph_nts = nt_list[nt_pos].snake_nts
        possible_states = nt_list[nt_pos].snake_states
        current_nt_state = sequence_list[nt_pos]

        nt_index = graph_nts.index(nt_pos)

        current_state = next((state for state in possible_states if state[nt_index] == current_nt_state), None)
        new_states = [state for state in possible_states if state != current_state]

        if new_states != []:
            new_state = random.choice(new_states)
            for i in range(len(graph_nts)):
                sequence_list[graph_nts[i]] = new_state[i]

    sequence_mutated = ''.join(sequence_list)
    mut_position = list(' ' * len(sequence))
    mut_position[nt_pos] = "#"

    if nt2_pos != None:
        mut_position[nt2_pos] = "#"
    
    
    if sim_options.oligo_state == "homodimer":
        ss1_inp, ss2_inp =  input_file.sec_struct.split("&")


    if (sim_options.oligo_state == "homodimer") and (nt2_pos != None) and (ss1_inp != ss2_inp) :
        seq1o, seq2o = sequence.split("&")
        seq1m, seq2m = sequence_mutated.split("&")

        diff_positions = sorted([nt_pos, nt2_pos])
        seq1_diff = diff_positions[0]
        seq2_diff = diff_positions[1] - len(seq1o) - 1

        seq1corr = seq1m[:seq2_diff] + seq2m[seq2_diff] + seq1m[seq2_diff + 1:]
        seq2corr = seq2m[:seq1_diff] + seq1m[seq1_diff] + seq2m[seq1_diff + 1:]

        sequence_mutated = seq1corr + "&" + seq2corr
    
    if (sim_options.oligo_state == "homodimer") and (ss1_inp == ss2_inp):

        seq1_old, seq2_old = sequence.split("&")
        seq1m, seq2m = sequence_mutated.split("&")

        if seq1m != seq1_old:
            sequence_mutated = seq1m + "&" + seq1m
        elif seq2m != seq2_old:
            sequence_mutated = seq2m + "&" + seq2m


    sequence_mutated = es.score_sequence(sequence_mutated, input_file, sim_options)
    sequence_mutated.get_replica_num(sequence_obj.replica_num)
    sequence_mutated.get_temp_shelf(sequence_obj.temp_shelf)

    return sequence_mutated


def find_diff_position(str1, str2):
    """
    Finds the first position where two strings differ.

    This function iterates through the characters of two strings of equal length and identifies the
    first index at which the characters differ. If the strings have different lengths or are identical,
    the function returns a specific message or -1, respectively.

    Args:
        str1 (str): The first string for comparison.
        str2 (str): The second string for comparison.

    Returns:
        int: The index of the first differing character.
        str: A message indicating if the strings are of different lengths or are identical.
    """

    if len(str1) != len(str2):
        return "The strings are of different lengths."

    for i in range(len(str1)):
        if str1[i] != str2[i]:
            return i

    return -1


def get_pk_struct(seq, ss_nopk, fc):
    """
    Get the pseudoknot structure of the sequence.

    Args:
    seq (str): RNA sequence.
    ss_nopk (str): Secondary structure of the RNA sequence without pseudoknots.
    fc (FoldCompound): RNA fold compound object.

    Returns:
    str: Secondary structure of the RNA sequence with pseudoknots.
    """

    constraints = ss_nopk.replace("(", "x").replace(")", "x")

    fc.hc_add_from_db(constraints)

    mfe_structure, mfe_energy = fc.mfe()

    l_ss_nopk = list(ss_nopk)
    l_mfe = list(mfe_structure)

    for i in range(0, len(l_ss_nopk)):
        if l_mfe[i] == "(":
            l_ss_nopk[i] = "["
        if l_mfe[i] == ")":
            l_ss_nopk[i] = "]"

    ss_pk = ''.join(l_ss_nopk)

    if "(" in mfe_structure:
        constraints = ss_pk.replace("(", "x").replace(")", "x").replace("[", "x").replace("]", "x")
        fc.hc_add_from_db(constraints)

        mfe_structure, mfe_energy = fc.mfe()
        l_ss_pk = list(ss_pk)
        l_mfe_pk2 = list(mfe_structure)

        for i in range(0, len(l_ss_pk)):
            if l_mfe_pk2[i] == "(":
                l_ss_pk[i] = "<"
            if l_mfe_pk2[i] == ")":
                l_ss_pk[i] = ">"

        ss_pk = ''.join(l_ss_pk)

        if "(" in mfe_structure:
            constraints = ss_pk.replace("(", "x").replace(")", "x").replace("[", "x").replace("]", "x").replace("<", "x").replace(">", "x")
            fc.hc_add_from_db(constraints)

            mfe_structure, mfe_energy = fc.mfe()
            l_ss_pk = list(ss_pk)
            l_mfe_pk2 = list(mfe_structure)

            for i in range(0, len(l_ss_pk)):
                if l_mfe_pk2[i] == "(":
                    l_ss_pk[i] = "{"
                if l_mfe_pk2[i] == ")":
                    l_ss_pk[i] = "}"

                ss_pk = ''.join(l_ss_pk)

    return ss_pk


def score_motifs(seq, sim_options):
    """
    Calculates the scoring of specific motifs within a given RNA sequence.

    This function assesses how well the given RNA sequence matches a set of specified motifs and assigns a score based on this assessment. Each motif has a predefined score, and the function calculates the total score by summing the individual motif scores for all motifs found in the sequence.

    Args:
        seq (str): The RNA sequence to be analyzed.
        sim_options (DesignOptions): A DesignOptions object with simulation settings, including motif definitions.

    Returns:
        float: The total score for all motifs found in the sequence.

    This function iterates through the specified motifs in the simulation options and checks if each motif is present in the RNA sequence. If a motif is found, its corresponding score is added to the total motif score. The final total score is returned as a float.
    """

    motif_score = 0
    for motif in sim_options.motifs:
        if sim_options.motifs[motif][0].search(seq):
            motif_score += sim_options.motifs[motif][1]
    return motif_score


def round_floats(obj):
    """
    Round float values in an object to 3 decimal places.

    Parameters
    ----------
    obj : float or dict
        If a float, the float is rounded to 3 decimal places.
        If a dictionary, all float values in the dictionary are rounded to 3 decimal places.

    Returns
    -------
    obj : float or dict
        The input object with float values rounded to 3 decimal places.
    """

    if isinstance(obj, float):
        return round(obj, 3)
    if isinstance(obj, dict):
        return {k: round_floats(v) for k, v in obj.items()}
    return obj


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
        self.snake = False
        self.snake_number = None
        self.snake_nts = None
        self.snake_states = None

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

    def add_allowed_l(self, list_allowed):
        """
        Sets the list of allowed letters for this nucleotide, considering any constraints.

        Parameters:
        list (list): The list of letters that are allowed for this nucleotide.
        """
        self.letters_allowed = list_allowed
