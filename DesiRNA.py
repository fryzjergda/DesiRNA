#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DesiRNA.py

Usage:
Run the script with command-line arguments to specify input files, parameters for sequence design, and configuration options. The script supports various options for detailed control over the RNA design process.

Example:
python DesiRNA.py [options]

Note:
This script uses several external libraries and modules, such as RNA for RNA secondary structure prediction, and pandas for data handling. Ensure these dependencies are installed and properly configured in your Python environment.
"""

import argparse as argparse
import sys
import random
import math
import csv
import time
import re
import os

from datetime import datetime
from pathlib import Path
from shutil import copy
from shutil import move

import numpy as np
import multiprocess as mp
import pandas as pd

import RNA

from utils import functions_classes as func
from utils.SimScore import SimScore
#from utils import sequence_utils as seq_utils
#from DesiRNA.utils import sequence_utils

def print_action_group_help(parser, action_group, include_all=False, print_usage=True):
    """
    Prints help information for a specific action group within an argparse parser.

    Parameters:
    parser (ArgumentParser): The argparse parser object.
    action_group (ArgumentGroup): The action group within the parser for which to print help.
    include_all (bool): If True, includes all actions in the usage. Defaults to False.
    print_usage (bool): If True, prints the usage information. Defaults to True.
    """

    formatter = parser._get_formatter()

    # Add usage for the specific action group
    if print_usage:
        actions = action_group._group_actions if not include_all else parser._actions
        formatter.add_usage(parser.usage, actions, parser._mutually_exclusive_groups)

    formatter.start_section(action_group.title)
    formatter.add_text(action_group.description)
    formatter.add_arguments(action_group._group_actions)
    formatter.end_section()
    print(formatter.format_help())


class StandardHelpAction(argparse.Action):
    """
    A custom argparse action to provide standard help options in a command-line interface.

    This action displays help for standard options when invoked.
    """

    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):
        """
        Initialize the StandardHelpAction.

        Parameters:
        option_strings (list): A list of option strings which should trigger this action.
        dest (str): The attribute to hold the result. Its value is argparse.SUPPRESS by default.
        default: The default value; argparse.SUPPRESS means suppress entirely.
        help (str): The help string for the action.
        """
        super(StandardHelpAction, self).__init__(option_strings=option_strings, dest=dest, default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        # Display only the standard options
        for action_group in parser._action_groups:
            if action_group.title == "Standard Options":
                print_action_group_help(parser, action_group)
        parser.exit()


class AdvancedHelpAction(argparse.Action):
    """
    A custom argparse action to provide advanced help options in a command-line interface.

    This action displays help for both standard and advanced options when invoked.
    """

    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):
        """
        Initialize the AdvancedHelpAction.

        Parameters:
        option_strings (list): A list of option strings which should trigger this action.
        dest (str): The attribute to hold the result. Its value is argparse.SUPPRESS by default.
        default: The default value; argparse.SUPPRESS means suppress entirely.
        help (str): The help string for the action.
        """
        super(AdvancedHelpAction, self).__init__(option_strings=option_strings, dest=dest, default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        # Display both standard and advanced options
        print_action_group_help(parser, parser._action_groups[-2], include_all=True)  # Standard Options with usage
        print_action_group_help(parser, parser._action_groups[-1], include_all=True, print_usage=False)  # Advanced Options without usage
        parser.exit()


def argument_parser():
    """
    Creates and configures an argument parser for the command line interface.

    This function sets up the parser with all necessary options and custom actions.

    Returns:
    ArgumentParser: Configured parser object.
    """

    parser = argparse.ArgumentParser(description=__doc__, prog='DesiRNA.py', add_help=False)

    standard_group = parser.add_argument_group('Standard Options')
    advanced_group = parser.add_argument_group('Advanced Options')

    standard_group.add_argument("-f", "--filename", required=True, dest="name",
                                help="Name of a file that contains secondary structures and constraints.")
    standard_group.add_argument("-R", "--replicas", required=False, dest="replicas", default=10, type=int,
                                help="Number of replicas. [default = 10]")
    standard_group.add_argument("-e", "--exchange", required=False, dest="exchange", default=100, type=int,
                                help="Frequency of replixa exchange attempt. [default = 100]")
    standard_group.add_argument("-t", "--timelimit", required=False, default=60, dest="timlim", type=int,
                                help="Timelimit for running the program [s]. [default = 60]")
    standard_group.add_argument("-s", "--steps", required=False, default=None, dest="steps", type=int,
                                help="Number of Replica Exchange steps after which the simulation ends. Overwrites the -t option. [default = None]")
    standard_group.add_argument("-r", "--results_number", required=False, dest="num_results", default=10, type=int,
                                help="Number of best results to be reported in the output. [default = 10]")
    advanced_group.add_argument("-sws", "--stop_when_solved", required=False, dest="sws", default='off', choices=['off', 'on'],
                                help="Stop after finding desired number of solutions (--results_number). [default = off]")
    advanced_group.add_argument("-p", "--param", required=False, dest="param", default='1999', choices=['2004', '1999'],
                                help="Turner energy parameter for calculating MFE. [default = 1999]")
    advanced_group.add_argument("-tmin", "--tmin", required=False, dest="t_min", default=10, type=float,
                                help="Minimal Replica Temperature. [default = 10]")
    advanced_group.add_argument("-tmax", "--tmax", required=False, dest="t_max", default=150, type=float,
                                help="Maximal Replica Temperature. [default = 150]")
    advanced_group.add_argument("-ts", "--tshelves", required=False, dest="tshelves", type=str, default='',
                                help="Custom temperature shelves for replicas in replica exchange simulation. Provide comma-separated values.")
    advanced_group.add_argument("-sf", "--scoring_function", required=False, dest="scoring_f", type=str, default='Ed-Epf:1.0',
                                help="Scoring functions and weights used to guide the design process, e.g. 'Ed-Epf:0.5,1-MCC:0.5'. \
                                Scoring functions to choose: Ed-Epf, 1-MCC, sln_Epf, Ed-MFE, 1-precision, 1-recall, Edef [default = Ed-Epf:1.0]")
    advanced_group.add_argument("-nd", "--negative_design", required=False, dest="subopt", default='off', choices=['off', 'on'],
                                help="Use negative design approach. [default = off]")
    standard_group.add_argument("-acgu", "--ACGU", required=False, dest="percs", default='off', choices=['off', 'on'],
                                help="Keep 'natural' ACGU content. If turned on the content will be A:15%%, C:30%%, G:30%%, U:15%%. [default = off]")
    advanced_group.add_argument("-acgu_content", "--ACGU_content", required=False, dest="acgu_content", default='', type=str,
                                help="Provide user defined ACGU content. Comma-separated values e.g., 10,40,40,10")
    advanced_group.add_argument("-o", "--avoid_oligomerization", required=False, dest="oligo", default='off', choices=['off', 'on'],
                                help="Designes sequneces that should not tend to oligomerize. Slows down the simulation. [default = off]")
    advanced_group.add_argument("-d", "--dimer", required=False, dest="dimer", default='off', choices=['off', 'on'],
                                help="Design of a homodimer complex, of two strands. [deafult = off, requires input file complying with RNA-RNA complex format]")
    advanced_group.add_argument("-tm", "--target_mutations", required=False, dest="pm", default='on', choices=['off', 'on'],
                                help="Targeted mutations. Targets mostly False Negativeas and False Positives. [default = on]")
    advanced_group.add_argument("-tm_perc_max", "--target_mutations_percentage_max", required=False, dest="tm_max", default=0.7, type=float,
                                help="Highest percentage of targeted mutations applied to lowest temperature replica. Percentage for replicas in between will be set evenly from 'tm_perc_max' to 'tm_perc_min'. Float from 0.0 to 1.0. [default = 0.7]")
    advanced_group.add_argument("-tm_perc_min", "--target_mutations_percentage_min", required=False, dest="tm_min", default=0.0, type=float,
                                help="Lowest percentage of targeted mutations applied to highest temperature replica. Percentage for replicas in between will be set evenly from 'tm_perc_max' to 'tm_perc_min'. Float from 0.0 to 1.0. [default = 0.0]")
    advanced_group.add_argument("-motifs", "--motif_sequences", required=False, dest="motifs", type=str, default='',
                                help="Sequence motifs along with their bonuses(-)/penalties(+). Provide comma-separated key,value,key,value sequence.")
    advanced_group.add_argument("-seed", "--seed_number", required=False, default=0, dest="in_seed", type=int,
                                help="User defined seed number for simulation. [default = 0]")
    advanced_group.add_argument("-re_seq", "--replicas_sequences", required=False, dest="diff_start_replicas", default='one', choices=['different', 'same'],
                                help="Choose wether replicas will start from the same random sequence or each replica will start from different random sequence. [default = same]")

    parser.add_argument('-h', action=StandardHelpAction, help='Show standard help message and exit.')
    parser.add_argument('-H', action=AdvancedHelpAction, help='Show advanced help message and exit.')

    args = parser.parse_args()

    if args.tshelves:
        temperature_shelves = [float(temp) for temp in args.tshelves.split(",")]

        # Check if the number of temperatures matches the number of replicas
        if len(temperature_shelves) != args.replicas:
            parser.error("The number of temperatures provided in -ts/--tshelves must match the number of replicas set by -R.")

    if args.motifs:
        iupac_re = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'U', 'W': '[AU]', 'S': '[GC]', 'M': '[AC]',
                    'K': '[GU]', 'R': '[AG]', 'Y': '[CU]', 'B': '[CGU]', 'D': '[AGU]', 'H': '[ACU]', 'V': '[ACG]', 'N': '[ACGU]', }

        motifs_list = args.motifs.split(',')
        try:
            motifs = {motifs_list[i]: (re.compile(''.join([iupac_re[ch]
                                                          for ch in motifs_list[i]])),
                                       float(motifs_list[i + 1])) for i in range(0, len(motifs_list), 2)}
        except:
            parser.error("Something is wrong with your motifs")
    else:
        motifs = {}

    if args.percs == 'off' and args.acgu_content:
        parser.error("The -acgu_content option is only applicable when -acgu is set to 'on'.")

    if args.pm == 'off' and (args.tm_max != 0.7 or args.tm_min != 0.0):
        parser.error("Targeted mutations must be 'on' to set -tm_perc_max and/or -tm_perc_min.")

    infile = args.name
    replicas = args.replicas
    timlim = args.timlim
    acgu_percs = args.percs
    t_max = args.t_max
    t_min = args.t_min
    oligo = args.oligo
    dimer = args.dimer
    param = args.param
    exchange_rate = args.exchange
    scoring_f = parse_scoring_functions(args.scoring_f)
    pm = args.pm
    steps = args.steps
    tshelves = args.tshelves
    in_seed = args.in_seed
    subopt = args.subopt
    diff_start_replicas = args.diff_start_replicas
    num_results = args.num_results
    acgu_content = args.acgu_content
    tm_max = args.tm_max
    tm_min = args.tm_min
    sws = args.sws
    
    return infile, replicas, timlim, acgu_percs, t_max, t_min, oligo, param, exchange_rate, scoring_f, pm, tshelves,\
        in_seed, subopt, diff_start_replicas, num_results, acgu_content, steps, tm_max, tm_min, motifs, dimer, sws


def read_input(infile):
    """
    Reads the input file and returns an InputFile object.

    Parameters:
    infile (str): Path to the input file.

    Returns:
    InputFile: An object containing data from the input file.
    """

    string = ""
    with open(infile, encoding='utf-8') as f:
        for line in f:
            string += line

    data = string.lstrip(">").rstrip("\n").split("\n>")
    data_dict = {}
    for i in range(0, len(data)):
        type_entry = data[i].split("\n")
        key, values = type_entry[0], type_entry[1:]
        data_dict[key] = values

    input_file = InputFile(data_dict['name'][0], data_dict['sec_struct'][0], data_dict['seq_restr'][0])

    if 'seed_seq' in data_dict:
        input_file.add_seed_seq(data_dict['seed_seq'][0])

    if 'alt_sec_struct' in data_dict:
        input_file.add_alt_sec_struct(data_dict['alt_sec_struct'])

    return input_file


def parse_scoring_functions(scoring_f_str):
    """
    Parses the scoring function string into a list of scoring function tuples.

    Parameters:
    scoring_f_str (str): A string representing scoring functions and their weights.

    Returns:
    list: A list of tuples where each tuple contains a scoring function and its weight.

    Raises:
    ValueError: If the scoring function format is invalid.
    """

    scoring_func = []
    for item in scoring_f_str.split(','):
        if ':' not in item:
            raise ValueError(f"Invalid scoring function format: {item}. Expected format: 'function:weight'")
        function, weight = item.split(':')
        scoring_func.append((function, float(weight)))
    return scoring_func


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


def random_sequence_generator(nt_list, input_file):
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

        allowed_choices = allowed_choice(nt1, nt_percentages)

        seq_l[nt1num] = random.choices(nt1, weights=allowed_choices)[0]
        seq_l[nt2num] = wc_pair(seq_l[nt1num])

    for i in range(0, len(seq_l)):
        if seq_l[i] not in ["A", "C", "G", "U"]:
            seq_l[i] = random.choice(nt_dictionary(seq_l[i]))

    result_sequence = ''.join(seq_l)

    return result_sequence


def initial_sequence_generator(nt_list, input_file):
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

    if acgu_percentages == 'off':

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

    elif acgu_percentages == 'on':
        for i in range(0, len(pair_list)):
            nt1 = nt_list[pair_list[i][0]].letters_allowed
            nt2 = nt_list[pair_list[i][1]].letters_allowed
            nt1.sort()
            nt2.sort()
            nt1num = nt_list[pair_list[i][0]].number
            nt2num = nt_list[pair_list[i][1]].number

            allowed_choices = allowed_choice(nt1, nt_percentages)

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


def get_rep_temps():
    """
    Generates the temperature shelves for each replica in the simulation.

    The function generates a list of temperatures using a geometric progression,
    with the first temperature being T_min and the last being T_max. The generated
    list of temperatures is returned.

    Returns:
        list: List of temperatures for each replica.
    """

    if replicas != 1:
        delta = (T_max - T_min) / (replicas - 1)
    else:
        delta = (T_max - T_min) / 2
    rep_temps = []

    T_curr = T_min

    for i in range(replicas):
        if replicas != 1:
            rep_temps.append(round(T_curr, 3))
            T_curr += delta
        else:
            T_curr += delta
            rep_temps.append(round(T_curr, 3))

    if replicas == 1:
        rep_temps = [T_max]

    # this piece of code is experimental, it may be used to generate non equally distributed temp shelves
    pot = 1.5
    rep_temps_mod = []

    for i in range(0, len(rep_temps)):
        if i == 0:
            rep_temps_mod.append(rep_temps[i])
        else:
            rep_temps_mod.append((1 - (1 - (rep_temps[i] / T_max)**pot)**(1 / pot)) * T_max)

    una = False
    if una == True:
        rep_temps = rep_temps_mod

    return rep_temps


def generate_initial_list(nt_list, input_file):
    """
    Generates an initial list of RNA sequences based on the input nucleotide list and input file.

    Parameters:
    nt_list (list): A list of Nucleotide objects.
    input_file (InputFile): The InputFile object containing data from the input file.

    Returns:
    list: A list of initialized sequence objects.
    """

    sequence = initial_sequence_generator(nt_list, input_file)

    seq_list = []

    for i in range(0, replicas):
        if diff_start_replicas == "different":
            sequence = initial_sequence_generator(nt_list, input_file)
        sequence_object = score_sequence(sequence)
        sequence_object.get_replica_num(i + 1)
        sequence_object.get_temp_shelf(rep_temps_shelfs[i])
        sequence_object.get_sim_step(0)
        seq_list.append(sequence_object)

    return seq_list


def generate_initial_list_random(nt_list, input_file):
    """
    Generate an initial list of random RNA sequences based on the input nucleotide list and input file.

    Parameters:
    nt_list (list): A list of Nucleotide objects.
    input_file (InputFile): The InputFile object containing data from the input file.

    Returns:
    list: A list of initialized random sequence objects.
    """

    seq_list = []
    for i in range(0, replicas):
        sequence = random_sequence_generator(nt_list, input_file)
        sequence_object = score_sequence(sequence)
        sequence_object.get_replica_num(i + 1)
        sequence_object.get_temp_shelf(rep_temps_shelfs[i])
        sequence_object.get_sim_step(0)
        seq_list.append(sequence_object)

    rand_sequences_txt = ""
    for i in range(0, len(seq_list)):
        rand_sequences_txt += str(vars(seq_list[i])) + "\n"

    with open(outname + '_random.csv', 'w', newline='\n', encoding='utf-8') as myfile:
        myfile.write(rand_sequences_txt)


def get_mutation_position(seq_obj, available_positions):
    """
    Determines a position in the sequence to mutate.

    Parameters:
    sequence_obj (ScoreSeq): A ScoreSeq object representing an RNA sequence.
    range_pos (list): A list of positions eligible for mutation.

    Returns:
    int: The position in the sequence selected for mutation.
    """

    if point_mutations == "off":
        mutation_position = random.choice(available_positions)
    elif point_mutations == "on":

        max_perc_prob = tm_max
        min_perc_prob = tm_min

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

        prob_shelfs = [round(i, 2) for i in np.linspace(max_perc_prob, min_perc_prob, num=len(rep_temps_shelfs))]

        shelf = rep_temps_shelfs.index(seq_obj.temp_shelf)

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


def mutate_sequence(sequence_obj, nt_list):
    """
    Mutate the given sequence and return the mutated sequence.

    Args:
    sequence_obj (ScoreSeq): A ScoreSeq object representing an RNA sequence.
    nt_list (list): List of nucleotide objects.

    Returns:
    ScoreSeq: A ScoreSeq object representing the mutated RNA sequence.
    """

    sequence = sequence_obj.sequence
    sequence_list = list(sequence)
    range_pos = []
    for i in range(0, len(sequence)):
        if len(nt_list[i].letters_allowed) != 1:
            range_pos.append(i)

    nt_pos = get_mutation_position(sequence_obj, range_pos)
    nt2_pos = None

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
        if acgu_percentages == "on":
            allowed_choices_nt1 = allowed_choice(available_mutations_nt1, nt_percentages)
            mutated_nt1 = random.choices(available_mutations_nt1, weights=allowed_choices_nt1)[0]
        elif acgu_percentages == "off":
            mutated_nt1 = random.choice(available_mutations_nt1)

        allowed_pairings_nt2 = can_pair(mutated_nt1)
        available_mutations_nt2 = list(set(allowed_mutations_nt2).intersection(allowed_pairings_nt2))

        allowed_mutations_nt2.sort()

        if acgu_percentages == "on":
            allowed_choices_nt2 = allowed_choice(available_mutations_nt2, nt_percentages)
            mutated_nt2 = random.choices(available_mutations_nt2, weights=allowed_choices_nt2)[0]
        elif acgu_percentages == "off":
            mutated_nt2 = random.choice(available_mutations_nt2)

        sequence_list[nt_pos] = mutated_nt1
        sequence_list[nt2_pos] = mutated_nt2

    sequence_mutated = ''.join(sequence_list)
    mut_position = list(' ' * len(sequence))
    mut_position[nt_pos] = "#"

    if nt2_pos != None:
        mut_position[nt2_pos] = "#"

    sequence_mutated = score_sequence(sequence_mutated)
    sequence_mutated.get_replica_num(sequence_obj.replica_num)
    sequence_mutated.get_temp_shelf(sequence_obj.temp_shelf)

    return sequence_mutated


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


def mc_delta(deltaF_o, deltaF_m, T_replica):
    """
    Determine whether a mutation is accepted based on the Metropolis criterion.

    Args:
    T_replica (float): Temperature of the replica.
    deltaF_o (float): Original free energy.
    deltaF_m (float): Mutated free energy.

    Returns:
    tuple: A tuple containing two boolean values indicating whether the mutation is accepted and whether the mutated scoring function is better.
    """

    accept_e = False
    if deltaF_m <= deltaF_o:
        accept = True
        accept_e = True
    else:
        diff = deltaF_m - deltaF_o
        p_dE = metropolis_score(T_replica, diff)
        rand_num = random.random()

        accept = p_dE > rand_num

    return accept, accept_e


def metropolis_score(temp, dE):
    """
    Calculate the Metropolis score.

    Args:
    temp (float): Temperature.
    dE (float): Energy difference.

    Returns:
    float: Metropolis score.
    """

    score = math.exp((-L / temp) * (dE))

    return score


def replica_exchange_attempt(T0, T1, dE0, dE1):
    """
    Attempt to exchange replicas.

    Args:
    T0 (float): Temperature of the first replica.
    T1 (float): Temperature of the second replica.
    dE0 (float): Energy of the first replica.
    dE1 (float): Energy of the second replica.

    Returns:
    tuple: A tuple containing two boolean values indicating whether the exchange is accepted and whether the exchange resulted in a better energy.
    """

    accept_e = False
    if dE1 <= dE0:
        accept = True
        accept_e = True
    else:
        rand_num = random.random()
        p = math.exp(L * (1 / T0 - 1 / T1) * (dE0 - dE1))
        accept = p > rand_num

    return accept, accept_e


def replica_exchange(replicas, stats_obj):
    """
    Perform the replica exchange.

    Args:
        replicas (list): List of replica objects.
        stats_obj (object): Statistics object.

    Returns:
        tuple: A tuple containing the list of replicas after exchange and the updated statistics object.
    """

    re_pairs = []

    num_shelfs = [i for i in range(1, len(replicas))]

    temps = []

    for i in range(0, len(replicas)):
        temps.append(replicas[i].temp_shelf)

    temps = sorted(temps)

    if stats_obj.global_step % 2 == 0:

        for i in range(0, len(num_shelfs) - 1, 2):
            re_pairs.append([i + 1, i + 2])
    else:
        for i in range(0, len(num_shelfs), 2):
            re_pairs.append([i, i + 1])

    replicas = sorted(replicas, key=lambda obj: obj.temp_shelf)

    for i in range(len(re_pairs)):
        rep_i = re_pairs[i][0]  # replica temporary number
        rep_j = re_pairs[i][1]  # replica temporary number

        T_i = replicas[rep_i].temp_shelf
        T_j = replicas[rep_j].temp_shelf

        dE_i = replicas[rep_i].scoring_function
        dE_j = replicas[rep_j].scoring_function

        accept, accept_e = replica_exchange_attempt(T_i, T_j, dE_i, dE_j)

        if accept == True:
            replicas[rep_i].get_temp_shelf(T_j)
            replicas[rep_j].get_temp_shelf(T_i)
            stats_obj.update_acc_re_step()
            if accept_e == True:
                stats_obj.update_acc_re_better_e()
        else:
            stats_obj.update_rej_re_step()

    replicas = sorted(replicas, key=lambda obj: obj.replica_num)

    return replicas, stats_obj


md = RNA.md()
md.compute_bpp = 0


def get_mfe_e_ss(seq):
    """
    This function calculates the minimum free energy (MFE) and ensemble free energy (EFE)
    of a given RNA sequence.
    Parameters:
    seq (str): The RNA sequence to analyze.
    Returns:
    (float, float): A tuple containing the MFE and EFE of the sequence.
    """
    fc = RNA.fold_compound(seq, md)

    if oligo_state in {"none", "avoid"}:
        structure, energy = fc.pf()
        structure = fc.mfe()[0]
        if pks == "on":
            structure = get_pk_struct(seq, structure, fc)
    elif oligo_state in {"homodimer", "heterodimer"}:
        seqa_len = len(seq.split("&")[0])
        structure_dim, energy = fc.mfe_dimer()
        structure = structure_dim[:seqa_len] + "&" + structure_dim[seqa_len:]
    return energy, structure, fc


def score_motifs(seq, motifs):
    """
    Calculates the scoring of specific motifs within a given RNA sequence.

    This function assesses how well the given RNA sequence matches a set of specified motifs and assigns a score based on this assessment. Each motif has a predefined score, and the function calculates the total score for all motifs present in the sequence.

    Parameters:
    seq (str): The RNA sequence to be analyzed.
    motifs (dict): A dictionary where keys are motifs (as strings) and values are their respective scores.

    Returns:
    float: The total score for all motifs found in the sequence. This score is a sum of individual motif scores for all motifs present in the sequence.
    """

    motif_score = 0
    for motif in motifs:
        if motifs[motif][0].search(seq):
            motif_score += motifs[motif][1]
    return motif_score


def score_sequence(seq):
    """
    Calculate various scoring metrics for a given RNA sequence.

    Parameters
    ----------
    seq : str
        RNA sequence to score.

    Returns
    -------
    scored_sequence : ScoreSeq object
        Object with various scoring metrics, including mfe, target energy,
        energy difference, precision, recall, mcc, and overall score.
        If the 'oligo' option is set, the oligomerization of the sequence
        is also calculated.

    Notes
    -----
    This function uses the ViennaRNA package to calculate the mfe
    and the energy of the target structure. The precision, recall,
    and MCC are calculated using the SimScore class. The overall score
    is a combination of these metrics.
    """

    if oligo_state != "homodimer":
        scored_sequence = func.ScoreSeq(sequence=seq)
    else:
        seq = seq.split("&")[0]
        seq = seq + "&" + seq
        scored_sequence = func.ScoreSeq(sequence=seq)
    pf_energy, mfe_structure, fold_comp = get_mfe_e_ss(seq)

    scored_sequence.get_Epf(pf_energy)
    scored_sequence.get_mfe_ss(mfe_structure)

    scored_sequence.get_edesired(fold_comp.eval_structure(input_file.sec_struct.replace("&", "")))

    scored_sequence.get_edesired_minus_Epf(scored_sequence.Epf, scored_sequence.edesired)

    ssc = SimScore(input_file.sec_struct.replace("&", "Ee"), scored_sequence.mfe_ss.replace("&", "Ee"))
    ssc.find_basepairs()
    ssc.cofusion_matrix()

    scored_sequence.get_precision(ssc.precision())
    scored_sequence.get_recall(ssc.recall())
    scored_sequence.get_mcc(ssc.mcc())

    for function, weight in scoring_f:
        if function == 'sln_Epf':
            scored_sequence.get_sln_Epf()
        if function == 'Ed-MFE':
            scored_sequence.get_MFE()
            scored_sequence.get_edesired_minus_MFE()
        if function == 'Edef':
            scored_sequence.get_ensemble_defect(input_file.sec_struct)

    scored_sequence.get_scoring_function(scoring_f)

    if input_file.alt_sec_struct != None:
        energies = [fold_comp.eval_structure(alt_dbn) for alt_dbn in input_file.alt_sec_structs]
        scored_sequence.get_edesired2(sum(energies) / len(energies))
        scored_sequence.get_edesired2_minus_Epf(scored_sequence.Epf, scored_sequence.edesired2)
        scored_sequence.get_scoring_function_w_alt_ss()

    if (subopt == "on") and (scored_sequence.mcc == 0):
        scored_sequence.get_subopt_e(func.get_first_suboptimal_structure_and_energy(seq, fold_comp)[1])
        scored_sequence.get_esubopt_minus_Epf(scored_sequence.Epf, scored_sequence.subopt_e)
        scored_sequence.get_scoring_function_w_subopt()

    if oligo_state == "homodimer" or oligo_state == "heterodimer":
        scored_sequence.get_scoring_function_oligomer(fold_comp)

    if oligo_state == "avoid":
        scored_sequence.get_scoring_function_monomer()

    if motifs:
        motif_score = score_motifs(seq, motifs)
        scored_sequence.update_scoring_function_w_motifs(motif_score)

    return scored_sequence


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
    elif isinstance(obj, dict):
        return {k: round_floats(v) for k, v in obj.items()}
    return obj


def single_replica_design(sequence_o, nt_list, worker_stats):
    """
    Perform the replica design for a single replica.

    Parameters
    ----------
    sequence_o : ScoreSeq object
        The original sequence to be mutated.
    nt_list : list
        List of possible nucleotides for mutation.
    worker_stats : Stats object
        Object to store the statistics of the worker.

    Returns
    -------
    best_score : ScoreSeq object
        The best scored sequence after performing the replica design.
    worker_stats : Stats object
        Updated statistics of the worker.
    """

    worker_stats.reset_mc_stats()

    for i in range(0, RE_attempt):
        sequence_m = mutate_sequence(sequence_o, nt_list)

        deltaF_o = sequence_o.scoring_function
        deltaF_m = sequence_m.scoring_function

        accept, accept_e = mc_delta(deltaF_o, deltaF_m, sequence_o.temp_shelf)

        if accept == True:
            sequence_o = sequence_m
            worker_stats.update_acc_mc_step()
            if accept_e == True:
                worker_stats.update_acc_mc_better_e()

        if accept == False:
            worker_stats.update_rej_mc_step()

    return sequence_o, worker_stats


def par_wrapper(args):
    """
    Wrapper function for parallel processing.

    Parameters
    ----------
    args : tuple
        Arguments to be passed to the function.

    Returns
    -------
    result : various
        The result of the function call.
    """
    args, seed = args[:-1], args[-1]
    random.seed(seed)

    return single_replica_design(*args)


def mutate_sequence_re(lst_seq_obj, nt_list, stats_obj):
    """
    This function mutates a list of sequence objects in parallel.

    Parameters:
    lst_seq_obj (list): A list of sequence objects to mutate.
    nt_list (list): A list of possible nucleotides to use in mutation.
    stats_obj (object): A statistical object to keep track of the mutation process.

    Returns:
    list: A list of mutated sequence objects.
    """

    with mp.Pool(replicas) as pool:
        # Generate a unique seed for each worker based on original_seed
        seeds = [original_seed + i for i in range(len(lst_seq_obj))]

        # inputs = [(seq_obj, nt_list, stats_obj) for seq_obj in lst_seq_obj]
        inputs = [(seq_obj, nt_list, stats_obj, seed) for seq_obj, seed in zip(lst_seq_obj, seeds)]
        lst_seq_obj_res_new = []
        workers_stats = []

        for result in pool.imap(par_wrapper, inputs):
            lst_seq_obj_res_new.append(result[0])
            workers_stats.append(result[1])

        attributes_to_update = ['acc_mc_step', 'acc_mc_better_e', 'rej_mc_step']

        stats_obj.update_step(RE_attempt)

        for stats in workers_stats:
            for attr in attributes_to_update:
                current_value = getattr(stats_obj, attr)
                additional_value = getattr(stats, attr)
                setattr(stats_obj, attr, current_value + additional_value)

        return lst_seq_obj_res_new, stats_obj


def initialize_simulation(input_file):
    """
    Initialize the simulation by setting up the necessary parameters and constraints.

    Args:
    input_file (InputFile): An object encapsulating input file parameters.

    Returns:
    list: List of nucleotides with applied constraints.
    """

    input_file.pairs = check_dot_bracket(input_file.sec_struct)  # check dotbracket correctness, assign as list of pairs
    input_file.set_target_pairs_tupl()

    if input_file.alt_sec_structs != None:
        print(input_file.alt_sec_structs)
        for i in range(0, len(input_file.alt_sec_structs)):
            print("Checking brackets for alternative structure ", i + 1)
            check_dot_bracket(input_file.alt_sec_structs[i])
        print("Alternative structures are ok.")

    check_seq_restr(input_file.seq_restr)
    check_length(input_file.sec_struct, input_file.seq_restr)
    nt_list = get_nt_list(input_file)
    check_input_logic(nt_list)

    return nt_list


def generate_sequences(nt_list, input_file):
    """
    Generates initial RNA sequences based on nucleotide constraints.

    Args:
    nt_list (list): A list of nucleotide constraints.
    input_file (InputFile): The input file object with sequence data.

    Returns:
    list: A list of initialized sequence objects.
    """
    sequence_score_list = generate_initial_list(nt_list, input_file)
    generate_initial_list_random(nt_list, input_file)

    return sequence_score_list


def handle_non_mutable_sequence(input_file, simulation_data, output_name):
    """
    Handles the scenario where the sequence cannot be mutated due to sequence constraints.
    Scores the sequence and writes the results to a CSV file, then exits the program.

    Args:
    input_file (InputFile): The input file object containing sequence data.
    simulation_data (list): The simulation data to be processed.
    output_name (str): The base name for the output file.
    """
    if all(char in "ACGU&" for char in input_file.seq_restr):
        print("Cannot mutate the sequence due to sequence constraints.\nScoring sequence.")
        sorted_results = sorted(round_floats(simulation_data), key=lambda d: (-d['mcc'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
        with open(outname + '_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = sorted_results[0].keys()  # header from keys of the first dictionary
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for data in sorted_results:
                writer.writerow(data)
        sys.exit()


def process_intermediate_results(simulation_data, oligo, output_name, num_results):
    """
    Processes intermediate results during the simulation by removing duplicates, 
    sorting the results, and writing them to a CSV file.

    Args:
    simulation_data (list): The simulation data to be processed.
    oligo (str): The oligo option to determine the sorting criteria.
    output_name (str): The base name for the output file.
    num_results (int): Number of top results to include in the output.

    Returns:
    list: A list of sorted intermediate results.
    """
    remove_duplicates = {item['sequence']: item for item in simulation_data}
    data_no_duplicates = list(remove_duplicates.values())

    if oligo != "off":
        sorted_results = sorted(round_floats(data_no_duplicates), key=lambda d: (-d['oligo_fraction'], -d['mcc'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
    else:
        sorted_results = sorted(round_floats(data_no_duplicates), key=lambda d: (-d['mcc'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)

    sorted_results = sorted_results[:num_results]
    mcc_zero_count = sum(1 for item in sorted_results if item['mcc'] == 0.0)

    with open(output_name + '_mid_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = sorted_results[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for data in sorted_results:
            writer.writerow(data)

    return mcc_zero_count


def generate_replica_csv(simulation_data, output_name):
    """
    Generates a CSV file containing simulation data for each replica, pivoted on specified metrics.

    Args:
    simulation_data (list): The simulation data to be processed.
    output_name (str): Name of the output CSV file.
    """
    df = pd.DataFrame(simulation_data)
    df.sort_values(['sim_step', 'replica_num'], inplace=True)
    df = df.drop_duplicates(subset=['sim_step', 'replica_num'])

    final_df = pd.DataFrame()
    header = ['']
    header2 = ['']

    for metric in ['sequence', 'Epf', 'edesired_minus_Epf', 'subopt_e', 'esubopt_minus_Epf', 'precision', 'recall', 'mcc', 'sln_Epf', 'temp_shelf', 'scoring_function']:
        sub_df = df[['sim_step', 'replica_num', metric]].pivot(index='sim_step', columns='replica_num', values=metric)
        sub_df.columns = [f'replica{col}' for col in sub_df.columns]

        header.extend([metric] * len(sub_df.columns) + [''])
        header2.extend(sub_df.columns.tolist() + [''])

        final_df = pd.concat([final_df, sub_df, pd.DataFrame(columns=[' '])], axis=1)

    with open(output_name + '_replicas.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow(header2)
        final_df.to_csv(csvfile, header=False, index=True)


def sort_trajectory(simulation_data):
    """
    Sorts the simulation data based on simulation step and replica number.

    Args:
    simulation_data (list): The simulation data to be sorted.

    Returns:
    list: Sorted simulation data.
    """
    return sorted(round_floats(simulation_data), key=lambda d: (d['sim_step'], d['replica_num']))


def generate_trajectory_csv(sorted_data, output_name):
    """
    Generates a CSV file from the sorted simulation data.

    Args:
    sorted_data (list): Sorted simulation data.
    output_name (str): Name of the output CSV file.
    """
    with open(output_name + '_traj.csv', 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = sorted_data[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for data in sorted_data:
            writer.writerow(data)


def generate_multifasta(sorted_data, output_name, infile, now):
    """
    Generates a multi-FASTA file from the sorted simulation data.

    Args:
    sorted_data (list): Sorted simulation data.
    output_name (str): Name of the output FASTA file.
    infile (str): Name of the input file.
    now (str): Timestamp or unique identifier.
    """
    fasta_txt = ""
    for item in sorted_data:
        fasta_txt += f">{infile}|{now}|{item['replica_num']}|{item['sim_step']}|{item['scoring_function']}\n{item['sequence']}\n"

    with open(output_name + '_multifasta.fas', 'w', encoding='utf-8') as fastafile:
        fastafile.write(fasta_txt)


def generate_best_fasta(simulation_data, num_results, output_name, infile, now):
    """
    Generates a FASTA file from the best sorted simulation data.

    Args:
    simulation_data (list): The simulation data to be processed.
    num_results (int): Number of top results to include in the FASTA file.
    output_name (str): Name of the output FASTA file.
    infile (str): Name of the input file.
    now (str): Timestamp or identifier for the sequence.
    """
    sorted_scores = sorted(round_floats(simulation_data), key=lambda d: (-d['scoring_function']), reverse=True)[:num_results]

    fasta_txt = "".join(f">{infile}|{now}|{score['replica_num']}|{score['sim_step']}|{score['scoring_function']}\n{score['sequence']}\n" for score in sorted_scores)

    with open(output_name + '_best_fasta.fas', 'w', encoding='utf-8') as fastafile:
        fastafile.write(fasta_txt)


def sort_and_filter_simulation_data(simulation_data, num_results, oligo):
    """
    Sorts and filters the simulation data based on specified criteria.

    Args:
    simulation_data (list): The simulation data to be processed.
    num_results (int): Number of top results to retain.
    oligo (str): The oligo option to determine the sorting criteria.

    Returns:
    list: Sorted and filtered list of simulation data.
    """
    remove_duplicates = {item['sequence']: item for item in simulation_data}
    data_no_duplicates = list(remove_duplicates.values())

    if oligo != "off":
        sorted_results = sorted(round_floats(data_no_duplicates), key=lambda d: (-d['oligo_fraction'], -d['mcc'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
    else:
        sorted_results = sorted(round_floats(data_no_duplicates), key=lambda d: (-d['mcc'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)

    return sorted_results[:num_results]


def generate_csv_from_data(sorted_data, output_name):
    """
    Generates a CSV file from sorted simulation data.

    Args:
    sorted_data (list): Sorted simulation data to be written to CSV.
    output_name (str): Name of the output CSV file.
    """
    fieldnames = sorted_data[0].keys()
    with open(output_name + '_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for data in sorted_data:
            writer.writerow(data)


def check_if_design_solved(sorted_results, input_file, oligo):
    """
    Checks the correctness of the results and generates the corresponding text output.

    Args:
    sorted_results (list): The sorted results to check.
    input_file (InputFile): The input file object containing sequence data.
    oligo (str): The oligo option to determine additional text.

    Returns:
    str: Text indicating the correctness of the results.
    """
    correct_count = sum(result['mcc'] == 0.0 for result in sorted_results[:10])
    correct_bool = correct_count > 0

    if oligo != "off":
        oligo_txt = f",oligo fraction: {sorted_results[0]['oligo_fraction']}" if correct_bool else ""
    else:
        oligo_txt = ""

    correct_result_txt = f">{input_file.name},{correct_bool},{correct_count},{sorted_results[0]['sequence']},{sorted_results[0]['mfe_ss']}{oligo_txt}"
    return correct_result_txt, correct_bool, correct_count


def write_best_str_file(correct_result_txt, output_name):
    """
    Writes the correctness text to a file.

    Args:
    correct_result_txt (str): The correctness text to be written to the file.
    output_name (str): The base name for the output file.
    """
    with open(output_name + '_best_str', 'w', newline='', encoding='utf-8') as myfile:
        myfile.write(correct_result_txt)


def generate_simulation_stats_text(stats, sorted_results, input_file, timlim, correct_bool, finish_time):
    """
    Generates text summarizing the statistics of the simulation.

    Args:
    stats (Stats): Object containing simulation statistics.
    sorted_results (list): The sorted results for correctness check.
    input_file (InputFile): Object containing input file data.
    timlim (int): Time limit for the simulation.
    correct_bool (bool): Indicates if the design was solved successfully.

    Returns:
    str: Text summarizing the simulation statistics.
    """
    sum_mc = stats.acc_mc_step + stats.rej_mc_step
    acc_perc = round(stats.acc_mc_step / sum_mc, 3) if sum_mc else 0
    sum_mc_metro = sum_mc - stats.acc_mc_better_e
    acc_metro = stats.acc_mc_step - stats.acc_mc_better_e
    sum_replica_att = stats.acc_re_step + stats.rej_re_step if stats.acc_re_step + stats.rej_re_step else 1
    formatted_time = time.strftime("%H:%M:%S", time.gmtime(finish_time))

    best_solution_txt = "\nDesign solved succesfully!\n\nBest solution:\n" if correct_bool else "\nDesign not solved!\n\nTarget structure:\n"
    best_solution_txt += f"{sorted_results[0]['sequence']}\nMFE Secondary Structure: {sorted_results[0]['mfe_ss']}\nPartition Function Energy: {round(sorted_results[0]['Epf'], 3)}\n1-MCC: {round(sorted_results[0]['mcc'], 3)}\n"

    stats_txt = f">{outname} time={timlim}s\nAcc_ratio={acc_perc}, Iterations={stats.step}, Accepted={stats.acc_mc_step}/{sum_mc}, Rejected={stats.rej_mc_step}/{sum_mc}\n" \
                f"Accepted Metropolis={acc_metro}/{sum_mc_metro}, Rejected Metropolis={sum_mc_metro - acc_metro}/{sum_mc_metro}\n" \
                f"Replica exchange attempts: {stats.global_step}\nReplica swaps attempts: {sum_replica_att}\nReplica swaps accepted: {stats.acc_re_step}\n" \
                f"Replica swaps rejected: {stats.rej_re_step}\nReplica exchange acc_ratio: {round(stats.acc_re_step / sum_replica_att, 3)}\n{best_solution_txt}\n\n" \
                f"Simulation time: {formatted_time}\n"

    return stats_txt


def write_stats_to_file(stats_txt, output_name):
    """
    Writes simulation statistics text to a file.

    Args:
    stats_txt (str): The simulation statistics text to be written.
    output_name (str): The base name for the output file.
    """
    with open(output_name + '_stats', 'w', newline='\n', encoding='utf-8') as myfile:
        myfile.write(stats_txt)


def move_results(file_list, directory, output_name):
    """
    Moves a list of files to the specified directory.

    Args:
    file_list (list): List of file names to be moved.
    directory (str): The target directory where files should be moved.
    output_name (str): The base name for the output files.
    """
    for file in file_list:
        move(output_name + file, directory + output_name + file)


def parse_and_output_results(simulation_data, input_file, stats, finish_time):
    """
    Parses simulation data and generates various outputs including CSV files, FASTA files, and statistics.

    Args:
    simulation_data (list): The simulation data to be processed.
    input_file (InputFile): The input file object containing sequence data.
    stats (Stats): Object containing simulation statistics.
    """
    # Generate necessary files and statistics
    sorted_simulation_trajectory = sort_trajectory(simulation_data)
    generate_trajectory_csv(sorted_simulation_trajectory, outname)
    generate_multifasta(sorted_simulation_trajectory, outname, infile, now)
    generate_replica_csv(simulation_data, outname)
    generate_best_fasta(simulation_data, num_results, outname, infile, now)

    sorted_results = sort_and_filter_simulation_data(simulation_data, num_results, oligo)
    generate_csv_from_data(sorted_results, outname)

    correct_result_txt, correct_bool, correct = check_if_design_solved(sorted_results[:10], input_file, oligo)
    write_best_str_file(correct_result_txt, outname)
    func.plot_simulation_data_combined(simulation_data, outname, input_file.alt_sec_struct, infile)

    stats_txt = generate_simulation_stats_text(stats, sorted_results, input_file, timlim, correct_bool, finish_time)
    print('\n' + stats_txt)
    write_stats_to_file(stats_txt, outname)

    # Move results to the specified directory
    files_to_move = ['_best_str', '_traj.csv', '_replicas.csv', '_best_fasta.fas', '_multifasta.fas', '_random.csv', '.command']
    move_results(files_to_move, "trajectory_files/", outname)


def run_functions():
    """
    This function runs the design process, which includes mutation, scoring, and selection.

    Parameters:
    sequence_o (str): The initial sequence.
    nt_list (list): A list of possible nucleotides to use in mutation.
    input_file (object): An object that encapsulates the input parameters.
    temps (list): A list of temperature values to use in the design process.
    steps (int): The number of steps to run in the design process.
    scoring_f (str): The scoring function to use.
    mutations (str): The type of mutations to use.
    exchange_rate (int): The rate at which to attempt exchanges between replicas.
    stats_o (object): A statistical object to keep track of the design process.

    Returns:
    list: A list of sequences that are the result of the design process.
    """

    nt_list = initialize_simulation(input_file)
    seqence_score_list = generate_sequences(nt_list, input_file)

    start_time = time.time()
    simulation_data = []

    for i in range(len(seqence_score_list)):
        simulation_data.append(vars(seqence_score_list[i]))

    handle_non_mutable_sequence(input_file, simulation_data, outname)

    stats = func.Stats()

    while time.time() - start_time < timlim:

        print('ETA', round((timlim - (time.time() - start_time)), 0), 'seconds', end='\r')
#        print(vars(stats))

        stats.update_global_step()
        seqence_score_list, stats = mutate_sequence_re(seqence_score_list, nt_list, stats)

        seqence_score_list, stats = replica_exchange(seqence_score_list, stats)

        for i in range(len(seqence_score_list)):
            seqence_score_list[i].get_sim_step(stats.step)
            simulation_data.append(vars(seqence_score_list[i]))

        if stats.global_step % 10 == 0:
            mcc_zero = process_intermediate_results(simulation_data, oligo, outname, num_results)
            if mcc_zero == num_results and sws == "on":
                break

        if stats.global_step == RE_steps:
            break
    finish_time = time.time() - start_time

    parse_and_output_results(simulation_data, input_file, stats, finish_time)


def get_outname(infile, replicas, timlim, acgu_percentages, pks, T_max, T_min, oligo, dimer, param, RE_attempt, scoring_f, point_mutations, alt_ss, tshelves, in_seed, subopt):
    """
    Generate an output name for a file based on the current parameters.

    The output name includes:
    - The name of the input file (excluding the extension)
    - The number of replicas
    - The energy
    - The time limit
    - The PKS setting
    - The ACGU percentages
    - The minimum and maximum temperatures
    - The parameter setting
    - The scoring function
    - The oligomerization setting
    - The dimer setting
    - The mutation setting
    - The point mutation setting

    Returns:
    -------
    outname : str
        The generated output name.
    """

    for i in range(0, len(scoring_f)):
        if scoring_f[i][0] not in available_scoring_functions:
            print(scoring_f[i][0], "is not an available option for scoring function. Check your command.")
            sys.exit()

    scoring_f_str = "_".join([f"{func}_{weight}" for func, weight in scoring_f])

    outname = infile.split(".")[0] + '_R' + str(replicas) + "_e" + str(RE_attempt) + "_t" + str(timlim) + "_pk" + str(pks) + "_ACGU" + str(acgu_percentages) +\
        "_Tmin" + str(T_min) + "_Tmax" + str(T_max) + "_p" + str(param) + "_SF" + scoring_f_str + "_O" + (str(oligo)) + "_D" + (str(dimer)) + "_PM" + (str(point_mutations))
    return outname


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


class InputFile:
    """
    This class represents an input file for the RNA sequence design problem.

    Attributes:
    -----------
    name : str
        The name of the input file.
    sec_struct : str
        The secondary structure in dot-bracket notation.
    seq_restr : str
        The sequence constraints.
    pairs : list
        The list of base pairs in the secondary structure.
    seed_seq : str
        The seed sequence for the design process.
    alt_sec_struct : str
        An alternative secondary structure for the design process.
    """

    def __init__(self, name, sec_struct, seq_restr):
        """
        Initializes the InputFile instance with the given name, secondary structure, and sequence constraints.

        Parameters:
        name (str): The name of the input file.
        sec_struct (str): The secondary structure in dot-bracket notation.
        seq_restr (str): The sequence constraints.
        """
        self.name = name
        self.sec_struct = sec_struct
        self.seq_restr = seq_restr
        self.pairs = []
        self.seed_seq = None
        self.alt_sec_struct = None
        self.alt_sec_structs = None
        self.target_pairs_tupl = {}

    def add_seed_seq(self, seed_seq):
        """
        Adds a seed sequence to the input file data.

        The seed sequence is used as a starting point for the RNA sequence design process.

        Parameters:
        seed_seq (str): The seed sequence to be added.
        """
        self.seed_seq = seed_seq

    def add_alt_sec_struct(self, alt_sec_structs):
        """
        Adds alternative secondary structures to the input file data.

        Alternative secondary structures are used in the RNA sequence design process to provide additional constraints or goals.

        Parameters:
        alt_sec_structs (list of str): A list of alternative secondary structures in dot-bracket notation.
        """
        self.alt_sec_struct = alt_sec_structs[0].strip()
        self.alt_sec_structs = [x.strip() for x in alt_sec_structs]

    def set_target_pairs_tupl(self):
        self.target_pairs_tupl = {tuple(pair) for pair in self.pairs}


def redirect_stderr_to_file(filename):
    """
    Redirects the standard error (stderr) stream to a file.

    Parameters:
    filename (str): The name of the file to which stderr will be redirected.

    Returns:
    tuple: A tuple containing the original stderr file descriptor and the opened log file object.
    """

    # Backup original stderr
    stderr_fileno = sys.stderr.fileno()
    stderr_save = os.dup(stderr_fileno)
    stderr_log = open(filename, 'w', encoding='utf-8')

    # Redirect stderr to our log file
    os.dup2(stderr_log.fileno(), stderr_fileno)

    return stderr_save, stderr_log


def restore_stderr(original_stderr):
    """
    Restores the standard error (stderr) stream from the original backed up version.

    Parameters:
    original_stderr (int): The file descriptor of the original stderr before redirection.
    """

    os.dup2(original_stderr, sys.stderr.fileno())


def print_filtered_log(filename, exclude_text):
    """
    Prints the content of a log file, excluding lines that contain specific text.

    Parameters:
    filename (str): The name of the log file to read.
    exclude_text (str): The text to exclude from printing.
    """

    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if exclude_text not in line:
                print(line, end='', file=sys.stderr)


if __name__ == "__main__":

    log_filename = "tmp_log.txt"
    original_stderr, stderr_log = redirect_stderr_to_file(log_filename)

    try:

        if len(sys.argv) == 1:
            print("DesiRNA")
            print("usage: DesiRNA.py [-h] [-H]")
            print("\nOptions:")
            print("  -h  Display basic usage and options.")
            print("  -H  Display detailed help, including advanced options.")
            sys.exit(1)

        print("DesiRNA")

        command = os.path.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:])

        script_path = os.path.dirname(os.path.abspath(__file__))

        now = datetime.now()
        now = now.strftime("%Y%m%d.%H%M%S")

        available_scoring_functions = ['Ed-Epf', '1-MCC', 'sln_Epf', 'Ed-MFE', '1-precision', '1-recall', 'Edef']

        infile, replicas, timlim, acgu_percentages, T_max, T_min, oligo, param, RE_attempt, scoring_f, point_mutations,  \
            tshelves, in_seed, subopt, diff_start_replicas, num_results, acgu_content, RE_steps, tm_max, tm_min, motifs, dimer, sws = argument_parser()

        input_file = read_input(infile)
        print(infile)

        if RE_steps != None:
            timlim = 100000000000000000

        if param == '1999':
            RNA.params_load(os.path.join(script_path, "rna_turner1999.par"))

        L = 504.12

        if acgu_content == '':
            nt_percentages = {"A": 15, "C": 30, "G": 30, "U": 15}
        else:
            acgu_l = [int(x) for x in acgu_content.split(',')]
            if sum(acgu_l) != 100:
                print("The ACGU content should sum up to 100, check your command.")
                sys.exit()
            nt_percentages = {"A": acgu_l[0], "C": acgu_l[1], "G": acgu_l[2], "U": acgu_l[3]}

        if input_file.alt_sec_struct != None:
            alt_ss = "on"
        else:
            alt_ss = "off"

        oligo_state = "none"

        if "&" in input_file.sec_struct and dimer == "off":
            too_much = input_file.sec_struct.count("&")
            if too_much != 1:
                print("\nToo much structures in the input. Can only design RNA complexes of max two sequences.\nPlease correct your input file.\n")
                sys.exit()
            oligo_state = "heterodimer"
        elif "&" in input_file.sec_struct and dimer == "on":
            oligo_state = "homodimer"
        elif "&" not in input_file.sec_struct and oligo == "on":
            oligo_state = "avoid"

        if set(input_file.sec_struct).issubset('.()&'):
            pks = "off"
        else:
            pks = "on"

        filename = os.path.basename(infile)

        outname = get_outname(filename, replicas, timlim, acgu_percentages, pks, T_max, T_min, oligo, dimer,
                              param, RE_attempt, scoring_f, point_mutations, alt_ss, tshelves, in_seed, subopt)

        if tshelves == '':
            rep_temps_shelfs = get_rep_temps()
        else:
            rep_temps_shelfs = [float(temp) for temp in tshelves.split(",")]

        if in_seed != 0:
            original_seed = 2137 + in_seed
        else:
            original_seed = random.random()

        random.seed(original_seed)
        outname = os.path.basename(outname)

        WORK_DIR = outname + "_" + str(now)
        Path(WORK_DIR).mkdir(parents=True, exist_ok=True)
        copy(infile, WORK_DIR + "/" + filename)
        Path(WORK_DIR + "/trajectory_files").mkdir(parents=True, exist_ok=True)
        os.chdir(WORK_DIR)

        with open(outname + ".command", 'w', encoding='utf-8') as f:
            print(command, file=f)

        run_functions()

    except Exception as e:
        # Close the log, restore stderr, and raise the exception
        stderr_log.close()
        restore_stderr(original_stderr)
        print("Exception occurred:", e)

    else:
        # Restore stderr
        restore_stderr(original_stderr)
        stderr_log.close()

        # Now, read the log and print messages that are not the warning
        print_filtered_log("../" + log_filename, "WARNING: pf_scale too large")

    os.remove("../" + log_filename)
