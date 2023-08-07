#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse as argparse
import sys
import random
import RNA
import math
import csv
import time
import random
import re
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime
import multiprocessing as mp
import pandas as pd
from pathlib import Path
from shutil import copy
from shutil import move

from utils import functions_classes as func
from utils.SimScore import SimScore
from utils.pcofold_dimer_multichain_energy import mfe_e_dimer

iupac_re = {'A' : 'A',
            'C' : 'C',
            'G' : 'G',
            'T' : 'T',
            'U' : 'U',
            'W' : '[AU]',
            'S' : '[GC]',
            'M' : '[AC]',
            'K' : '[GU]',
            'R' : '[AG]',
            'Y' : '[CU]',
            'B' : '[CGU]',
            'D' : '[AGU]',
            'H' : '[ACU]',
            'V' : '[ACG]',
            'N' : '[ACGU]',
            }

def argument_parser():
    """
    Parses command-line arguments.

    Returns:
        A Namespace containing the parsed command-line arguments.
    """
    
    parser = argparse.ArgumentParser(description=__doc__, prog='DesiRNA.py')
    parser.add_argument("-f", "--filename", required=True, dest="name",
                        help="Name of a file that contains secondary structures and constraints.")
    parser.add_argument("-R", "--replicas", required=False, dest="replicas", default=10, type=int,
                            help="Number of replicas. [default = 10]")
    parser.add_argument("-e", "--exchange", required=False, dest="exchange", default=100, type=int,
                            help="Frequency of replixa exchange attempt. [default = 100]")
    parser.add_argument("-t", "--timelimit", required=False, default=60, dest="timlim", type=int,
                            help="Timelimit for running the program [s]. [default = 60]")
    parser.add_argument("-s", "--steps", required=False, default=None, dest="steps", type=int,
                            help="Number of Replica Exchange steps after which the simulation ends. Overwrites the -t option. [default = None]")


    parser.add_argument("-r", "--results_number", required=False, dest="num_results", default=10, type=int,
                            help="Number of best results to be reported in the output. [default = 10]")

    parser.add_argument("-p", "--param", required=False, dest="param", default='1999', choices=['2004','1999'],
                            help="Turner energy parameter for calculating MFE. [default = 1999]")

    parser.add_argument("-tmin", "--tmin", required=False, dest="t_min", default=10, type=float,
                            help="Minimal Replica Temperature. [default = 10]")
    parser.add_argument("-tmax", "--tmax", required=False, dest="t_max", default=150, type=float,
                            help="Maximal Replica Temperature. [default = 150]")
    parser.add_argument("-ts", "--tshelves", required=False, dest="tshelves", type=str, default='',
                            help="Custom temperature shelves for replicas in replica exchange simulation. Provide comma-separated values.")

    parser.add_argument("--motifs", required=False, dest="motifs", type=str, default='',
                        help="Sequence motifs along with their bonuses(-)/penalties(+). Provide comma-separated key,value,key,value sequence.")


    parser.add_argument("-sf", "--scoring_function", required=False, dest="scoring_f", default='ed-mfe', choices=['ed-mfe','1-mcc', 'sln_mfe'], nargs='+',
                            help="Scoring function used to guide the design process. [default = ed-mfe]")




    parser.add_argument("-nd", "--negative_design", required=False, dest="subopt", default='off', choices=['off','on'],
                            help="Use negative design approach. [default = off]")


                            
    parser.add_argument("-acgu", "--ACGU", required=False, dest="percs", default='off', choices=['off','on'],
                            help="Keep 'natural' ACGU content. If turned on the content will be A:15%%, C:30%%, G:30%%, U:15%%. [default = off]")
    parser.add_argument("-acgu_content", "--ACGU_content", required=False, dest="acgu_content", default='', type=str,
                            help="Provide user defined ACGU content. Comma-separated values e.g., 10,40,40,10")
    
#    parser.add_argument("-pk", "--PK", required=False, dest="pks", default='off', choices=['off','on'],
#                            help="Design of pseudoknotted structures. [default = off, turns on automatically if pseudoknot is detected in the input file]")
#    parser.add_argument("-tr", "--treplicaexchange", required=False, dest="t_re", default=10, type=float,
#                            help="Temperature of replica exchange attempt. [default = 10]")

    parser.add_argument("-o", "--oligomerization", required=False, dest="oligo", default='off', choices=['off','on'],
                            help="Check if the designed sequence tends to oligomerize. Slows down the simulation. [default = off]")
#    parser.add_argument("-d", "--dimer", required=False, dest="dimer", default='off', choices=['off','on'],
#                            help="Design of a RNA complex, of two strands. [deafult = off, turns on automatically is '&' detected in the input file]")

#    parser.add_argument("-m", "--mutations", required=False, dest="mutations", default='one', choices=['one','multi'],
#                            help="More mutatuions pre one MC step in higher temperature replicas. Slows down the simulation. [default = off]")    
    parser.add_argument("-tm", "--target_mutations", required=False, dest="pm", default='on', choices=['off','on'],
                            help="Targeted mutations. Targets mostly False Negativeas and False Positives. [default = on]")

    parser.add_argument("-tm_perc_max", "--target_mutations_percentage_max", required=False, dest="tm_max", default=0.7, type=float,
                            help="Highest percentage of targeted mutations applied to lowest temperature replica. Percentage for replicas in between will be set evenly from 'tm_perc_max' to 'tm_perc_min'. Float from 0.0 to 1.0. [default = 0.7]")
    parser.add_argument("-tm_perc_min", "--target_mutations_percentage_min", required=False, dest="tm_min", default=0.0, type=float,
                            help="Lowest percentage of targeted mutations applied to highest temperature replica. Percentage for replicas in between will be set evenly from 'tm_perc_max' to 'tm_perc_min'. Float from 0.0 to 1.0. [default = 0.0]")


#    parser.add_argument("-a", "--alt_ss", required=False, dest="alt_ss", default='off', choices=['off','on'],
#                            help="Design of sequences folding into two structures. [default = off, turns on automatically if alternative structure is detected in the input file]")
    parser.add_argument("-seed", "--seed", required=False, default=0, dest="in_seed", type=int,
                            help="User defined seed number for simulation. [default = 0]")
                            
    parser.add_argument("-re_seq", "--replicas_sequences", required=False, dest="diff_start_replicas", default='one', choices=['different','same'],
                            help="Choose wether replicas will start from the same random sequence or each replica will start from different random sequence. [default = same]")
#    parser.add_argument("-rand_seed", "--random_seed", required=False, dest="rand_seed", default='on', choices=['off','on'],
#                            help="Use random seed number for simulation. [default = on]")


                            
    args = parser.parse_args() 
    
    if args.tshelves:
        temperature_shelves = [float(temp) for temp in args.tshelves.split(",")]

        # Check if the number of temperatures matches the number of replicas
        if len(temperature_shelves) != args.replicas:
            parser.error("The number of temperatures provided in -ts/--tshelves must match the number of replicas set by -R.")

    if args.motifs:

        motifs_list = args.motifs.split(',')
        try:
            motifs = {motifs_list[i]:(re.compile(''.join([iupac_re[ch]
                                                          for ch in motifs_list[i]])),
                                      float(motifs_list[i+1])) for i in range(0,len(motifs_list),2)}
        except:
            parser.error("Something is wrong with your motifs")
    else:
        motifs = {}
    
    if args.percs == 'off' and args.acgu_content:
        parser.error("The -acgu_content option is only applicable when -acgu is set to 'on'.")    
    
    if args.pm == 'off' and ('tm_max' in vars(args) or 'tm_min' in vars(args)):
        parser.error("Targeted mutations must be 'on' to set -tm_perc_max and/or -tm_perc_min.")

    infile = args.name
    replicas = args.replicas
    timlim = args.timlim
    acgu_percs = args.percs
#    pks = args.pks
    t_max = args.t_max
    t_min = args.t_min
    oligo = args.oligo
#    dimer = args.dimer
    param = args.param
    exchange_rate = args.exchange
    scoring_f = args.scoring_f
#    mutations = args.mutations
    pm = args.pm
    steps = args.steps
#    alt_ss = args.alt_ss
    tshelves = args.tshelves
#    t_re = args.t_re
    in_seed = args.in_seed
#    rand_seed = args.rand_seed

    subopt = args.subopt
    diff_start_replicas = args.diff_start_replicas
    num_results = args.num_results
    acgu_content = args.acgu_content
    tm_max = args.tm_max
    tm_min = args.tm_min

    
    return infile, replicas, timlim, acgu_percs,  t_max, t_min, oligo,  param, exchange_rate, scoring_f, pm,   tshelves,\
           in_seed, subopt, diff_start_replicas, num_results, acgu_content, steps, tm_max, tm_min, motifs


def read_input(infile):
    """
    Reads the input file and returns an InputFile object.

    Returns:
        input (InputFile): An InputFile object containing the data from the input file.
    """
    
    string = ""
    with open(infile) as f:
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
        input_file.add_alt_sec_struct(data_dict['alt_sec_struct'][0])
    
    return input_file


def check_dot_bracket(ss):
    """
    Checks if the given secondary structure string is in the correct dot-bracket notation.
    It ensures that every opening bracket has a matching closing bracket, and vice versa.

    Args:
        ss (str): The secondary structure string in dot-bracket notation.

    Returns:
        pairs_list (list): A list of pairs, where each pair is a list of two indices representing a base pair.

    Raises:
        SystemExit: If the secondary structure string is not in the correct dot-bracket notation.
    """
    
    
    db_list = [['(',')'],['[',']'],['<','>'],['{','}'],['A','a'],['B','b'],['C','c'],['D','d'],['E','e']]
    allowed_characters = '()[]<>{}AaBbCcDdEe.&'
    
    for i in range(0, len(ss)):
        if ss[i] not in allowed_characters:
            sys.exit("Not allowed characters in structures. Check input file.")    
    
    stack_list =[]
    pairs_list =[]

    # stack-pop for all versions of brackets form the db_list    
    for i in range(0, len(db_list)):
        for c, s in enumerate(ss):
            if s == db_list[i][0]:
                stack_list.append(c)
            elif s == db_list[i][1]:
                if len(stack_list) == 0:
                    sys.exit("There is no opening bracket for nt position "+str(c+1)+'-'+ss[c])
                elif s == db_list[i][1]:
                    pairs_list.append([stack_list.pop(), c])
        if len(stack_list) > 0:
            err = stack_list.pop()
            sys.exit("There is no closing bracket for nt position "+str(err)+'-'+ss[err])

    pairs_list_clean = [x for x in pairs_list if x != []]

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
            sys.exit("Not allowed characters in sequence restraints. Check input file.")


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
        sys.exit("Secondary structure and sequence restraints are of different length. Check input file.")


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
        obj = Nucleotide(number = i)
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

    Args:
        nt (str): The nucleotide character.

    Returns:
        nt_all (list): A list of nucleotide characters.
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

    Args:
        nt_list (list): A list of Nucleotide objects.

    Raises:
        SystemExit: If the logic of the input is incorrect.
    """
    
    for i in range(0, len(nt_list)):
        if nt_list[i].pairs_with != None:
            if_pair(nt_list[i], nt_list[nt_list[i].pairs_with])

    
def if_pair(nt1, nt2):
    """
    Checks if two nucleotide characters can form a base pair.

    Args:
        nt1 (str): The first nucleotide character.
        nt2 (str): The second nucleotide character.

    Raises:
        SystemExit: If the two nucleotide characters cannot form a base pair.
    """
    
    pair_dict = {'A': ['U'], 'U': ['G','A'], 'G': ['U','C'], 'C': ['G'] }
    
    store = False
    for i in range(0, len(nt1.letters)):
        nt1l = nt1.letters[i]
        for j in range(0, len(nt2.letters)):
            nt2l = nt2.letters[j]
            if nt2l in pair_dict[nt1l]:
                store =  True

    if store == True:
        pass
    else:
        sys.exit("Wrong restraints in the input file. Nucleotide "\
        + str(nt1.number+1)+" "+str(nt1.letters)+", cannot pair with nucleotide "\
        + str(nt2.number+1)+ " "+str(nt2.letters) )


def wc_pair(nt1):
    """
    Returns the Watson-Crick pair of a given nucleotide character.

    Args:
        nt1 (str): The nucleotide character.

    Returns:
        nt2 (str): The Watson-Crick pair of the nucleotide character.
    """

    pair_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G' }
    
    nt2 = pair_dict[nt1]

    return nt2


def can_pair(nt):
    """
    Returns a list of nucleotide characters that can form a base pair with a given nucleotide character.

    Args:
        nt (str): The nucleotide character.

    Returns:
        pairing_l (list): A list of nucleotide characters.
    """
    
    pair_dict = {'A': ['U'], 'U': ['G','A'], 'G': ['U','C'], 'C': ['G'] }
    
    pairing_l = pair_dict[nt]

    return pairing_l


def random_sequence_generator(nt_list,input_file):
    """
    Generates a random sequence that satisfies the sequence restraints.

    Args:
        nt_list (list): A list of Nucleotide objects.
        input_file (InputFile): The InputFile object containing the data from the input file.

    Returns:
        result_sequence (str): The randomly generated sequence.
    """
    
    seq_l = list(input_file.seq_restr).copy()
    pair_list = sorted(input_file.pairs.copy())

    
    for i in range(0 , len(nt_list)):
        if nt_list[i].pairs_with == None:
            random.choice(nt_dictionary(seq_l[i]))
    
    for i in range(0, len(pair_list)):
        nt1 = nt_list[pair_list[i][0]].letters_allowed
        nt2 = nt_list[pair_list[i][1]].letters_allowed
        nt1.sort()
        nt2.sort()
        nt1num =nt_list[pair_list[i][0]].number 
        nt2num =nt_list[pair_list[i][1]].number

        allowed_choices = allowed_choice(nt1, nt_percentages)

        seq_l[nt1num] = random.choices(nt1, weights = allowed_choices)[0]
        seq_l[nt2num] = wc_pair(seq_l[nt1num])

    for i in range(0, len(seq_l)):
        if seq_l[i] not in ["A","C","G","U"]:
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
    for i in range(0 , len(nt_list)):
        if nt_list[i].pairs_with == None and "A" in nt_list[i].letters:
            seq_l[i] = "A"
    
    # Loop boosting - G (or if not possible = U) as the first of nonbonded nts, unless its (.( - single bulge, then stay A
    for i in range(1, len(nt_list)-1):
        if (nt_list[i].pairs_with == None and nt_list[i-1].pairs_with != None) and (nt_list[i].pairs_with == None and nt_list[i+1].pairs_with == None):
            if "G" in nt_list[i].letters:
                seq_l[i] = "G"
            elif "U" in nt_list[i].letters:
                seq_l[i] = "U"
    

    if acgu_percentages == 'off':    

    # paired nucleotides handling, if possible GC, if not AU
        for i in range(0, len(pair_list)):
            nt1 = nt_list[pair_list[i][0]].letters_allowed
            nt2 = nt_list[pair_list[i][1]].letters_allowed
            nt1num =nt_list[pair_list[i][0]].number 
            nt2num =nt_list[pair_list[i][1]].number
            if ("C" in nt1 and "G" in nt1) and ("C" in nt2 and "G" in nt2):
                seq_l[nt1num] = random.choice(["C","G"])
                seq_l[nt2num] = wc_pair(seq_l[nt1num])
            elif "C" in nt1 and "G" in nt2:
                seq_l[nt1num] = "C" 
                seq_l[nt2num] = "G"
            elif "G" in nt1 and "C" in nt2:
                seq_l[nt1num] = "G" 
                seq_l[nt2num] = "C"
            elif ("A" in nt1 and "U" in nt1) and ("A" in nt2 and "U" in nt2):
                seq_l[nt1num] = random.choice(["A","U"])
                seq_l[nt2num] = wc_pair(seq_l[nt1num])
            elif "A" in nt1 and "U" in nt2:
                seq_l[nt1num] = "A" 
                seq_l[nt2num] = "U"
            elif "U" in nt1 and "A" in nt2:
                seq_l[nt1num] = "U" 
                seq_l[nt2num] = "A"
    
    
    # weightened choice here!
    
    elif acgu_percentages == 'on':
        for i in range(0, len(pair_list)):
            nt1 = nt_list[pair_list[i][0]].letters_allowed
            nt2 = nt_list[pair_list[i][1]].letters_allowed
            nt1.sort()
            nt2.sort()
            nt1num =nt_list[pair_list[i][0]].number 
            nt2num =nt_list[pair_list[i][1]].number
        
            allowed_choices = allowed_choice(nt1, nt_percentages)

            seq_l[nt1num] = random.choices(nt1, weights = allowed_choices)[0]

            seq_l[nt2num] = wc_pair(seq_l[nt1num])
        


    # random choice (from allowed letters) of whatever left    
    for i in range(0, len(seq_l)):
        if seq_l[i] not in ["A","C","G","U"]:
            seq_l[i] = random.choice(nt_dictionary(seq_l[i]))	
            


    result_sequence = ''.join(seq_l)

    if input_file.seed_seq:
        result_sequence = input_file.seed_seq

    return result_sequence


def allowed_choice(allowed, percs):
    """
    Determines the allowed mutations based on nucleotide percentages.

    Args:
        available_mutations (list): List of available mutations.
        nt_percentages (dict): Dictionary with the percentage of each nucleotide.

    Returns:
        list: List of weights for the available mutations.
    """
    
    percs_allowed = percs.copy()
    
    #allowed=["C","U"]
    for nt in percs:
        if nt not in allowed:
            percs_allowed[nt] = 0
    
    
    list_percs_allowed = list(percs_allowed.values())
 
    list_percs_allowed = [i for i in list_percs_allowed if i != 0]
    
    return list_percs_allowed


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
        delta = (T_max-T_min)/(replicas-1)
    else:
        delta = (T_max-T_min)/2
    #rep_temps = np.arange(T_min, T_max, delta)
    rep_temps = []

    T_curr = T_min
 
    for i in range(replicas):
        if replicas != 1:
            rep_temps.append(round(T_curr,3))
            T_curr +=delta
        else:
            T_curr +=delta
            rep_temps.append(round(T_curr,3))

    if replicas ==1:
        rep_temps = [T_max]

    pot= 1.5
    rep_temps_mod = []
    
    for i in range(0, len(rep_temps)): 
        if i == 0:
            rep_temps_mod.append(rep_temps[i])
        else:
            rep_temps_mod.append((1-(1-(rep_temps[i]/T_max)**pot)**(1/pot))*T_max)

    una =False
    if una == True:
        rep_temps = rep_temps_mod    
    
    return rep_temps


def generate_initial_list(nt_list, input_file):
    
    sequence = initial_sequence_generator(nt_list, input_file)
    
    seq_list = []
    
    for i in range(0,replicas):
        if diff_start_replicas == "different":
            sequence = initial_sequence_generator(nt_list, input_file)
        sequence_object = score_sequence(sequence)
        sequence_object.get_replica_num(i+1)
        sequence_object.get_temp_shelf(rep_temps_shelfs[i])
        sequence_object.get_sim_step(0)
        seq_list.append(sequence_object)        
    
    
    init_sequences_txt =""
    for i in range(0, len(seq_list)):
        init_sequences_txt += str(vars(seq_list[i]))+"\n"
    
    return seq_list


def generate_initial_list_random(nt_list, input_file):
    """
    Generate an initial list of random RNA sequences.

    Args:
        nt_list (list): List of nucleotide objects.
        input_file (object): An InputFile object that contains the secondary structure of the RNA.

    Returns:
        list: List of ScoreSeq objects each representing a random RNA sequence.
    """
    
    seq_list = []
    for i in range(0,replicas):
        sequence = random_sequence_generator(nt_list, input_file)
        deltaF = delta_MFE_EOS(input_file.sec_struct, sequence)
        score = metropolis_score(310, deltaF)
        seq_list.append([sequence, score])

    seq_list.sort(key = lambda x: x[1], reverse = True)


    result_list = seq_list[:replicas]

    with open(outname+'_random.csv', 'w', newline='\n') as myfile:
         wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
         wr.writerows(result_list)

    return result_list


def get_seq_to_mutate(sequence_list, step):
    """
    Select a sequence to mutate from the list of nucleotides.

    Args:
        nt_list (list): List of nucleotide objects.

    Returns:
        list: List of selected nucleotides to mutate.
    """

    sequence = sequence_list[step]
    
    return sequence


def mutate_sequence_more(sequence_obj, nt_list):
    """
    Mutate a list of sequences in parallel.

    Args:
        lst_seq_obj (list): List of ScoreSeq objects each representing an RNA sequence.
        nt_list (list): List of nucleotide objects.
        stats_obj (object): A Stats object containing statistics of the simulation.

    Returns:
        tuple: A tuple containing the list of mutated ScoreSeq objects and the updated Stats object.
    """
    
    sequence_mutated_re = sequence_obj

    part_size = len(rep_temps_shelfs) // 3

    rep_temps_shelfs_sliced = [rep_temps_shelfs[i*part_size:i*part_size + part_size if i != 2 else None] for i in range(3)]

    if sequence_obj.temp_shelf in rep_temps_shelfs_sliced[0]:
        
        sequence_mutated_re = mutate_sequence(sequence_mutated_re, nt_list)
    elif sequence_obj.temp_shelf in rep_temps_shelfs_sliced[1]:
        for i in range(0,2):
            sequence_mutated_re = mutate_sequence(sequence_mutated_re, nt_list)
    else:
        for i in range(0,3):
            sequence_mutated_re = mutate_sequence(sequence_mutated_re, nt_list)

    return sequence_mutated_re


def get_mutation_position(seq_obj, available_positions):
    """
    Determine the position of a mutation in the sequence.

    Args:
        sequence_obj (object): A ScoreSeq object representing an RNA sequence.
        nt_list (list): List of nucleotide objects.

    Returns:
        int: Position of the mutation in the sequence.
    """

    if point_mutations == "off":
        mutation_position = random.choice(available_positions)
    elif point_mutations == "on":

        max_perc_prob = tm_max
        min_perc_prob = tm_min
        
        pair_list_mfe = check_dot_bracket(seq_obj.mfe_ss)
        
        query_structure = {tuple(pair) for pair in pair_list_mfe}

        false_negatives = target_pairs_tupl - query_structure
        false_negatives = [item for tup in false_negatives for item in tup]
        false_negatives = [item for item in false_negatives if item in available_positions]

        false_positives = query_structure - target_pairs_tupl
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
            false_cases = expand_cases(false_cases, len(seq_obj.sequence)-1)
            range_pos = random.choices([false_cases, available_positions], weights=[mutat_point_prob, 1-mutat_point_prob])[0]

        mutation_position = random.choice(range_pos)
        
    return mutation_position


def expand_cases(cases, max_value, range_expansion=3):
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
        sequence_obj (object): A ScoreSeq object representing an RNA sequence.
        nt_list (list): List of nucleotide objects.

    Returns:
        object: A ScoreSeq object representing the mutated RNA sequence.
    """
    
    sequence = sequence_obj.sequence
    sequence_list = list(sequence)
    range_pos =[]
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
    elif nt_list[nt_pos].pairs_with != None:
        nt1 = nt_list[nt_pos]
        available_mutations_nt1 = nt1.letters_allowed.copy()
        if sequence_list[nt_pos] in available_mutations_nt1:
        
            available_mutations_nt1.remove(sequence_list[nt_pos])
        available_mutations_nt1.sort()
        
        nt2 = nt_list[nt1.pairs_with]
        allowed_mutations_nt2 = nt2.letters_allowed.copy()    
        nt2_pos = nt2.number
        if acgu_percentages == "on":        
            allowed_choices_nt1 = allowed_choice(available_mutations_nt1, nt_percentages)
            mutated_nt1 = random.choices(available_mutations_nt1, weights= allowed_choices_nt1)[0]
        elif acgu_percentages == "off": 
            mutated_nt1 = random.choice(available_mutations_nt1)        
        
        allowed_pairings_nt2 = can_pair(mutated_nt1)
        available_mutations_nt2 = list(set(allowed_mutations_nt2).intersection(allowed_pairings_nt2))

        allowed_mutations_nt2.sort()
        
        if acgu_percentages == "on":
            allowed_choices_nt2 = allowed_choice(available_mutations_nt2, nt_percentages)
            mutated_nt2 = random.choices(available_mutations_nt2, weights= allowed_choices_nt2)[0]
        elif acgu_percentages == "off":
            mutated_nt2 = random.choice(available_mutations_nt2)

        sequence_list[nt_pos] = mutated_nt1
        sequence_list[nt2_pos] = mutated_nt2	

        
    sequence_mutated = ''.join(sequence_list)
    mut_position = list(' '*len(sequence))
    mut_position[nt_pos] = "#"

    if nt2_pos != None:
        mut_position[nt2_pos] = "#"

    sequence_mutated = score_sequence(sequence_mutated)
    sequence_mutated.get_replica_num(sequence_obj.replica_num)
    sequence_mutated.get_temp_shelf(sequence_obj.temp_shelf)

    return sequence_mutated


def get_pk_struct(seq, ss_nopk):
    """
    Get the pseudoknot structure of the sequence.

    Args:
        seq (str): RNA sequence.
        ss_nopk (str): Secondary structure of the RNA sequence without pseudoknots.

    Returns:
        str: Pseudoknot structure of the RNA sequence.
    """

    constraints = ss_nopk.replace("(","x").replace(")","x")

    fc = RNA.fold_compound(seq)
    
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
        constraints = ss_pk.replace("(","x").replace(")","x").replace("[","x").replace("]","x")        
        fc = RNA.fold_compound(seq)
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
            constraints = ss_pk.replace("(","x").replace(")","x").replace("[","x").replace("]","x").replace("<","x").replace(">","x")
            fc = RNA.fold_compound(seq)
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


def delta_MFE_EOS(struct, seq):
    """
    Calculate the difference in Minimum Free Energy (MFE) and Ensemble Of Structure (EOF) of the given structure and sequence.

    Args:
        struct (str): Secondary structure of the RNA sequence.
        seq (str): RNA sequence.

    Returns:
        float: The difference in MFE and EOF.
    """
    
    F0 = RNA.energy_of_struct(seq, struct)
    F1 = RNA.pf_fold(seq)[1]
    
        
    delta_F = F0 - F1
    
    return delta_F


def mc_delta(deltaF_o, deltaF_m,  T_replica):
    """
    Determine whether to accept the mutation based on the Metropolis criterion.

    Args:
        deltaF_o (float): Original scoring function value.
        deltaF_m (float): Mutated scoring function value.
        T_replica (float): Temperature of the replica.

    Returns:
        tuple: A tuple containing two boolean values indicating whether the mutation is accepted and whether the mutated scoring function is better.
    """

    accept_e = False
    if deltaF_m <= deltaF_o:
        accept = True
        accept_e = True
    else:
        diff= deltaF_m - deltaF_o
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
    
    
    score = math.exp((-L/temp) * (dE))
    
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
        p = math.exp(L*(1/T0 -1/T1)*(dE0-dE1))
#        p = math.exp(L*(1/T_re)*(dE0-dE1)) # bardziej rownomierne wymiany replik
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

    sec_struct = input_file.sec_struct
    re_pairs = []

    num_shelfs = [i for i in range(1, len(replicas))]

    temps = []
    
    for i in range(0, len(replicas)):
        temps.append(replicas[i].temp_shelf)
    
    temps = sorted(temps)
    
    if stats_obj.global_step % 2 ==0:
    
        for i in range(0, len(num_shelfs)-1, 2):
            re_pairs.append([i+1, i+2])
    else:
        for i in range(0, len(num_shelfs), 2):
            re_pairs.append([i, i+1])
    
    replicas = sorted(replicas, key=lambda obj: obj.temp_shelf)


    for i in range(len(re_pairs)):
        rep_i = re_pairs[i][0] #replica temporary number
        rep_j = re_pairs[i][1] #replica temporary number
    
        T_i = replicas[rep_i].temp_shelf
        T_j = replicas[rep_j].temp_shelf
        seq_i =replicas[rep_i].sequence
        seq_j =replicas[rep_j].sequence
        
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


def get_mfe_e_ss(seq):
    """
    This function calculates the minimum free energy (MFE) and ensemble free energy (EFE) 
    of a given RNA sequence.

    Parameters:
    seq (str): The RNA sequence to analyze.

    Returns:
    (float, float): A tuple containing the MFE and EFE of the sequence.
    """

    if dimer == "off":
        pf_struct = RNA.pf_fold(seq)
        structure_nopk = pf_struct[0]
        energy = pf_struct[1]
        
        if (any(char in ",{}|" for char in structure_nopk)) or (structure_nopk.count("(") != structure_nopk.count(")")):
            structure_nopk = RNA.fold(seq)[0]
        
        if pks == "off":
            structure = structure_nopk
    
        elif pks == "on":
            structure = get_pk_struct(seq, structure_nopk)            
    
    elif dimer == "on":
        dimer_struct = mfe_e_dimer(seq)
        structure = dimer_struct[3]
        energy = dimer_struct[0]

#    if oligo == "on":
#        pf_struct = RNA.fold(seq)
#        structure_nopk = pf_struct[0]
#        energy = pf_struct[1]

    return energy, structure


def score_motifs(seq, motifs):
    """Returns the overall bonus/penalty for sequence motifs
       present/absent in the sequence"""
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
    
    scored_sequence = func.ScoreSeq(sequence = seq)


    mfe_energy, mfe_structure = get_mfe_e_ss(seq)
    
    scored_sequence.get_mfe(mfe_energy)
    scored_sequence.get_mfe_ss(mfe_structure)


    scored_sequence.get_edesired(RNA.energy_of_struct(seq, input_file.sec_struct.replace("&","")))

    scored_sequence.get_edesired_minus_mfe(scored_sequence.mfe, scored_sequence.edesired)


    ssc = SimScore(input_file.sec_struct.replace("&","Ee"), scored_sequence.mfe_ss.replace("&","Ee"))
    ssc.find_basepairs()
    ssc.cofusion_matrix()

    scored_sequence.get_precision(ssc.precision())
    scored_sequence.get_recall(ssc.recall())
    scored_sequence.get_mcc(ssc.mcc())


    scored_sequence.get_scoring_function(scoring_f)

    if 'sln_mfe' in scoring_f:
        scored_sequence.get_sln_mfe()
    
    if input_file.alt_sec_struct != None:
        scored_sequence.get_edesired2(RNA.energy_of_struct(seq, input_file.alt_sec_struct))
        scored_sequence.get_edesired2_minus_mfe(scored_sequence.mfe, scored_sequence.edesired2)    
        scored_sequence.get_scoring_function_w_alt_ss()

    if (subopt == "on") and (scored_sequence.mcc == 0):        
        scored_sequence.get_subopt_e(func.get_first_suboptimal_structure_and_energy(seq)[1])
        scored_sequence.get_esubopt_minus_mfe(scored_sequence.mfe, scored_sequence.subopt_e)
        scored_sequence.get_scoring_function_w_subopt()
        

    if oligo == "on":
        scored_sequence.get_oligomerization()
        scored_sequence.get_scoring_function_w_oligo()

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
    
#    best_score = ''
#    in_sequence = sequence_o
    worker_stats.reset_mc_stats()
    
    for i in range(0, RE_attempt):
        sequence_m = mutate_sequence(sequence_o, nt_list)
#        if mutations == 'multi':
#            sequence_m = mutate_sequence(sequence_m, nt_list)
#            sequence_m = mutate_sequence(sequence_m, nt_list)
#            sequence_m = mutate_sequence(sequence_m, nt_list)

        deltaF_o = sequence_o.scoring_function
        deltaF_m = sequence_m.scoring_function

        accept, accept_e = mc_delta(deltaF_o, deltaF_m, sequence_o.temp_shelf)

#        if accept == True and oligo == "on":
#           if sequence_m.oligomerization == True:
#                accept = False
        
#        if accept == True and (sequence_m.scoring_function < sequence_o.scoring_function):
#            best_score = sequence_m
            
        if accept ==True:
            sequence_o = sequence_m 
            worker_stats.update_acc_mc_step()
            if accept_e == True:
                worker_stats.update_acc_mc_better_e()
        '''    
        if accept == True and (sequence_o.scoring_function < best_score.scoring_function) and (worker_stats.global_step % 1000 !=0):
            best_score = sequence_o
        elif accept == True and (sequence_o.scoring_function >= best_score.scoring_function) and (worker_stats.global_step % 1000 ==0): 
            print("kupka")
            best_score = sequence_o
        '''

        if accept == False:
            worker_stats.update_rej_mc_step()
            
#        if best_score == '':
#            best_score = in_sequence

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
        
    return single_replica_design(*args);


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
        
        #inputs = [(seq_obj, nt_list, stats_obj) for seq_obj in lst_seq_obj]
        inputs = [(seq_obj, nt_list, stats_obj, seed) for seq_obj, seed in zip(lst_seq_obj, seeds)]
        lst_seq_obj_res_new = []
        workers_stats =[]

        for result in pool.imap(par_wrapper,inputs):
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
    
    input_file.pairs = check_dot_bracket(input_file.sec_struct)  #check dotbracket correctness, assign as list of pairs 
    check_seq_restr(input_file.seq_restr)
    check_length(input_file.sec_struct, input_file.seq_restr)
    nt_list = get_nt_list(input_file)
    check_input_logic(nt_list)
    
    
    global target_pairs_tupl
    target_pairs_tupl = {tuple(pair) for pair in input_file.pairs}
    
    seqence_score_list = generate_initial_list(nt_list, input_file)
    seqence_score_list_random = generate_initial_list_random(nt_list, input_file)

    start_time = time.time()
    
    
    seqence_score_list_new = seqence_score_list.copy()

    simulation_data=[]    

    for i in range(len(seqence_score_list)):
        simulation_data.append(vars(seqence_score_list[i]))
                
    stats = func.Stats()
    
    while time.time() - start_time < timlim:

        print('ETA', round((timlim - (time.time() - start_time)), 0), 'seconds', end='\r')
        print(vars(stats))
        
        stats.update_global_step()

        seqence_score_list, stats = mutate_sequence_re(seqence_score_list, nt_list, stats)


        seqence_score_list, stats = replica_exchange(seqence_score_list, stats)

        for i in range(len(seqence_score_list)):
            seqence_score_list[i].get_sim_step(stats.step)
            simulation_data.append(vars(seqence_score_list[i]))

        if stats.global_step % 10 == 0:
            remove_duplicated_results = {item['sequence']: item for item in simulation_data}
            simulation_data_noduplicates = list(remove_duplicated_results.values())
            if oligo == "on":
                sorted_results_mid = sorted(round_floats(simulation_data_noduplicates), key=lambda d: (-d['oligomerization'],-d['mcc'], -d['edesired_minus_mfe'], -d['mfe']), reverse = True)
            else:
                sorted_results_mid = sorted(round_floats(simulation_data_noduplicates), key=lambda d: (-d['mcc'], -d['edesired_minus_mfe'], -d['mfe']), reverse = True)
            sorted_results_mid = sorted_results_mid[:num_results]
#            sorted_results_mid = round_floats(simulation_data)
            with open(outname+'_mid_results.csv', 'w', newline='') as csvfile:
                fieldnames = sorted_results_mid[0].keys()  # header from keys of the first dictionary
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

                writer.writeheader()
                for data in sorted_results_mid:
                    writer.writerow(data)
        
        if stats.global_step == RE_steps:
            break
        
        
    sorted_dicts = sorted(round_floats(simulation_data), key=lambda d: (d['sim_step'], d['replica_num']))

    with open(outname+'_traj.csv', 'w', newline='') as csvfile:
        fieldnames = sorted_dicts[0].keys()  # header from keys of the first dictionary
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for data in sorted_dicts:
            writer.writerow(data)

    fasta_txt = ""

    for i in range(0, len(sorted_dicts)):
        fasta_txt += ">"+infile+"|"+now+"|"+str(sorted_dicts[i]["replica_num"])+"|"+str(sorted_dicts[i]["sim_step"])+"|"+str(sorted_dicts[i]["scoring_function"])+"\n"
        fasta_txt += str(sorted_dicts[i]["sequence"])+"\n"
        
    with open(outname+'_multifasta.fas', 'w') as fastafile:
        fastafile.write(fasta_txt)

    
    df = pd.DataFrame(simulation_data)

    df.sort_values(['sim_step', 'replica_num'], inplace=True)

    
    final_df = pd.DataFrame()

    header = ['']
    header2 = ['']
    
    df = df.drop_duplicates(subset=['sim_step', 'replica_num'])


    for metric in ['sequence', 'mfe', 'edesired_minus_mfe', 'subopt_e', 'esubopt_minus_mfe', 'precision', 'recall', 'mcc', 'sln_mfe','temp_shelf', 'scoring_function']:
    # For each metric, create a sub-DataFrame and pivot it
        sub_df = df[['sim_step', 'replica_num', metric]].pivot(index='sim_step', columns='replica_num', values=metric)
        sub_df.columns = [f'replica{col}' for col in sub_df.columns]
        header.extend([metric] * len(sub_df.columns))
        header2.extend([f'replica{col}' for col in range(1, len(sub_df.columns) + 1)])
        header.append('')
        header2.append('')
        # Concatenate the pivoted sub-DataFrame to the final DataFrame
        final_df = pd.concat([final_df, sub_df, pd.DataFrame(columns=[' '])], axis=1)  # Add an empty column here

    with open(outname+'_replicas.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow(header2)
    final_df.to_csv(outname+'_replicas.csv',mode='a', header=False)

    
    sorted_scores = sorted(round_floats(simulation_data), key=lambda d: (-d['scoring_function']), reverse = True)
    sorted_scores = sorted_scores[:num_results]
    
    best_fasta_txt = ""
    for i in range(0, len(sorted_scores)):
        best_fasta_txt += ">"+infile+"|"+now+"|"+str(sorted_scores[i]["replica_num"])+"|"+str(sorted_scores[i]["sim_step"])+"|"+str(sorted_scores[i]["scoring_function"])+"\n"
        best_fasta_txt += str(sorted_scores[i]["sequence"])+"\n"
        
    with open(outname+'_best_fasta.fas', 'w') as fastafile:
        fastafile.write(best_fasta_txt)

    remove_duplicated_results = {item['sequence']: item for item in simulation_data}
    simulation_data_noduplicates = list(remove_duplicated_results.values())
    if oligo == "on":
        sorted_results = sorted(round_floats(simulation_data_noduplicates), key=lambda d: (-d['oligomerization'], -d['mcc'], -d['edesired_minus_mfe'], -d['mfe']), reverse = True)
    else:
        sorted_results = sorted(round_floats(simulation_data_noduplicates), key=lambda d: (-d['mcc'], -d['edesired_minus_mfe'], -d['mfe']), reverse = True)
    sorted_results = sorted_results[:num_results]
    
    with open(outname+'_results.csv', 'w', newline='') as csvfile:
        fieldnames = sorted_results[0].keys()  # header from keys of the first dictionary
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for data in sorted_results:
            writer.writerow(data)


    sorted_results = sorted_results[:10]

    correct = 0
    correct_bool = False    
    for i in range(0,len(sorted_results)):
        if sorted_results[i]['mcc'] ==0.0:
            correct +=1
            correct_bool = True
    
    if oligo == "on":
        if sorted_results[0]['oligomerization'] == True:
            oligo_txt = ",oligomerize: yes"
        elif sorted_results[0]['oligomerization'] == False:
            oligo_txt = ",oligomerize: no"
        correct_result_txt = ">"+infile+","+str(correct_bool)+","+str(correct)+","+sorted_results[0]['sequence']+','+sorted_results[0]['mfe_ss']+","+input_file.sec_struct+oligo_txt+'\n'
    else:
        correct_result_txt = ">"+infile+","+str(correct_bool)+","+str(correct)+","+sorted_results[0]['sequence']+','+sorted_results[0]['mfe_ss']+","+input_file.sec_struct+'\n'
    
#    print("\n\nDesign solved: ",correct_bool)
    
    with open(outname+'_best_str', 'w', newline='') as myfile:
        myfile.write(correct_result_txt)

    func.plot_simulation_data_combined(simulation_data, outname,  input_file.alt_sec_struct, infile)


    traj_txt = ""
    
    sum_mc = stats.acc_mc_step + stats.rej_mc_step
    acc_perc  = round(stats.acc_mc_step/sum_mc, 3)
    acc_perc_metro = round((stats.acc_mc_step-stats.acc_mc_better_e)/sum_mc, 3)


    sum_mc_metro = sum_mc-stats.acc_mc_better_e
    acc_metro = stats.acc_mc_step-stats.acc_mc_better_e

    if correct_bool == True:
        best_solution_txt = "\nDesign solved succesfully!\n\nBest solution:\n"+str(sorted_results[0]['sequence'])+"\n"+str(sorted_results[0]['mfe_ss'])\
                            +"\nMFE:"+str(round(sorted_results[0]['mfe'],3))+"\n"
    else:
        best_solution_txt = "\nDesign not solved!"+"\n\nTarget structure:\n"+str(input_file.sec_struct)+"\n\nClosest solution:\n"+str(sorted_results[0]['sequence'])+"\n"+str(sorted_results[0]['mfe_ss'])\
                           +"\nMFE: "+str(round(sorted_results[0]['mfe'],3))+"\n1-MCC: "+str(round(sorted_results[0]['mcc'],3))+"\n"
    
    stats_txt = ">"+outname+" time="+str(timlim)+"s\n"+'Acc_ratio= '+str(acc_perc)+', Iterations='+str(stats.step)+', Accepted='\
                        +str(stats.acc_mc_step)+'/'+str(sum_mc)+', Rejected='+str(stats.rej_mc_step)+'/'+str(sum_mc) + '\n'\
                        +'Accepted Metropolis='+str(acc_metro)+'/'+str(sum_mc_metro)+', Rejected Metropolis='+str(sum_mc_metro-acc_metro)+'/'+str(sum_mc_metro) + '\n'\
                        +"Replica exchange attempts: "+str(stats.global_step)+"\nReplica swaps attempts: "+str(stats.acc_re_step+stats.rej_re_step)+"\nReplica swaps accepted: "+str(stats.acc_re_step)\
                        +"\nReplica swaps rejected: "+str(stats.rej_re_step) + "\nReplica exchange acc_ratio: "\
                        + str( round(stats.acc_re_step/(stats.acc_re_step+stats.rej_re_step) , 3))+'\n'+best_solution_txt
    
    print('\n' + stats_txt)

    with open(outname+'_stats', 'w', newline='\n') as myfile:
        myfile.write(stats_txt)
        
    files_to_move = ['_best_str', '_traj.csv', '_replicas.csv', '_best_fasta.fas','_multifasta.fas','_random.csv', '.command']
    
    for i in range(0, len(files_to_move)):
        move(outname+files_to_move[i], "trajectory_files/"+(outname+files_to_move[i]))
    
    
    
    
    
    
    
    
def get_outname(infile, replicas, timlim, acgu_percentages, pks, T_max, T_min, oligo, dimer, param, RE_attempt, scoring_f, point_mutations,  alt_ss, tshelves, in_seed, subopt):
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
    
    outname = infile.split(".")[0]+'_R'+str(replicas)+"_e"+str(RE_attempt)+"_t"+str(timlim)+"_pk"+str(pks)+"_ACGU"+str(acgu_percentages)+\
                "_Tmin"+str(T_min)+"_Tmax"+str(T_max)+"_p"+str(param)+"_SF"+str("_".join(scoring_f))+"_O"+(str(oligo))+"_D"+(str(dimer))+"_PM"+(str(point_mutations))
    return outname


class Nucleotide:
    """
    This class represents a nucleotide in a sequence.

    Attributes:
    -----------
    number : int
        The position of the nucleotide in the sequence.
    letters : list
        The possible letters that this nucleotide could be (A, C, G, or U).
    pairs_with : int
        The position of the nucleotide that this one pairs with.
    pair_letters : list
        The possible letters that the paired nucleotide could be.
    letters_allowed : list
        The letters that are allowed for this nucleotide, considering constraints.
    """

    def __init__(self, number):
        self.number = number
        self.letters = []
        self.pairs_with = None
        self.pair_letters = []    
        self.letters_allowed = None
               
    def add_letter(self, letter):
        #self.letters.append(letter)
        self.letters = letter
    
    def add_pair(self, pair):
        self.pairs_with = pair
    
    def add_pairing_l(self, pairing_l):
        self.pair_letters += pairing_l
        self.pair_letters = list(set(self.pair_letters))
    
    def add_allowed_l(self, list):
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
        self.name = name
        self.sec_struct = sec_struct
        self.seq_restr = seq_restr 
        self.pairs = []
        self.seed_seq = None
        self.alt_sec_struct = None

    def add_seed_seq(self, seed_seq):
        self.seed_seq = seed_seq

    def add_alt_sec_struct(self, alt_sec_struct):
        self.alt_sec_struct = alt_sec_struct



if __name__ == '__main__':


    print("DesiRNA")

    command = os.path.basename(sys.argv[0])+' '+' '.join(sys.argv[1:])
    
    script_path = os.path.dirname(os.path.abspath(__file__))


    now = datetime.now()
    now = now.strftime("%Y%m%d.%H%M%S")
    
    infile, replicas, timlim, acgu_percentages,  T_max, T_min, oligo,  param, RE_attempt, scoring_f, point_mutations,  \
    tshelves, in_seed, subopt, diff_start_replicas, num_results, acgu_content, RE_steps, tm_max, tm_min, motifs = argument_parser()
    
    input_file = read_input(infile)

    if RE_steps != None:
        timlim = 100000000000000000
    
    
    
    if param == '1999':
        RNA.params_load(os.path.join(script_path, "rna_turner1999.par"))

    L = 504.12
    
    if acgu_content == '':
        nt_percentages = {"A":15, "C":30, "G":30,"U":15}
    else:
        acgu_l = [int(x) for x in acgu_content.split(',')]
        if sum(acgu_l) != 100:
            print("The ACGU content should sum up to 100, check your command.")
            quit()        
        nt_percentages = {"A":acgu_l[0], "C":acgu_l[1], "G":acgu_l[2],"U":acgu_l[3]}
    

    if input_file.alt_sec_struct != None:
        alt_ss = "on"
        scoring_f = "alt"
    else:
        alt_ss = "off"    

    if "&" in input_file.sec_struct:
        too_much = input_file.sec_struct.count("&")
        if too_much == 1:
            dimer = "on"
        else:
            print("\nToo much structures in the input. Can only design dimers.\nPlease correct your input file.\n")
            quit()
    else:
        dimer = "off"

    if set(input_file.sec_struct) != set('.()'):
#    if ("[" or "<" or "{") in input_file.sec_struct:
        pks = "on"
    else:
        pks = "off"
    
    outname = get_outname(infile, replicas, timlim, acgu_percentages, pks, T_max, T_min, oligo, dimer, 
                          param, RE_attempt, scoring_f, point_mutations,  alt_ss, tshelves, in_seed, subopt)

    if tshelves == '':
        rep_temps_shelfs = get_rep_temps()
    else: 
        rep_temps_shelfs = [float(temp) for temp in tshelves.split(",")]

    with open(outname+".command", 'w') as f:
        print(command, file=f)
    
    if in_seed != 0:
        original_seed = 2137+in_seed
    else:
        original_seed = random.random()

#    if rand_seed == 'on':
#        original_seed = random.random()
    random.seed(original_seed)

    WORK_DIR = outname+"_"+str(now)
    Path(WORK_DIR).mkdir(parents=True, exist_ok=True)
    copy(infile, WORK_DIR+"/"+infile)
    move(outname+".command", WORK_DIR+"/"+(outname+".command"))
    Path(WORK_DIR+"/trajectory_files").mkdir(parents=True, exist_ok=True)
    os.chdir(WORK_DIR)


    
    run_functions()
    
    
    
