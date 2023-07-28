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
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime


from utils import functions_classes as func
from utils.SimScore import SimScore
from utils.pcofold_dimer_multichain_energy import mfe_e_dimer

def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__, prog='2DesignRNA')
    parser.add_argument("-f", "--filename", required=True, dest="name",
                        help="Name of a file that contains secondary structures and constraints.")
    parser.add_argument("-R", "--replicas", required=False, dest="replicas", default=10, type=int,
                            help="Number of replicas. [default = 10]")
    parser.add_argument("-e", "--exchange", required=False, dest="exchange", default=10, type=int,
                            help="TBA. [default = 10]")
    parser.add_argument("-t", "--timelimit", required=False, default=3600, dest="timlim", type=int,
                            help="Timelimit for running the program [s]. [default = 60s]")
    parser.add_argument("-n", "--number_of_sequences", required=False, dest="number_seq", default=10, type=int,
                            help="Number of best sequences to return. [default = 10]")
    parser.add_argument("-acgu", "--ACGU", required=False, dest="percs", default='off', choices=['off','on'],
                            help="TBA")
    parser.add_argument("-pk", "--PK", required=False, dest="pks", default='off', choices=['off','on'],
                            help="TBA")
    parser.add_argument("-tmin", "--tmin", required=False, dest="t_min", default=0.05, type=float,
                            help="Minimal Replica Temperature. [default = 2]")
    parser.add_argument("-tmax", "--tmax", required=False, dest="t_max", default=30, type=float,
                            help="Maximal Replica Temperature. [default = 620]")

    parser.add_argument("-o", "--oligomerization", required=False, dest="oligo", default='off', choices=['off','on'],
                            help="TBA")
    parser.add_argument("-d", "--dimer", required=False, dest="dimer", default='off', choices=['off','on'],
                            help="TBA")

    parser.add_argument("-p", "--param", required=False, dest="param", default='1999', choices=['2004','1999'],
                            help="TBA")
    parser.add_argument("-sf", "--scoring_function", required=False, dest="scoring_f", default='dmt', choices=['dmt','mcc', 'mfe', 'mix', 'mix2'],
                            help="TBA")
    parser.add_argument("-m", "--mutations", required=False, dest="mutations", default='one', choices=['one','multi'],
                            help="TBA")    
    parser.add_argument("-pm", "--point_mutations", required=False, dest="pm", default='off', choices=['off','on'],
                            help="TBA")

    parser.add_argument("-s", "--steps", required=False, default=100000000, dest="steps", type=int,
                            help="Number of steps")


                            
    args = parser.parse_args() 
    
    infile = args.name
    replicas = args.replicas
    timlim = args.timlim
    num_results = args.number_seq
    acgu_percs = args.percs
    pks = args.pks
    t_max = args.t_max
    t_min = args.t_min
    oligo = args.oligo
    dimer = args.dimer
    param = args.param
    exchange_rate = args.exchange
    scoring_f = args.scoring_f
    mutations = args.mutations
    pm = args.pm
    steps = args.steps
    
    
    return infile, replicas, timlim, num_results, acgu_percs, pks, t_max, t_min, oligo, dimer, param, exchange_rate, scoring_f, mutations, pm, steps


def read_input():
    
    string = ""
    with open(infile) as f:
        for line in f:
            string += line
            
    data = string.lstrip(">").rstrip("\n").split("\n>")
    data_dict = {}
    for i in range(0, len(data)):
        type = data[i].split("\n")
        key, values = type[0], type[1:]
        data_dict[key] = values
    
    if 'seed_seq' in data_dict:
        input = InputFile(data_dict['name'][0], data_dict['sec_struct'][0], data_dict['seq_restr'][0], data_dict['seed_seq'][0])
    else:
        input = InputFile(data_dict['name'][0], data_dict['sec_struct'][0], data_dict['seq_restr'][0], None)
    
    return input

'''
def psudoknot_specifier(structure):
    counter=0
    OB_list = ['(','[','<','{','A','B','C','D','E']
    for open_bracket in OB_list:
        if open_bracket in structure:
            counter+=1
    if counter >1:
        flag=1
    else:
        flag=0
    return (flag)
'''

def check_dot_bracket(ss):
    
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
    allowed_characters = 'ACGUWSMKRYBDHVN-&'
    for i in range(0, len(restr)):
        if restr[i] not in allowed_characters:
            sys.exit("Not allowed characters in sequence restraints. Check input file.")

def check_length(ss, restr):
    if len(ss) != len(restr):
        sys.exit("Secondary structure and sequence restraints are of different length. Check input file.")

def get_nt_list(input):

    pair_list = input.pairs
    restr_seq = input.seq_restr
    
    list_of_nt = []
    
    for i in range(0, len(restr_seq)):
        obj = Nucleotide(number = i)
        obj.add_letter(nt_dictionary(restr_seq[i]))
        list_of_nt.append(obj)

    for i in range(0, len(pair_list)):
        open = pair_list[i][0]
        close = pair_list[i][1]
        list_of_nt[open].add_pair(close)
        list_of_nt[close].add_pair(open)
        nt_open = list_of_nt[open]
        nt_close = list_of_nt[close]

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
    
    for i in range(0, len(nt_list)):
        if nt_list[i].pairs_with != None:
            if_pair(nt_list[i], nt_list[nt_list[i].pairs_with])
    
def if_pair(nt1, nt2):
    
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
    pair_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G' }
    
    nt2 = pair_dict[nt1]
    return nt2


def can_pair(nt):
    
    pair_dict = {'A': ['U'], 'U': ['G','A'], 'G': ['U','C'], 'C': ['G'] }
    
    pairing_l = pair_dict[nt]
    return pairing_l



def random_sequence_generator(nt_list,input):

    seq_l = list(input.seq_restr).copy()
    pair_list = sorted(input.pairs.copy())

    
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
            #print(nt_dictionary(seq_l[i]))
            seq_l[i] = random.choice(nt_dictionary(seq_l[i]))  


    result_sequence = ''.join(seq_l)

#    mutated_nucleotide_count = {nt: result_sequence.count(nt) for nt in ['A', 'U', 'C', 'G']}
    
    
    return result_sequence    

    
def initial_sequence_generator(nt_list, input):

    seq_l = list(input.seq_restr).copy()
    pair_list = sorted(input.pairs.copy())


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
#    result_sequence = "CCAGUGGAGGAAACAGGAAACAUUGCAGGAAACAGG"

    if input.seed_seq:
        result_sequence = input.seed_seq


    return result_sequence


def allowed_choice(allowed, percs):
    
    percs_allowed = percs.copy()
    
    #allowed=["C","U"]
    for nt in percs:
        if nt not in allowed:
            percs_allowed[nt] = 0
    
    
    list_percs_allowed = list(percs_allowed.values())
 
    list_percs_allowed = [i for i in list_percs_allowed if i != 0]
    
    return list_percs_allowed

def get_rep_temps():

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
    
    
    return rep_temps

def generate_initial_list(nt_list, input):
    
    sequence = initial_sequence_generator(nt_list, input)
    
#    T_min = 1
#    T_max = 313
    '''
    if replicas != 1:
        delta = (T_max-T_min)/(replicas-1)
    else:
        delta = (T_max-T_min)/2
    #rep_temps = np.arange(T_min, T_max, delta)
    rep_temps = []

    T_curr = T_min
 
    for i in range(replicas):
        if replicas != 1:
            rep_temps.append(round(T_curr))
            T_curr +=delta
        else:
            T_curr +=delta
            rep_temps.append(round(T_curr))
    print(rep_temps)
    '''

    seq_list = []
    
    for i in range(0,replicas):
        sequence_object = score_sequence(sequence)
        sequence_object.get_replica_num(i+1)
        sequence_object.get_temp_shelf(rep_temps_shelfs[i])
        sequence_object.get_sim_step(0)
        
        print(vars(sequence_object))
        seq_list.append(sequence_object)        
    
    print(seq_list)
    
    init_sequences_txt =""
    for i in range(0, len(seq_list)):
        init_sequences_txt += str(vars(seq_list[i]))+"\n"
    
    print(init_sequences_txt)
    
    return seq_list


def generate_initial_list_random(nt_list, input):

    seq_list = []
    for i in range(0,replicas):
        sequence = random_sequence_generator(nt_list, input)
        deltaF = delta_MFE_EOS(input.sec_struct, sequence)
        score = metropolis_score(310, deltaF)
        seq_list.append([sequence, score])

    seq_list.sort(key = lambda x: x[1], reverse = True)


    result_list = seq_list[:num_results]

    with open(outname+'_random.csv', 'w', newline='\n') as myfile:
         wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
         wr.writerows(result_list)

    return result_list

    


def get_seq_to_mutate(sequence_list, step):

#    print(sequence_list)
    sequence = sequence_list[step]
#    score = sequence_list[step][1]
#    rep_num = sequence_list[step][2]
    
    return sequence

def mutate_sequence_re(sequence_obj, nt_list):
    
    sequence_mutated_re = sequence_obj
#    rep_temp_shelfs = get_rep_temps()
    print(rep_temps_shelfs, "rep_temp_shelfs")

    part_size = len(rep_temps_shelfs) // 3

    rep_temps_shelfs_sliced = [rep_temps_shelfs[i*part_size:i*part_size + part_size if i != 2 else None] for i in range(3)]
    print(rep_temps_shelfs_sliced)

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

    if point_mutations == "off":
        mutation_position = random.choice(available_positions)
    elif point_mutations == "on":

        max_perc_prob = 0.7
        min_perc_prob = 0.0

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

        if len(false_negatives) == 0:
            range_pos = available_positions
        else:
            range_pos = random.choices([false_cases, available_positions], weights=[mutat_point_prob, 1-mutat_point_prob])[0]


        mutation_position = random.choice(range_pos)
        
    return mutation_position


def mutate_sequence(sequence_obj, nt_list):
    
    sequence = sequence_obj.sequence
    sequence_list = list(sequence)
    range_pos =[]
    for i in range(0, len(sequence)):
        if len(nt_list[i].letters_allowed) != 1:
            range_pos.append(i)


#    nt_pos = random.choice(range_pos)
    nt_pos = get_mutation_position(sequence_obj, range_pos)
    mmm = " "*len(sequence_obj.sequence)
    mmm = mmm[:nt_pos]+"*"+mmm[nt_pos:]
    print(mmm)
    print(sequence_obj.sequence)
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


def run_cons_maker(seq,s1,s2,s3):

    cons_struct = '.'*len(seq)
    
    name = outname+'_final_P'+str(replicas)+"_t"+str(timlim)+"_pk"+str(pks)

    cons_temp = open(name+'.cons','w')
    cons_temp.write(s1+'\n'+s2+'\n'+s3)
    cons_temp.close()    


    cmd = 'RNA_Consensus '+" "+name+".cons .5  | grep -A1 'consensus_pairs_unique_single_line' | tail -1 >"+name+".SS"
    os.system(cmd)

    cons_result = open(name+".SS",'r')
    cons_struct= cons_result.read().rstrip()
    cons_result.close()

    
    print(seq,cons_struct, "cons_maker")

    return cons_struct


def get_pk_struct(seq, ss_nopk):


    print(ss_nopk, "nopk")

    constraints = ss_nopk.replace("(","x").replace(")","x")

    fc = RNA.fold_compound(seq)
    
    fc.hc_add_from_db(constraints)

    mfe_structure, mfe_energy = fc.mfe()
    print(mfe_structure, "mfe pk")

    
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
        print(constraints, "cons pk2")
        fc = RNA.fold_compound(seq)
        fc.hc_add_from_db(constraints)

        mfe_structure, mfe_energy = fc.mfe()
        print(mfe_structure, "mfe pk2")
        l_ss_pk = list(ss_pk)
        l_mfe_pk2 = list(mfe_structure)
        
        for i in range(0, len(l_ss_pk)):
            if l_mfe_pk2[i] == "(":
                l_ss_pk[i] = "<"
            if l_mfe_pk2[i] == ")":
                l_ss_pk[i] = ">"        

        ss_pk = ''.join(l_ss_pk)
    
    print(ss_pk,"mfe_wpk")


#    return ipknot_struct.replace('\n','')
    return ss_pk


def delta_MFE_EOS_w_psk(seq):

#    dotknot_struct=run_dotknot(seq)
#    hotknot_struct=run_Hotknot(seq)
    ipknot_struct=run_ipknot(seq)
 
#    cons_struct=run_cons_maker(seq,dotknot_struct,hotknot_struct,ipknot_struct)

    F0 = energy_of_pk(seq, input.sec_struct)
    F1 = energy_of_pk(seq, ipknot_struct)
    print(F0, "\n", F1)

    
#    delta_F = abs(F0 - F1)
    delta_F = abs(F0 - F1)


    return delta_F, ipknot_struct



def delta_MFE_EOS(struct, seq):
    
    F0 = RNA.energy_of_struct(seq, struct)
    F1 = RNA.pf_fold(seq)[1]
    
        
    delta_F = F0 - F1
    
    return delta_F



def mc_delta(deltaF_o, deltaF_m,  T_replica):

#    mut_score = metropolis_standard_score( deltaF_m) 
    accept_d = False
    if deltaF_m <= deltaF_o:
        accept = True
        accept_d = True
#        accept = False
    else:
        diff= deltaF_m - deltaF_o
        p_dE = metropolis_score(T_replica, diff)
        rand_num = random.random()

        accept = p_dE > rand_num

#        if accept:
#            acc[0]+=1
#        else:
#            acc[1]+=1
#    print(accept)

    
    return accept, accept_d


def metropolis_score(temp, dE):
    
    score = math.exp((-L/temp) * (dE))
    
    return score



def metropolis_standard_score(dE):

    score = math.exp((-504.12/310) * (dE))

    return score



def replica_exchange_attempt(T0, T1, dE0, dE1):

    
    if dE1 <= dE0:
        accept = True
    else:
        rand_num = random.random()
#        p = math.exp(L*(1/T0 -1/T1)*(dE0-dE1))
#        T_re = 10
        p = math.exp(L*(1/T_re)*(dE0-dE1)) # bardziej rownomierne wymiany replik
        #p = math.exp((1/T_re)*(dE0-dE1))
#        print(L*(1/T0 -1/T1))
#        print(1/T_re)
#        print(p)
        #quit()
        accept = p > rand_num   


    return accept


def replica_exchange(replicas, c):

    sec_struct = input.sec_struct
    re_pairs = []
    stats = [0,0]

    num_shelfs = [i for i in range(1, len(replicas))]
    

    temps = []
    
    for i in range(0, len(replicas)):
        temps.append(replicas[i].temp_shelf)
    
    temps = sorted(temps)
    
    if c % (2*RE_attempt) ==0:
    
        for i in range(0, len(num_shelfs)-1, 2):
            re_pairs.append([i+1, i+2])
    else:
        for i in range(0, len(num_shelfs), 2):
            re_pairs.append([i, i+1])
    
    replicas = sorted(replicas, key=lambda obj: obj.temp_shelf)


    for i in range(len(re_pairs)):
        rep_i = re_pairs[i][0] #replica temporary number
        rep_j = re_pairs[i][1] #replica temporary number
        print(rep_i, rep_j)
    
        T_i = replicas[rep_i].temp_shelf
        T_j = replicas[rep_j].temp_shelf
        seq_i =replicas[rep_i].sequence
        seq_j =replicas[rep_j].sequence
        
        dE_i = replicas[rep_i].scoring_function
        dE_j = replicas[rep_j].scoring_function
        
        accept = replica_exchange_attempt(T_i, T_j, dE_i, dE_j)
        
        if accept == True:
            replicas[rep_i].get_temp_shelf(T_j)
            replicas[rep_j].get_temp_shelf(T_i)
            stats[0] += 1
        else:
            stats[1] += 1
        print("Replica exchange:" , accept, replicas[rep_i].replica_num, replicas[rep_j].replica_num, T_i, T_j)
    replicas = sorted(replicas, key=lambda obj: obj.replica_num)


    return replicas, stats



def get_mfe_e_ss(seq):


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


    return energy, structure

def score_sequence(seq):

    
    scored_sequence = func.ScoreSeq(sequence = seq)


    mfe_energy, mfe_structure = get_mfe_e_ss(seq)
    
    print(mfe_energy, mfe_structure)
    scored_sequence.get_mfe_e(mfe_energy)
    scored_sequence.get_mfe_ss(mfe_structure)


    scored_sequence.get_target_e(RNA.energy_of_struct(seq, input.sec_struct.replace("&","")))
    print(scored_sequence.target_e, "target_energy")

    scored_sequence.get_d_mfe_target(scored_sequence.mfe_e, scored_sequence.target_e)

    if dimer== "off":
        scored_sequence.get_subopt_ss(func.get_first_suboptimal_structure_and_energy(seq)[0])
        scored_sequence.get_subopt_e(func.get_first_suboptimal_structure_and_energy(seq)[1])
    elif dimer == "on":
        scored_sequence.get_subopt_ss(mfe_structure)
        scored_sequence.get_subopt_e(mfe_energy)


    scored_sequence.get_d_mfe_subopt(scored_sequence.mfe_e, scored_sequence.subopt_e)    
    
    ssc = SimScore(input.sec_struct.replace("&","Ee"), scored_sequence.mfe_ss.replace("&","Ee"))
    ssc.find_basepairs()
    ssc.cofusion_matrix()

    scored_sequence.get_precision(ssc.precision())
    scored_sequence.get_recall(ssc.recall())
    scored_sequence.get_mcc(ssc.mcc())
    
    scored_sequence.get_reb()
    scored_sequence.get_web()
    
    scored_sequence.get_reb_subopt()
#    scored_sequence.get_web_subopt()
    scored_sequence.get_d_mfe_subopt_norm()
    
    scored_sequence.get_score()
    
    scored_sequence.get_scoring_function(scoring_f)

    if oligo == "on":
        scored_sequence.get_oligomerization()

    return scored_sequence



def round_floats(obj):
    """Round float values in an object to 3 decimal places."""
    if isinstance(obj, float):
        return round(obj, 3)
    elif isinstance(obj, dict):
        return {k: round_floats(v) for k, v in obj.items()}
    return obj


def run_functions():
    
#    input = read_input()
    input.pairs = check_dot_bracket(input.sec_struct)  #check dotbracket correctness, assign as list of pairs 
    check_seq_restr(input.seq_restr)
    check_length(input.sec_struct, input.seq_restr)
    nt_list = get_nt_list(input)
    check_input_logic(nt_list)
    
    
    global target_pairs_tupl
    target_pairs_tupl = {tuple(pair) for pair in input.pairs}
    
    seqence_score_list = generate_initial_list(nt_list, input)
    seqence_score_list_random = generate_initial_list_random(nt_list, input)

    start_time = time.time()
    print(seqence_score_list)
    
    c =0
    graph= ""
    
    acc_ratio = [0,0,0] # accepted all, rejceted, accepted default(lower energy)
    
    seqence_score_list_new = seqence_score_list.copy()

    position =0
        
    simulation_data=[]    

    pp = 0

    replica_step = 1
    
    replica_stats = [0,0,0]
    
    for i in range(len(seqence_score_list)):
        simulation_data.append(vars(seqence_score_list[i]))
                

    while time.time() - start_time < timlim:
        c +=1

#        print(c,end="\r")
        print("\n",c, "step")
        print(position, "position")

#        sequence_o =  get_seq_to_mutate(seqence_score_list,position)
        sequence_o = seqence_score_list[position]

#        sequence_m = mutate_sequence(sequence_o, nt_list)
        
        
        if mutations == "one": 
            sequence_m = mutate_sequence(sequence_o, nt_list)
        else:
#            sequence_m = mutate_sequence(sequence_o, nt_list)

            sequence_m = mutate_sequence_re(sequence_o, nt_list)

        sequence_o.get_sim_step(replica_step)
        sequence_m.get_sim_step(replica_step)

        print(vars(sequence_o), "original")
        print(vars(sequence_m), "mutant")

        deltaF_o = sequence_o.scoring_function
        deltaF_m = sequence_m.scoring_function


        print(deltaF_o,"deltaF_o")
        print(deltaF_m, "deltaF_m")



        accept, accept_def = mc_delta(deltaF_o, deltaF_m, sequence_o.temp_shelf)
        
        if oligo == "on":
            if sequence_m.oligomerization == True:
                accept = False
#        score_mutant = [mutated_sequence, score]

#        score_mutant, accept, acc_ratio = scoring_function(seq_to_mutate, mutated_sequence, rep_num, acc_ratio)
#        score_mutant.append(rep_num)
        
#        print(seqence_score_list,"seq score")
#        if c== 301:
#            quit()
        

        if accept == True:
            print(True)
            print(acc_ratio[0])
            print(replica_step, "replica_setp")
            #if replica_step ==1:
            #    simulation_data.append(vars(sequence_o))
            #    print("kupka")
            #    quit()
            simulation_data.append(vars(sequence_m))
            """
            Ap = score_mutant[0].count("A")
            Cp = score_mutant[0].count("C")
            Gp = score_mutant[0].count("G")
            Up = score_mutant[0].count("U")
            """
            #lenseq= len(score_mutant[0])
            #score_mutant.append(["A:"+str(round(Ap/lenseq,3)), "C:"+str(round(Cp/lenseq,3)),"G:"+str(round(Gp/lenseq,3)),"U:"+str(round(Up/lenseq,3))])
            '''
            if pks == 'on':
                score_mutant.append(cons_m)
            '''
#            print(seqence_score_list)
#            print(seqence_score_list[position])
            seqence_score_list[position] = sequence_m
            if replica_step % RE_attempt == 0 and position == replicas-1:
                print("\nAttempting Replica Exchange\n")
                
                seqence_score_list, curr_stats = replica_exchange(seqence_score_list,replica_step)
                replica_stats[0] +=1
                replica_stats[1] +=curr_stats[0]
                replica_stats[2] +=curr_stats[1]
#            print(seqence_score_list)
#            quit()
#            sequence_m.get_deltafm(deltaF_m)
#            if c==1:
#                simulation_data.append(vars(sequence_o))
#            simulation_data.append(vars(sequence_m))
#            sequence_o = sequence_m
            acc_ratio[0]+= 1
            position+=1
            if accept_def == True:
                acc_ratio[2] +=1
            if position == replicas:
                position = 0
                replica_step += 1
            if "[" in sequence_m.mfe_ss:
                pp+=1          


        else:
            print(False)
            print(acc_ratio[0])
            acc_ratio[1]+=1
                    
        #if replica_step % RE_attempt == 0:
        #    print("\nAttempting Replica Exchange\n")
        #    seqence_score_list = replica_exchange(seqence_score_list,replica_step, input.sec_struct)        
        

        if accepted_steps == acc_ratio[0]:
            break
        
    sum_mc = acc_ratio[0] + acc_ratio[1]
        

    
    # convert objects to dictionaries using vars(), and sort by 'mcc' property
#    print(simulation_data[0])

    # write dictionaries to csv file
    if len(simulation_data) > replicas*3:
        simulation_data = simulation_data[:-replicas*2]

    print(simulation_data, "s data")   
    sorted_dicts = sorted(round_floats(simulation_data), key=lambda d: (d['sim_step'], d['replica_num']))

    with open(outname+'_traj.csv', 'w', newline='') as csvfile:
        fieldnames = sorted_dicts[0].keys()  # header from keys of the first dictionary
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for data in sorted_dicts:
            writer.writerow(data)


    import pandas as pd
    import numpy as np


    fasta_txt = ""

    for i in range(0, len(sorted_dicts)):
        fasta_txt += ">"+infile+"|"+now+"|"+str(sorted_dicts[i]["replica_num"])+"|"+str(sorted_dicts[i]["sim_step"])+"|"+str(sorted_dicts[i]["scoring_function"])+"\n"
        fasta_txt += str(sorted_dicts[i]["sequence"])+"\n"
        
    with open(outname+'_multifasta.fas', 'w') as fastafile:
        fastafile.write(fasta_txt)


    
    # Convert list of dicts to pandas DataFrame

# Convert list of dicts to DataFrame
    
    df = pd.DataFrame(simulation_data)

# Sort DataFrame
    df.sort_values(['sim_step', 'replica_num'], inplace=True)

# Initialize an empty DataFrame for the final result
    
    final_df = pd.DataFrame()

    header = ['']
    header2 = ['']
    
#    duplicates = df[df.duplicated(subset=['sim_step', 'replica_num'], keep=False)]
#    print(duplicates.sort_values(by=['sim_step', 'replica_num']))
    df = df.drop_duplicates(subset=['sim_step', 'replica_num'])


    for metric in ['mfe_e', 'd_mfe_target', 'subopt_e', 'd_mfe_subopt', 'precision', 'recall', 'mcc', 'reb','temp_shelf', 'web', 'reb_subopt', 'd_mfe_subopt_norm', 'score', 'scoring_function']:
    # For each metric, create a sub-DataFrame and pivot it
        sub_df = df[['sim_step', 'replica_num', metric]].pivot(index='sim_step', columns='replica_num', values=metric)
        sub_df.columns = [f'replica{col}' for col in sub_df.columns]
        header.extend([metric] * len(sub_df.columns))
        header2.extend([f'replica{col}' for col in range(1, len(sub_df.columns) + 1)])
        header.append('')
        header2.append('')
    # Concatenate the pivoted sub-DataFrame to the final DataFrame
        #final_df = pd.concat([final_df, sub_df], pd.DataFrame(columns=[' '])], axis=1))
        final_df = pd.concat([final_df, sub_df, pd.DataFrame(columns=[' '])], axis=1)  # Add an empty column here

    with open(outname+'_replicas.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow(header2)
# Save final DataFrame to csv
    final_df.to_csv(outname+'_replicas.csv',mode='a', header=False)

    
    # Save to CSV

    sorted_scores = sorted(round_floats(simulation_data), key=lambda d: (-d['scoring_function']), reverse = True)
    sorted_scores = sorted_scores[:100]
    
    best_fasta_txt = ""
    for i in range(0, len(sorted_scores)):
        best_fasta_txt += ">"+infile+"|"+now+"|"+str(sorted_scores[i]["replica_num"])+"|"+str(sorted_scores[i]["sim_step"])+"|"+str(sorted_scores[i]["scoring_function"])+"\n"
        best_fasta_txt += str(sorted_scores[i]["sequence"])+"\n"
        
    with open(outname+'_best_fasta.fas', 'w') as fastafile:
        fastafile.write(best_fasta_txt)

    sorted_results = sorted(round_floats(simulation_data), key=lambda d: (-d['mcc'], -d['d_mfe_target'], -d['mfe_e']), reverse = True)
    sorted_results = sorted_results[:100]
    
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
    
    correct_result_txt = ">"+infile+","+str(correct_bool)+","+str(correct)+","+sorted_results[0]['sequence']+','+sorted_results[0]['mfe_ss']+","+input.sec_struct+'\n'
    
    print("\nDesign solved: ",correct_bool)
    
    with open(outname+'_best_str', 'w', newline='') as myfile:
        myfile.write(correct_result_txt)
    
    ####

#    sortedlist = sorted(seqence_score_list, key=lambda x: x[1], reverse=True)
    
#    seqence_score_list = sortedlist
    
#    for i in range(0,len(sortedlist)):
#        if "e" in str(sortedlist[i][1]):
#            sortedlist[i][1] = '{:.3e}'.format(sortedlist[i][1])
#        else:
#            sortedlist[i][1] = str(round(sortedlist[i][1],5))

    if sum_mc ==0:
        sum_mc =1
    
    acc_perc  = round(acc_ratio[0]/sum_mc, 3)
    acc_perc_metro = round((acc_ratio[0]-acc_ratio[2])/sum_mc, 3)

#    with open(outname+'_final_P'+str(replicas)+"_t"+str(timlim)+"_pk"+str(pks)+'.csv', 'w', newline='\n') as myfile:
#    with open(outname+'_final.csv', 'w', newline='\n') as myfile:
#         myfile.write(">"+outname+" time="+str(timlim)+"s\n"+'Acc_ratio= '+str(acc_perc)+', Iterations='+str(c)+', Accepted='\
#                         +str(acc_ratio[0])+'/'+str(sum_mc)+', Rejected='+str(acc_ratio[1])+'/'+str(sum_mc)+'\n"'+input.sec_struct+'"\n"'+input.seq_restr+'"\n')
#         wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
#         wr.writerows(seqence_score_list)

#    compare_best_seq_to_struct(seqence_score_list, input.sec_struct, '')

#    print(">"+outname+" time="+str(timlim)+"s\n"+'Acc_ratio= '+str(acc_perc)+', Iterations='+str(c)+', Accepted='\
#                         +str(acc_ratio[0])+'/'+str(sum_mc)+', Rejected='+str(acc_ratio[1])+'/'+str(sum_mc))

#    outname_plot = outname+'_final_P'+str(replicas)+"_t"+str(timlim)+"_pk"+str(pks) 
    func.plot_simulation_data(simulation_data, outname, scoring_f)
    func.plot_simulation_data_re(simulation_data, outname)
    func.plot_simulation_data_combined(simulation_data, outname, scoring_f)


    traj_txt = ""
    
#    for i in range(0, len(simulation_data)):
#        for key in simulation_data[i]:
#            print(simulation_data[i][key], type(simulation_data[i][key]))
#            if isinstance(simulation_data[i][key], float):
#                simulation_data[i][key] = round(simulation_data[i][key], 3) 
#        traj_txt += str(simulation_data[i])+"\n"
        
#    with open(outname+'_traj.csv', 'w', newline='\n') as myfile:
#        myfile.write(str(traj_txt).replace(":",",").replace("}","").replace("{",""))
    
#    print(">"+outname+" time="+str(timlim)+"s\n"+'Acc_ratio= '+str(acc_perc)+', Iterations='+str(c)+', Accepted='\
#                              +str(acc_ratio[0])+'/'+str(sum_mc)+', Rejected='+str(acc_ratio[1])+'/'+str(sum_mc))
#    print(replica_stats,"re stats")
    
#    print(str(round(replica_stats[1]/(replica_stats[1]+replica_stats[2]), 3)))
    if replica_stats[1]+replica_stats[2] == 0:
        replica_stats[1] = 999
    sum_mc_metro = sum_mc-acc_ratio[2]
    acc_metro = acc_ratio[0]-acc_ratio[2]
    stats_txt = ">"+outname+" time="+str(timlim)+"s\n"+'Acc_ratio= '+str(acc_perc)+', Iterations='+str(c)+', Accepted='\
                        +str(acc_ratio[0])+'/'+str(sum_mc)+', Rejected='+str(acc_ratio[1])+'/'+str(sum_mc) + '\n'\
                        +'Accepted Metropolis='+str(acc_metro)+'/'+str(sum_mc_metro)+', Rejected Metropolis='+str(sum_mc_metro-acc_metro)+'/'+str(sum_mc_metro) + '\n'\
                        +"Replica exchange attempts: "+str(replica_stats[0])+"\nReplica swaps accepted: "+str(replica_stats[1])\
                        +"\nReplica swaps rejected: "+str(replica_stats[2]) + "\nReplica exchange acc_ratio: "\
                        + str( round(replica_stats[1]/(replica_stats[1]+replica_stats[2]) , 3))+'\n'
    
    print('\n' + stats_txt)
    with open(outname+'_stats', 'w', newline='\n') as myfile:
        myfile.write(stats_txt)
        '''
        myfile.write(">"+outname+" time="+str(timlim)+"s\n"+'Acc_ratio= '+str(acc_perc)+', Iterations='+str(c)+', Accepted='\
                        +str(acc_ratio[0])+'/'+str(sum_mc)+', Rejected='+str(acc_ratio[1])+'/'+str(sum_mc) + '\n')
        myfile.write("Replica exchange attempts: "+str(replica_stats[0])+"\nReplica swaps accepted: "+str(replica_stats[1])\
                        +"\nReplica swaps rejected: "+str(replica_stats[2]) + "\nReplica exchange acc_ratio: "+ str( round(replica_stats[1]/(replica_stats[1]+replica_stats[2]) , 3))+'\n' )
        '''
        
        
def get_outname():

    outname = infile.split(".")[0]+'_R'+str(replicas)+"_e"+str(RE_attempt)+"_t"+str(timlim)+"_pk"+str(pks)+"_ACGU"+str(acgu_percentages)+\
                "_Tmin"+str(T_min)+"_Tmax"+str(T_max)+"_p"+str(param)+"_SF"+str(scoring_f)+"_O"+(str(oligo))+"_D"+(str(dimer))+"_M"+(str(mutations))+"_PM"+(str(point_mutations))
#    print(outname)
#    quit()
    return outname


def compare_best_seq_to_struct(seq_score_list, input_ss, suffix):
    
    name = outname
    out_fa = name+'_best_seq'+suffix+'.seq'
    out_str = name+suffix+'.best_str'
    

    num_solved = 0

    for i in range(len(seq_score_list)-1, -1,-1):
    
        best_seq = seq_score_list[i][0]
        score = seq_score_list[i][1]
    
        with open(out_fa, 'w')as myfile:
            myfile.write(best_seq)


        cmd = "RNAfold -p -d2 < %s" % (out_fa)
        rnafold = os.popen(cmd).read().splitlines()[1].split(' ')[0]
        #if pks == 'on':
        #    rnafold = seq_score_list[i][-1]
 
        same = False
    
        if rnafold == input_ss:
            same = True
            num_solved += 1
        else:
            same = False
        
    
    if num_solved > 0:
        same = True
    else:
        same = False
        
    with open(out_str, 'w') as file:
        file.write(">"+outname+","+str(same)+","+str(score)+","+str(num_solved)+","+best_seq+","+input_ss+","+rnafold+"\n")
    
        
#    cmd = "rm *ps *.seq *efn2  *ct *SS *.cons  *ipknot *ipknot.fas"
    cmd = "rm *ps *.seq"
    os.system(cmd)



class Nucleotide:
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
    def __init__(self, name, sec_struct, seq_restr, seed_seq):
        self.name = name
        self.sec_struct = sec_struct
        self.seq_restr = seq_restr 
        self.seed_seq = seed_seq
        self.pairs = []





if __name__ == '__main__':


    print("DesiRNA")
    print(RNA.__version__)
#    print(os.path.basename(sys.argv[0])+' '+' '.join(sys.argv[1:]))
    command = os.path.basename(sys.argv[0])+' '+' '.join(sys.argv[1:])
    
#    script_path = os.path.dirname(__file__)
    script_path = os.path.dirname(os.path.abspath(__file__))

    print(script_path)

    now = datetime.now()
    now = now.strftime("%Y%m%d.%H%M%S")
    
    infile, replicas, timlim, num_results, acgu_percentages, pks, T_max, T_min, oligo, Dimer, param, RE_attempt, scoring_f, mutations, point_mutations, accepted_steps = argument_parser()

    if param == '1999':
#        RNA.params_load_RNA_Turner1999()
        RNA.params_load(script_path+"/rna_turner1999.par")


    T_re = 10

    rep_temps_shelfs = get_rep_temps()
    print(rep_temps_shelfs)
    
    L = 504.12
    #L = 5000
    #L = t_lambda
    
    
    print(acgu_percentages)

    nt_percentages = {"A":15, "C":30, "G":30,"U":15}
    
    input = read_input()


    if "&" in input.sec_struct:
        dimer = "on"
    else:
        dimer = "off"


    if "[" in input.sec_struct:
        pks = "on"
    else:
        pks = "off"
    
    outname = get_outname()
    
    print(command)
    with open(outname+".command", 'w') as f:
        print(command, file=f)

    
    random.seed(2137)
#    if num_results != replicas:
#        num_results = replicas
    
    run_functions()
    
    
    