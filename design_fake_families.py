#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse as argparse
import random
import sys

def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__, prog='design_fake_families.py')
    parser.add_argument("-f", "--filename", required=True, dest="name",
                        help="Name of a file that contains secondary structures and constraints.")
                            
    args = parser.parse_args() 
    
    infile = args.name

    return infile


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
    
    input = InputFile(data_dict['name'][0], data_dict['sec_struct'][0], data_dict['seq_restr'][0], data_dict['number'][0])

    return input


def check_dot_bracket(ss):
    
    db_list = [['(',')'],['[',']'],['<','>'],['{','}'],['A','a'],['B','b'],['C','c'],['D','d'],['E','e']]
    allowed_characters = '()[]<>{}AaBbCcDdEe.'
    
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
    allowed_characters = 'ACGUWSMKRYBDHVN-'
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
            'U': ['U']
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
    

def initial_sequence_generator(nt_list, input):

    seq_l = list(input.seq_restr).copy()
    pair_list = sorted(input.pairs.copy())

    # paired nucleotides handling, if possible GC, if not AU
    for i in range(0, len(pair_list)):
        nt1 = nt_list[pair_list[i][0]].letters_allowed
        nt2 = nt_list[pair_list[i][1]].letters_allowed
        nt1num =nt_list[pair_list[i][0]].number 
        nt2num =nt_list[pair_list[i][1]].number
        seq_l[nt1num] = random.choice(["C","G","A","U"])
        seq_l[nt2num] = wc_pair(seq_l[nt1num])

    # random choice (from allowed letters) of whatever left    
    for i in range(0, len(seq_l)):
        if seq_l[i] not in ["A","C","G","U"]:
            seq_l[i] = random.choice(nt_dictionary(seq_l[i]))	

    # introduce gaps
    perc_of_gaps = (random.randrange(0,40))/100
    num_of_gaps = int(len(seq_l)*perc_of_gaps)
    for i in range(0, num_of_gaps):
        position = random.randrange(0, len(seq_l))
        seq_l[position] = "-"
    
    result_sequence = ''.join(seq_l)
    
    return result_sequence

    
def convert_to_stk(stk_out):

    result = open(infile.replace("input","output"),"w")
    result.write(stk_out)
    result.close()


def to_stk(stk_list):

    lenline = 0
    new_list = []
    for i in range(0, len(stk_list)):
        if len(stk_list[i]) == 2:
            if len(stk_list[i][0]) > lenline:
                lenline = len(stk_list[i][0])
            str = "&"+stk_list[i][0]+' '+stk_list[i][1]
            new_list.append(str)
        elif len(stk_list[i]) ==4:
            if len(stk_list[i][0]) > lenline:
                lenline = len(stk_list[i][0])
            if len(stk_list[i][2]) > lenline:
                lenline = len(stk_list[i][2])
            str1 = "&"+stk_list[i][0]+' '+stk_list[i][1]
            str2 = "&"+stk_list[i][2]+' '+stk_list[i][3]
            new_list.append(str1)
            new_list.append(str2)

    newstring = ""

    for i in range(0, len(new_list)):
        if new_list[i][0:4] == "&#_S":
            newstring += new_list[i].replace("&","")+"\n"
        elif (new_list[i][0:5] == "&#=GF") or (new_list[i][0:5] == "&#=GS"):
            newstring += new_list[i].replace("&","").replace(" ","   ")+"$\n&"
        elif len(new_list[i]) != 0:
            splitted = new_list[i].split(" ")
            name =  splitted[0].replace("&","")
            sequence = splitted[1]
            spaces = " "
            if len(name) < lenline:
                spaces = (lenline-len(name))*" "
                name += spaces
            newstring += "$"+name.replace("_PP"+spaces,spaces+"PP ")+" " + sequence + "\n"

    stk = newstring.rstrip().replace("$\n&$","\n\n").replace("$","").replace("&","") + "\n//"
    return stk


def generate_initial_list(nt_list, input):
    
    seq_list = []
    for i in range(0,input.number):
        sequence = initial_sequence_generator(nt_list, input)
        seq_list.append(sequence)

    result_list = []
    for i in range(0, len(seq_list)):
        no = str(i+1)
        name = infile.replace("input","seq."+no)
        result_list.append([name, seq_list[i]])

    result_list.append(["#=GC_SS_cons", input.sec_struct])
    stk_out = to_stk(result_list)
    convert_to_stk(stk_out)

    return result_list
    
    
def get_random_seq(sequence_list):
    
    position = random.randrange(num_results)
    sequence = sequence_list[position][0]
    score = sequence_list[position][1]
    prob = sequence_list[position][2]
    return sequence, score, prob, position


def run_functions():
    input = read_input()
    input.pairs = check_dot_bracket(input.sec_struct)  #check dotbracket correctness, assign as list of pairs 
    check_seq_restr(input.seq_restr)
    check_length(input.sec_struct, input.seq_restr)
    nt_list = get_nt_list(input)
    check_input_logic(nt_list)
    seqence_score_list = generate_initial_list(nt_list, input)
    
    print("\nGenerated fake family: " + infile.replace("input","output") + "\n")

class Nucleotide:
    def __init__(self, number):
        self.number = number
        self.letters = []
        self.pairs_with = None
        self.pair_letters = []    
        self.letters_allowed = None
               
    def add_letter(self, letter):
        self.letters = letter
    
    def add_pair(self, pair):
        self.pairs_with = pair
    
    def add_pairing_l(self, pairing_l):
        self.pair_letters += pairing_l
        self.pair_letters = list(set(self.pair_letters))
    
    def add_allowed_l(self, list):
        self.letters_allowed = list
        
        
class InputFile:
    def __init__(self, name, sec_struct, seq_restr, number):
        self.name = name
        self.sec_struct = sec_struct
        self.seq_restr = seq_restr 
        self.pairs = []
        self.number = int(number)


if __name__ == '__main__':

    infile  = argument_parser()
    run_functions()
    