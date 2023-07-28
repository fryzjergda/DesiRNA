#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import random
import os

def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__, prog='generate_random_sequences.py', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-n", "--num", required=False, dest="num", default = 100, type = int,
                        help="Number of sequences to generate.")
    parser.add_argument("-min", "--minimum", required=False, dest="min", default = 10, type = int,
                        help="Minimal length of generated sequence.")
    parser.add_argument("-max", "--maximum", required=False, dest="max", default = 100, type = int,
                        help="Maximal length of generated sequence.")
    parser.add_argument("-s", "--starting_number", required=False, dest="start", default = 1, type = int,
                        help="Maximal length of generated sequence.")


    args = parser.parse_args() 

    number = args.num
    min = args.min
    max = args.max
    starting_num = args.start - 1

    return number, min, max, starting_num
   

def random_proportions():
    
    low_cut = 10
    high_cut = 35
    
    probs = []
    for i in range(0,4):
        prob = random.randrange(low_cut,high_cut)
        probs.append(prob)


    return probs

def predict_ss(sequence_ungapped):

    create_tempfile(sequence_ungapped)
    cmd = "RNAfold -p -d2 --noLP < seq.temp" 
    ss = os.popen(cmd).read().splitlines()[2].split(' ')[0]

    return ss


def create_tempfile(segment):


    temp = open('seq.temp','w')
    temp.write(">temp\n")
    temp.write(segment)
    temp.close()


def to_stk(stk_list):

    lenline = 0
    new_list = []
    for i in range(0, len(stk_list)):
        if len(stk_list[i]) == 2:
            if len(stk_list[i][0]) > lenline:
                lenline = len(stk_list[i][0])
            str = "&"+stk_list[i][0]+' '+stk_list[i][1]
            #print(str)
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
            #print(name, len(name))
            #print(lenline)
            sequence = splitted[1]
            spaces = " "
            if len(name) < lenline:
                spaces = (lenline-len(name))*" "
                name += spaces
            newstring += "$"+name.replace("_PP"+spaces,spaces+"PP ")+" " + sequence + "\n"

    stk = newstring.rstrip().replace("$\n&$","\n\n").replace("$","").replace("&","")+ "\n//" #.replace("_"," ").replace("^","_") + "\n//"
    return stk


def convert_to_stk(stk_out):

    result = open(outname+".sto","w")
    result.write(stk_out)
    result.close()


def write_input(name, seq, ss):
    
    
    number = random.randrange(10, 1000)
    restr = "N" * len(seq)
    
    result = open(name+".input", "w")
    result.write(">name\n")
    result.write(name+"\n")
    result.write(">sec_struct\n")
    result.write(ss+"\n")
    result.write(">seq_restr\n")
    result.write(restr+"\n")
    result.write(">number\n")
    result.write(str(number)+"\n")
    
    
    
def generate_sequences():
    
    seq_list = []
    letters = ["A","G","C","U"]
    
    for i in range(1, num_seq+1):
        name = "RF_rand_" + str(i+starting_num)
        seq_len = random.randrange(min_len,max_len)
        probs = random_proportions()
        sequence = ""
        for k in range(0, seq_len):
            letter = random.choices(letters, probs)
            sequence += ''.join(letter)
        
        predicted = predict_ss(sequence)
        name_pred = "#=GR_"+name+"_predicted_SS"
        data = [name, sequence, name_pred, predicted]
        seq_list.append(data)
        write_input(name, sequence, predicted)    
        print("Step " + str(i) + "/" + str(num_seq), end = "\r", flush=True)
    return seq_list



if __name__ == '__main__':

    outname = "test"
    num_seq, min_len, max_len, starting_num = argument_parser()

    sequence_list = generate_sequences()
    
    stk_out = to_stk(sequence_list)
    
    convert_to_stk(stk_out)


    print("\nSuccesfully generated " + str(num_seq)+" random family design inputs!\nJOB GOOD!\n")    