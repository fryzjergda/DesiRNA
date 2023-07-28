#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import RNA
import sys
import math
import pandas as pd
from utils.pcofold_dimer_multichain_energy import *
from utils.di_nuc_analyser import *
"""
Functions:
1 check_string: Checks whether the string is single or multi-chain.
2 add_ampersand: Maps the predicted secondary structure of a multi-chain sequence to the given sequence.
3 multi_chain_energy_reporter: Reports the energy and predicted secondary structure of a multi-chain sequence.
4 dimer_multi_chain_energy_structure_flags: Reports the structure and energy of dimers and multi-chains, adding the output as additional columns to the given output CSV file. It can also be added to the input CSV file.
5 get_first_suboptimal_structure_and_energy: Reports the most probable suboptimal secondary structure and its minimum free energy (MFE).
"""
import RNA


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec



import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

def plot_simulation_data_combined_(simulation_data, outname, scoring_f):
    # Define line colors
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'teal', 'magenta', 'cyan', 'gray', 'black', 'lime', 'pink', 'olive', 'navy']

    # Create subplots
    fig, axs = plt.subplots(4, 3, figsize=(24, 16))
    plt.subplots_adjust(top=0.9, bottom=0.05, left=0.1, right=0.95, hspace=0.4, wspace=0.3)  # Adjust spacing

    gs = GridSpec(4, 3, figure=fig, width_ratios=[1, 1, 1.3])
    fig.suptitle(outname, fontsize=16, fontweight='bold')  # Add title

    # Get unique replica numbers
    replica_nums = set(data['replica_num'] for data in simulation_data)
    num_replicas = len(replica_nums)

    # Group data by metric
    metrics = [
        'mfe_e', 'target_e', 'd_mfe_target', 'd_mfe_subopt', 'precision', 'recall',
        'mcc', 'score', 'reb', 'web', 'reb_subopt', 'd_mfe_subopt_norm', 'scoring_function'
    ]

    for i, metric_name in enumerate(metrics):
        # Create the first y-axis for the current metric
        ax = axs[i // 3, i % 3]
        ax.set_xlabel('Accepted Structures')
        ax.set_ylabel(metric_name)
        ax.set_title(metric_name)

        # Iterate over replicas and plot data
        for j, replica_num in enumerate(replica_nums):
            replica_data = [data for data in simulation_data if data['replica_num'] == replica_num]

            # Extract values for the current replica and metric
            metric_values = [data[metric_name] for data in replica_data]
            accepted_structures = range(len(metric_values))

            # Plot data for the current replica
            ax.plot(accepted_structures, metric_values, color=colors[j], label=f'Replica {replica_num}')

        ax.legend()

    # Save the plot to a file
    plt.savefig(outname + '.png')


import matplotlib.pyplot as plt
import matplotlib.cm as cm

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import matplotlib.pyplot as plt

def get_replica_color(replica_num, num_replicas):
    # Define a list of ten colors ranging from blue to red
    colors = ['#0000FF', '#00BFFF', '#00FF00', '#7CFC00', '#FFFF00',
              '#FFD700', '#FFA500', '#FF0000', '#8B0000', '#B22222']
    
    # Calculate the index of the color for the given replica
    if num_replicas == 1:
        color_index = 0
    else:
        color_index = int((replica_num - 1) * (len(colors) - 1) / (num_replicas - 1))
    
    return colors[color_index]

def plot_simulation_data_combined(simulation_data, outname, scoring_f):
    # Define the metrics to plot
    metrics = ['mfe_e', 'precision', 'reb', 'target_e', 'recall', 'web',
               'd_mfe_target', 'mcc', 'temp_shelf', 'd_mfe_subopt',
               'score', 'scoring_function']

    # Create subplots
    num_plots = len(metrics)
    num_cols = 3
    num_rows = (num_plots + num_cols - 1) // num_cols
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(24, 16))

    # Flatten the axs if there's only one row
    if num_rows == 1:
        axs = [axs]

    scoring_function_titles = {
    'dmt': 'Scoring Function - d_mfe_target',
    'mcc': 'Scoring Function - MCC',
    'mix': 'Scoring Function - Score',
    'mix2': 'Scoring Function - Score2',
    'mfe': 'Scoring Function - MFE'}



    # Get the number of replicas
    num_replicas = len(set(data['replica_num'] for data in simulation_data))

    # Iterate over replicas and plot data
    for i, metric in enumerate(metrics):
        row = i // num_cols
        col = i % num_cols
        ax = axs[row][col]

        for replica_num in range(1, num_replicas + 1):
            replica_data = [data[metric] for data in simulation_data if data['replica_num'] == replica_num]
            accepted_structures = range(len(replica_data))
            color = get_replica_color(replica_num, num_replicas)  # Get color for the replica
            ax.plot(accepted_structures, replica_data, color=color, label=f"Replica {replica_num}")


        ax.set_xlabel('Accepted Structures')
        ax.set_ylabel(metric.capitalize().replace('_', ' '))
        ax.set_title(metric.capitalize().replace('_', ' '))
        ax.legend(loc='upper right')


        if metric == 'mcc':
            ax.set_title('1 - MCC')
        if metric == 'precision':
            ax.set_title('1 - Precision')
        if metric == 'recall':
            ax.set_title('1 - Recall')

        if metric == 'scoring_function':
            if scoring_f in scoring_function_titles:
                ax.set_title(scoring_function_titles[scoring_f])


    # Add overall title to the plot
    fig.suptitle(outname, fontsize=16, fontweight='bold')

    # Adjust the layout
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save the plot to a file
    plt.savefig(outname + '_replicas.png')

# Example usage:
# plot_simulation_data_combined(simulation_data, 'combined_metrics_plot', 'dmt')




def plot_simulation_data_re(simulation_data, outname):
    replica_nums = set(data['replica_num'] for data in simulation_data)
    num_replicas = len(replica_nums)

    # Create the figure and subplots with fixed layout
    fig, axs = plt.subplots(4, 3, figsize=(24, 16))
    plt.subplots_adjust(top=0.9, bottom=0.05, left=0.1, right=0.95, hspace=0.4, wspace=0.3)  # Adjust spacing

    gs = GridSpec(4, 3, figure=fig, width_ratios=[1, 1, 1.3])

    fig.suptitle(outname, fontsize=16, fontweight='bold')  # Add title


    # Define line colors for each replica
    colors = plt.cm.tab10(range(num_replicas))

    # Iterate over replicas and plot data
    for i, replica_num in enumerate(replica_nums):
        replica_data = [data for data in simulation_data if data['replica_num'] == replica_num]

        # Extract values for the current replica
        accepted_structures = range(len(replica_data))
        scoring_function_values = [data['scoring_function'] for data in replica_data]
        temp_shelf_values = [data['temp_shelf'] for data in replica_data]

        # Calculate the subplot position
        row = i % 4   # Determine the row index
        col = i // 4  # Determine the column index

        # Create the first y-axis for scoring function values
        ax1 = axs[row, col]
        ax1.plot(accepted_structures, scoring_function_values, color=colors[i], label='Scoring Function')
        ax1.set_xlabel('Accepted Structures')
        ax1.set_ylabel('Scoring Function')
        ax1.set_title(f'Replica {replica_num}')
        ax1.legend()
#        ax1.set_ylim([-1, 2])  # Set the fixed y-axis limits for scoring function values

        # Create the second y-axis for temp shelf values
        ax2 = ax1.twinx()
        ax2.plot(accepted_structures, temp_shelf_values, color=colors[i], linestyle='dashed', label='Temp Shelf')
        ax2.set_ylabel('Temp Shelf')
        ax2.set_ylim([0, 320])  # Set the fixed y-axis limits for temp shelf values
        ax2.legend()

        # Adjust the spacing between subplots
        fig.subplots_adjust(hspace=0.5)

    # Remove unused subplots if there are fewer replicas
    for j in range(num_replicas, 12):
        row = j % 4   # Determine the row index
        col = j // 4  # Determine the column index
        fig.delaxes(axs[row, col])

    fig.suptitle(outname, fontsize=16, fontweight='bold')
#    fig.tight_layout(rect=[0, 0, 1, 0.95])

    
    plt.savefig(outname+'_re.png')





def plot_simulation_data(simulation_data, outname, scoring_f):
    # Extract values over accepted structures
    accepted_structures = range(len(simulation_data))
    mfe_e_values = [data['mfe_e'] for data in simulation_data]
    target_e_values = [data['target_e'] for data in simulation_data]
    d_mfe_target_values = [data['d_mfe_target'] for data in simulation_data]
    d_mfe_subopt_values = [data['d_mfe_subopt'] for data in simulation_data]
    precision_values = [data['precision'] for data in simulation_data]
    recall_values = [data['recall'] for data in simulation_data]
    mcc_values = [data['mcc'] for data in simulation_data]
    score_values = [data['score'] for data in simulation_data]
    reb_values = [data['reb'] for data in simulation_data]
    web_values = [data['web'] for data in simulation_data]
    reb_subopt_values = [data['reb_subopt'] for data in simulation_data]
    delta_mfe_subopt_norm_values = [data['d_mfe_subopt_norm'] for data in simulation_data]
    scoring_function_values = [data['scoring_function'] for data in simulation_data]
    
    # Define line colors
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'teal', 'magenta', 'cyan', 'gray', 'black', 'lime', 'pink', 'olive', 'navy']

    # Create subplots
    fig, axs = plt.subplots(4, 3, figsize=(24, 16))


    plt.subplots_adjust(top=0.9, bottom=0.05, left=0.1, right=0.95, hspace=0.4, wspace=0.3)  # Adjust spacing

    gs = GridSpec(4, 3, figure=fig, width_ratios=[1, 1, 1.3])

    fig.suptitle(outname, fontsize=16, fontweight='bold')  # Add title



    axs[0, 0].plot(accepted_structures, mfe_e_values, color=colors[0])
    axs[0, 0].set_xlabel('Accepted Structures')
    axs[0, 0].set_ylabel('MFE Energy')
    axs[0, 0].set_title('MFE Energy')

    axs[1, 0].plot(accepted_structures, target_e_values, color=colors[1])
    axs[1, 0].set_xlabel('Accepted Structures')
    axs[1, 0].set_ylabel('Target Energy')
    axs[1, 0].set_title('Target Energy')

    axs[2, 0].plot(accepted_structures, d_mfe_target_values, color=colors[2])
    axs[2, 0].set_xlabel('Accepted Structures')
    axs[2, 0].set_ylabel('Delta MFE to Target')
    axs[2, 0].set_title('Delta MFE to Target')

    axs[3, 0].plot(accepted_structures, d_mfe_subopt_values, color=colors[3])
    axs[3, 0].set_xlabel('Accepted Structures')
    axs[3, 0].set_ylabel('Delta MFE to Subopt')
    axs[3, 0].set_title('Delta MFE to Subopt')

    axs[0, 1].plot(accepted_structures, precision_values, color=colors[4])
    axs[0, 1].set_xlabel('Accepted Structures')
    axs[0, 1].set_ylabel('Precision')
    axs[0, 1].set_title('Precision')

    axs[1, 1].plot(accepted_structures, recall_values, color=colors[5])
    axs[1, 1].set_xlabel('Accepted Structures')
    axs[1, 1].set_ylabel('Recall')
    axs[1, 1].set_title('Recall')

    axs[2, 1].plot(accepted_structures, mcc_values, color=colors[6])
    axs[2, 1].set_xlabel('Accepted Structures')
    axs[2, 1].set_ylabel('MCC')
    axs[2, 1].set_title('MCC')

    axs[3, 1].plot(accepted_structures, score_values, color=colors[7])
    axs[3, 1].set_xlabel('Accepted Structures')
    axs[3, 1].set_ylabel('Score')
    axs[3, 1].set_title('Score')

    axs[0, 2].plot(accepted_structures, reb_values, color=colors[8])
    axs[0, 2].set_xlabel('Accepted Structures')
    axs[0, 2].set_ylabel('REB')
    axs[0, 2].set_title('REB')

    axs[1, 2].plot(accepted_structures, web_values, color=colors[9])
    axs[1, 2].set_xlabel('Accepted Structures')
    axs[1, 2].set_ylabel('WEB')
    axs[1, 2].set_title('WEB')

    axs[2, 2].plot(accepted_structures, reb_subopt_values, color=colors[10])
    axs[2, 2].set_xlabel('Accepted Structures')
    axs[2, 2].set_ylabel('REB Subopt')
    axs[2, 2].set_title('REB Subopt')
    
    axs[3, 2].plot(accepted_structures, scoring_function_values, color=colors[11])
    if scoring_f == 'dmt':
        axs[3, 2].set_title('Current Scoring Function: d_mfe_target')
    elif scoring_f == 'mcc':
        axs[3, 2].set_title('Current Scoring Function: MCC')
    elif scoring_f == 'mix':
        axs[3, 2].set_title('Current Scoring Function: Score')

#    axs[3, 2].set_title('Current Scoring Function')
    axs[3, 2].set_xlabel('Accepted Structures')
    axs[3, 2].set_ylabel('Current Scoring Function')



    # Save the plot to a file
    plt.savefig(outname + '.png')





class ScoreSeq:
    def __init__(self, sequence):
        self.sequence = sequence
#        self.mfe_ss = None
#        self.mfe_e = None
#        self.target_e = None
#        self.d_mfe_target = None
#        self.subopt_e = 0
#        self.d_mfe_subopt = 0

    def get_replica_num(self, rep_num):
        self.replica_num = rep_num

    def get_temp_shelf(self, temp):
        self.temp_shelf = temp

    def get_sim_step(self, step):
        self.sim_step = step
        
    def get_mfe_e(self, mfe):
        self.mfe_e = mfe

    def get_mfe_ss(self, ss):
        self.mfe_ss = ss

    def get_target_e(self, e_target):
        self.target_e = e_target

    def get_d_mfe_target(self, mfe, e_target):
        self.d_mfe_target = e_target -mfe
    
    def get_d_mfe_target_pk(self, mfe, e_target):
        self.d_mfe_target = abs(e_target -mfe)

    def get_subopt_e(self, e_subopt):
        self.subopt_e = e_subopt

    def get_subopt_ss(self, ss_subopt):
        self.subopt_ss = ss_subopt

    def get_d_mfe_subopt(self, mfe, e_subopt):
        self.d_mfe_subopt = e_subopt -mfe

    def get_oligomerization(self):
        mfe_oligo = energy_of_oligomer(self.sequence)
        self.oligomerization = if_oligomer(self.mfe_e, mfe_oligo)

    def get_precision(self, precision):
        self.precision = 1-precision
        
    def get_recall(self, recall):
        self.recall = 1-recall
    
    def get_mcc(self, mcc):
        self.mcc = 1-mcc


    def get_deltafm(self, delta):
        self.deltaF_m = delta 

    def get_reb(self):
        self.reb = -(self.mfe_e + 0.31*len(self.sequence) - 7.05)/(-0.049*len(self.sequence) - 20.743)

    def get_web(self):
        self.web = (self.mcc**4)*self.reb

    def get_reb_subopt(self):
        self.reb_subopt = -(self.subopt_e + 0.31*len(self.sequence) - 7.05)/(-0.049*len(self.sequence) - 20.743)

    def get_d_mfe_subopt_norm(self):
        self.d_mfe_subopt_norm = (self.d_mfe_subopt/3)

#    def get_web_subopt(self):
#        self.web_subopt = (self.mcc**4)*self.reb_subopt*0.25
    

    def get_score(self):
        self.score = self.d_mfe_target +self.mcc + self.web
#        self.score = -self.mcc
#        self.score = self.d_mfe_target -self.mcc + self.web
#        self.score = self.d_mfe_target

    def get_scoring_function(self, scoring_f):
        if scoring_f == 'dmt':
            self.scoring_function = self.d_mfe_target
        elif scoring_f == 'mcc':
            self.scoring_function = self.mcc
        elif scoring_f == 'mfe':
            self.scoring_function = self.mfe_e
        elif scoring_f == 'mix':
            self.scoring_function = self.d_mfe_target +self.mcc + self.web
        elif scoring_f == 'mix2':
            self.scoring_function = self.d_mfe_target + 8*self.mcc 
            
def get_first_suboptimal_structure_and_energy(sequence):
    """
    Accept a single RNA sequence and report its first possible secondary structure and energy.
    """
    RNA.cvar.uniq_ML = 1
    subopt_data = {'sequence': sequence}
    subopt_list = []

    def print_subopt_result(structure, energy, data):
            if structure != None:
                subopt_list.append([structure, energy])

    a = RNA.fold_compound(sequence)

    for i in range(100,5000,100):
        a.subopt_cb(i, print_subopt_result, subopt_data)
        if len(subopt_list) >=2:
            break
        subopt_list = []        
        
    subopt_list.sort(key=lambda x: x[1])
    
    if len(subopt_list) ==0:
        structure1 = "."*len(sequence)
        energy1 = 0
    else:
        structure1 = subopt_list[1][0]
        energy1 = subopt_list[1][1]

    return structure1, energy1



def di_nuc_analyser_penalty_aplier(input_csv_file):
    """
    Accepts an input CSV file, calculates the number of di_nucleotides, applies penalties,
    and saves the updated results in the input CSV file.

    The input CSV file should have 3 columns without a header. The first column is considered
    as the designed sequences, and the other two columns are scores.

    Args:
        input_csv_file (str): Path to the input CSV file.

    Returns:
        None
    """
    di_nuc = DiNucCounter(input_csv_file)
    di_nuc.add_di_nuc_df()
    di_nuc.di_nuc_calculator()
    di_nuc.df_to_file()
    DiNuc_penalty_applier(input_csv_file)

def dimer_multi_chain_energy_structure_flags(input_csv_file_sequences, ouput_csv_file, target_structure):
    """
    Calculates energy values and structural flags for monomers, dimers, and multi-chain systems based on input sequences.

    Args:
        input_csv_file_sequences (str): Path to the input CSV file containing sequences, the sequences must be under a column named 'sequence'.
        ouput_csv_file (str): Path to the output CSV file to save the results.
        target_structure (str): Target structure for energy calculations.

    Returns:
        None
    """
    df = pd.read_csv(input_csv_file_sequences, delimiter=",")
    chain_flag = True
    monomer_energy = []
    dimer_energy = []
    dimer_structure = []
    dimerization_list = []
    Multichain_list = []
    seq1_seq2_energy = []
    seq1_seq2_struct = []
    seq1_energy = []
    seq1_struct = []
    seq2_energy = []
    seq2_struct = []
    multi_chain_flag = []
    multi_chain_eval = []
    multi_chain_delta = []
    seq11_energy_l = []
    seq11_struct = []
    seq22_energy_l = []
    seq22_struct = []
    for seq1 in df["sequence"]:
        if "&" not in seq1:
            chain_flag = False
            seq2 = seq1
            dme = DimerMonomerEnergy(seq1, seq2, target_structure)
            mono_struct, mono = dme.monomer_energy(seq1)
            dimer_struct, dimer = dme.dimer_energy()
            dimerization = dme.dimerization_flag()
            monomer_energy.append(round(mono, 3))
            dimer_energy.append(round(dimer, 3))
            dimer_structure.append(dimer_struct)
            dimerization_list.append(dimerization)
        else:
            seq2 = seq1.split("&")[1]
            seq1 = seq1.split("&")[0]
            dme = DimerMonomerEnergy(seq1, seq2, target_structure)
            seq_structure, seq_energy = dme.monomer_energy(seq1)
            seq1_energy.append(round(seq_energy, 3))
            seq1_struct.append(seq_structure)
            seq_structure, seq_energy = dme.monomer_energy(seq2)
            seq2_energy.append(round(seq_energy, 3))
            seq2_struct.append(seq_structure)
            dimer_struct, dimer = dme.dimer_energy()
            seq1_seq2_energy.append(round(dimer, 3))
            seq1_seq2_struct.append(dimer_struct)
            multichain = dme.multi_chain_flag()
            multi_chain_flag.append(multichain)
            multi_chain_eval_energy = dme.multi_chain_eval()
            multi_chain_eval.append(multi_chain_eval_energy)
            multi_chain_delta.append(round((dimer - multi_chain_eval_energy), 3))
            dme = DimerMonomerEnergy(seq1, seq1, target_structure)
    if chain_flag:
        df["seq1_energy"] = seq1_energy
        df["seq2_energy"] = seq2_energy
        df["seq1_seq2_energy"] = seq1_seq2_energy
        df["multi_chain_eval"] = multi_chain_eval
        df["delta"] = multi_chain_delta
        df["seq1_struct"] = seq1_struct
        df["seq2_struct"] = seq2_struct
        df["seq1_seq2_struct"] = seq1_seq2_struct
        df["multi_chain_flag"] = multi_chain_flag
    else:
        df["Monomer_energy"] = monomer_energy
        df["Dimer_energy"] = dimer_energy
        df["Dimer_structure"] = dimer_structure
        df["Dimerization"] = dimerization_list
    df.to_csv(ouput_csv_file, index=False)

def check_string(string):
    """
    Check if the given sequence is multichain sequence or a single one.

    Args:
        string (str): The input sequence to check.

    Returns:
        tuple: A tuple containing a boolean value indicating if multiple sequences are present
               and a list of locations where the '&' character is found.

    """
    if "&" in string:
        locations = [index for index, char in enumerate(string) if char == "&"]
        return True, locations
    return False, []


def add_ampersand(string, indices):
    """
    Add ampersands to the exact locations of the predicted secondary structure.

    Args:
        string (str): The input string to modify.
        indices (list): A list of indices indicating the positions to add ampersands.

    Returns:
        str: The modified string with ampersands added.

    """
    for index in indices:
        if index < len(string):
            string = string[:index] + "&" + string[index:]
    return string


def multi_chain_energy_reporter(multi_chain_sequence):
    """
    Accept either single or multi-chain sequences and return the predicted structure and energy.

    Args:
        multi_chain_sequence (str): The input sequence to analyze.

    Returns:
        tuple: A tuple containing the predicted structure and minimum free energy.

    """
    fc = RNA.fold_compound(multi_chain_sequence)
    structure, mfe = fc.mfe_dimer()
    has_ampersand, ampersand_locations = check_string(multi_chain_sequence)
    structure = add_ampersand(structure, ampersand_locations)
    return structure, mfe


def main():
    """Main function for testing each individual function separately."""
    # sequence = "CCGCCuaACACUGCCAAUGCCGGUCCCAAGCCCGGAU&AAAAGUGGAGGGGGCGG&CCGCCuaACACUGCCAAUGCCGGUCCCAAGCCCGGAUAAAAGUGGAGGGGGCGG"
    # result = multi_chain_energy_reporter(sequence)
    # print(result)
    #dimer_multi_chain_energy_structure_flags("init.csv", "tmp_init.csv", "..(((((......)))))..........")
    #di_nuc_analyser_penalty_aplier("init.csv")

    sequence = "CCGCCuaACACUGCCAAUGCCGGUCCCAAGCCCGGAUAAAAGUGGAGGGGGCGG"
    subopt_structure1,subopt_energy1 = get_first_suboptimal_structure_and_energy(sequence)
if __name__ == "__main__":
    main()
