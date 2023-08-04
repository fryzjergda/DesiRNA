#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import RNA
import sys
import math
import pandas as pd
from utils.pcofold_dimer_multichain_energy import *
from utils.di_nuc_analyser import *
import RNA


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


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

def plot_simulation_data_combined(simulation_data, outname, scoring_f, alt):
    # Define the metrics to plot
    if alt == None:
        metrics = ['mfe_e', 'precision', 'd_mfe_subopt', 'target_e', 'recall','temp_shelf', 'd_mfe_target', 'mcc','scoring_function']
    else:
        metrics = ['mfe_e', 'precision', 'd_alt_mfe_target', 'target_e', 'recall','temp_shelf', 'd_mfe_target', 'mcc','scoring_function']

    # Create subplots
    num_plots = len(metrics)
    num_cols = 3
    num_rows = (num_plots + num_cols - 1) // num_cols
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(24, 16))

    # Flatten the axs if there's only one row
    if num_rows == 1:
        axs = [axs]

    scoring_function_titles = {
    'dmt': 'Scoring Function - dE MFE target',
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


        ax.set_xlabel('RE steps')
        ax.set_ylabel(metric.capitalize().replace('_', ' '))
        ax.set_title(metric.capitalize().replace('_', ' '))
        ax.legend(loc='upper right')


        if metric == 'mcc':
            ax.set_title('1 - MCC')
        if metric == 'precision':
            ax.set_title('1 - Precision')
        if metric == 'recall':
            ax.set_title('1 - Recall')
        if metric == 'mfe_e':
            ax.set_title('MFE E')
        if metric == 'temp_shelf':
            ax.set_title('Temperature shelf')
        if metric == 'target_e':
            ax.set_title('Target E')
        if metric == 'd_mfe_target':
            ax.set_title('dE MFE target')
        if (metric == 'd_mfe_subopt') and (alt == None):
            ax.set_title('dE MFE subopt')
        if (metric == 'd_alt_mfe_target') and (alt != None):
            ax.set_title('dE altMFE target')


        if metric == 'scoring_function':
            if scoring_f in scoring_function_titles:
                ax.set_title(scoring_function_titles[scoring_f])

    if num_plots % num_cols != 0:
        axs[-1, -1].axis('off')

    # Add overall title to the plot
    fig.suptitle(outname, fontsize=16, fontweight='bold')

    # Adjust the layout
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save the plot to a file
    plt.savefig(outname + '_replicas.png')





def _plot_simulation_data_combined(simulation_data, outname, scoring_f):
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


class Stats:
    def __init__(self):
        self.global_step = 0
        self.step = 0
        self.acc_mc_step = 0
        self.acc_mc_better_e = 0
        self.rej_mc_step = 0
        self.acc_re_step = 0
        self.acc_re_better_e = 0
        self.rej_re_step = 0
        
    def update_global_step(self):
        self.global_step += 1
    
    def update_step(self, steps):
        self.step += steps
    
    def update_acc_mc_step(self):
        self.acc_mc_step += 1

    def update_acc_mc_better_e(self):
        self.acc_mc_better_e += 1

    def update_rej_mc_step(self):
        self.rej_mc_step += 1

    def reset_mc_stats(self):
        self.acc_mc_step = 0
        self.acc_mc_better_e = 0
        self.rej_mc_step = 0

    def update_acc_re_step(self):
        self.acc_re_step += 1

    def update_acc_re_better_e(self):
        self.acc_re_better_e += 1

    def update_rej_re_step(self):
        self.rej_re_step += 1
        
    def get_tot_mc_steps(self):
        self.tot_mc_steps = self.acc_mc_step + self.rej_mc_step
    
    def get_mc_acc_ratio(self):
        self.mc_acc_ratio = round(self.acc_mc_step/self.tot_mc_steps, 3)
    
    def get_mc_rand_acc_ratio(self):
        self.mc_rand_acc_ratio = round(self.acc_mc_better_e/self.acc_mc_step, 3)
    
    
    def get_tot_re_steps(self):
        self.tot_re_steps = self.acc_re_step + self.rej_re_step
 
    def get_re_acc_ratio(self):
        self.re_acc_ratio = round(self.acc_re_step/self.tot_re_steps, 3)

    def get_re_rand_acc_ratio(self):
        self.re_rand_acc_ratio = round(self.acc_re_better_e/self.acc_re_step, 3)

    


class ScoreSeq:
    def __init__(self, sequence):
        self.sequence = sequence
#        self.mfe_ss = None
#        self.mfe_e = None
#        self.target_e = None
#        self.d_mfe_target = None
#        self.subopt_e = 0
#        self.d_mfe_subopt = 0
        self.reb = 0
        self.web = 0
        self.reb_subopt = 0
        self.web_subopt = 0
        self.d_mfe_subopt_norm =0
        self.score = 0
        self.recall = 0
        self.precision = 0
        self.d_mfe_subopt = 0

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
    def get_alt_mfe_e(self, altmfe):
        self.alt_mfe_e = altmfe
    
    def get_alt_target_e(self, e_alt_target):
        self.alt_target_e = e_alt_target

    
    def get_d_alt_mfe_target(self, alt_mfe, e_alt_target):
        self.d_alt_mfe_target = e_alt_target -alt_mfe

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
            self.scoring_function = self.d_mfe_target + 2*self.mcc 
        elif scoring_f == 'alt':
            self.scoring_function = self.d_mfe_target + self.d_alt_mfe_target

    def get_scoring_function_w_subopt(self):
        self.scoring_function = self.d_mfe_target -self.d_mfe_subopt 
    
    def get_scoring_function_w_oligo(self):
        if self.oligomerization == True:
            self.scoring_function = 300*self.d_mfe_target
        else:
            self.scoring_function = self.d_mfe_target
            
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
