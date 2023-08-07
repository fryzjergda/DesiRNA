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
from matplotlib import ticker

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

def plot_simulation_data_combined(simulation_data, outname, alt, infile):
    # Define the metrics to plot
#    if alt == None:
    metrics = ['scoring_function', 'edesired_minus_mfe','mfe', 'edesired', 'mcc']
    #else:
     #   metrics = ['mfe', 'precision', 'd_alt_mfe_target', 'edesired', 'recall','temp_shelf', 'edesired_minus_mfe', 'mcc','scoring_function']
    

    # Create subplots
    num_plots = len(metrics)
    num_cols = 1
    num_rows = num_plots+1 #+ num_cols - 1) // num_cols

    width_mm = 160  
    height_mm = (20 / 8) * width_mm  

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(width_mm / 25.4, height_mm / 25.4))

    # Flatten the axs if there's only one row
    if num_rows == 1:
        axs = [axs]


    # Get the number of replicas
    num_replicas = len(set(data['replica_num'] for data in simulation_data))
    
    handles, labels = [], []

    # Iterate over replicas and plot data
    for i, metric in enumerate(metrics):
        row = i #// num_cols
#        col = i % num_cols
        ax = axs[row]#[col]
        ax.tick_params(axis='both', labelsize=7)

        
        for replica_num in range(1, num_replicas + 1):
            replica_data = [data[metric] for data in simulation_data if data['replica_num'] == replica_num]
            accepted_structures = range(len(replica_data))
            color = get_replica_color(replica_num, num_replicas)  # Get color for the replica
            ax.plot(accepted_structures, replica_data, color=color, label=f"Replica {replica_num}")
            line, = ax.plot(accepted_structures, replica_data, color=color, label=f"Replica {replica_num}")
            if i == 0:  # Only collect handles and labels from the first subplot
                handles.append(line)
                labels.append(f"Replica {replica_num}")
                ylabel_pos = ax.yaxis.label.get_position()

#        ax.set_xlabel('Simulatiion steps')
#        ax.set_ylabel(metric.capitalize().replace('_', ' '))
#        ax.set_title(metric.capitalize().replace('_', ' '))
#        ax.legend(loc='upper right')
        ax.yaxis.label.set_position(ylabel_pos)

        if metric == 'mcc':
            ax.set_ylabel('1 - MCC', fontsize=12)
            ax.set_xlabel('Simulatiion steps', fontsize=12)
        if metric == 'precision':
            ax.set_title('1 - Precision')
        if metric == 'recall':
            ax.set_title('1 - Recall')
        if metric == 'mfe':
            ax.set_ylabel('MFE', fontsize=12)
        if metric == 'temp_shelf':
            ax.set_title('Temperature shelf')
        if metric == 'edesired':
            ax.set_ylabel('E desired', fontsize=12)
        if metric == 'edesired_minus_mfe':
            ax.set_ylabel('E desired - MFE', fontsize=12)
        if (metric == 'esubopt_minus_mfe') and (alt == None):
            ax.set_title('dE MFE subopt')
        if (metric == 'd_alt_mfe_target') and (alt != None):
            ax.set_title('dE altMFE target')


        if metric == 'scoring_function':
            ax.set_ylabel('Scoring function', fontsize=12)


    if num_plots % num_cols != 0:
        axs[-1, -1].axis('off')


    max_width = 0
    for ax in axs:
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.draw()  # Update the figure to get the correct tick label sizes
        tick_width = max([t.get_window_extent().width for t in ax.yaxis.get_ticklabels()])
        max_width = max(max_width, tick_width)

# Set the position of the y-labels based on the maximum width of the y-tick labels
    for ax in axs:
        ax.yaxis.labelpad = max_width   # Adjust the padding as needed


    # Add overall title to the plot
    fig.suptitle(infile.split(".")[0], fontsize=16, fontweight='bold')


    num_replica_columns = min(num_replicas, 4)

    legend_ax = axs[-1]
    legend_ax.legend(handles, labels, loc='center', ncol=num_replica_columns, fontsize=8, frameon=False)
    # Get the position of the last subplot
    last_subplot_pos = axs[-2].get_position()
    # Calculate the position for the legend based on the last subplot
    legend_pos = [last_subplot_pos.x0, last_subplot_pos.y0 - 0.15, last_subplot_pos.width, 0.1]  # Adjust as needed

    # Set the position for the legend
    legend_ax.set_position(legend_pos)
    legend_ax.axis('off')
    
#    plt.subplots_adjust(bottom=0.25)

    # Adjust the layout to make room for the legend
    fig.tight_layout(rect=[0, 0.0, 1, 0.95])
#    plt.subplots_adjust(bottom=0.1) 
#    plt.subplots_adjust(hspace=0.3, bottom=0.15)

    # Save the plot to a file
    plt.savefig(outname + '_replicas.png', dpi=400)





def _plot_simulation_data_combined(simulation_data, outname, scoring_f):
    # Define the metrics to plot
    metrics = ['mfe', 'precision', 'sln_mfe', 'edesired', 'recall', 'web',
               'edesired_minus_mfe', 'mcc', 'temp_shelf', 'esubopt_minus_mfe',
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
    'dmt': 'Scoring Function - edesired_minus_mfe',
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
        self.scoring_function = 0
        self.replica_num = None
        self.temp_shelf = None
        self.sim_step = 0
        self.edesired_minus_mfe =0
        self.mfe =0
        self.edesired = 0
        self.mcc = 0
        self.mfe_ss = None
        self.subopt_e = 0
        self.esubopt_minus_mfe = 0
        self.sln_mfe = 0
#        self.score = 0
        self.recall = 0
        self.precision = 0
        self.edesired2 = 0
        self.edesired2_minus_mfe =0
        self.dimer_mfe = 0    

    def get_replica_num(self, rep_num):
        self.replica_num = rep_num

    def get_temp_shelf(self, temp):
        self.temp_shelf = temp

    def get_sim_step(self, step):
        self.sim_step = step
        
    def get_mfe(self, mfe):
        self.mfe = mfe

    def get_mfe_ss(self, ss):
        self.mfe_ss = ss

    def get_edesired(self, e_target):
        self.edesired = e_target

    def get_edesired_minus_mfe(self, mfe, e_target):
        self.edesired_minus_mfe = e_target -mfe

    def get_edesired2(self, e_target):
        self.edesired2 = e_target

    def get_edesired2_minus_mfe(self, mfe, e_target):
        self.edesired2_minus_mfe = e_target -mfe
        
    def get_subopt_e(self, e_subopt):
        self.subopt_e = e_subopt

    def get_esubopt_minus_mfe(self, mfe, e_subopt):
        self.esubopt_minus_mfe = e_subopt -mfe

    def get_oligomerization(self):
        mfe_oligo = energy_of_oligomer(self.sequence)
        self.oligomerization = if_oligomer(self.mfe, mfe_oligo)

    def get_precision(self, precision):
        self.precision = 1-precision
        
    def get_recall(self, recall):
        self.recall = 1-recall
    
    def get_mcc(self, mcc):
        self.mcc = 1-mcc

    def get_sln_mfe(self):
        self.sln_mfe = -(self.mfe + 0.31*len(self.sequence) - 7.05)/(-0.049*len(self.sequence) - 20.743)
    

    def get_scoring_function(self, scoring_f):
        self.scoring_function = 0
        if 'ed-mfe' in scoring_f:
            self.scoring_function += self.edesired_minus_mfe
        if '1-mcc' in scoring_f:
            self.scoring_function += self.mcc*10
        if 'sln_mfe' in scoring_f:
            self.scoring_function += self.sln_mfe

    def get_scoring_function_w_alt_ss(self):        
        self.scoring_function = self.scoring_function + self.edesired2_minus_mfe
        
    def get_scoring_function_w_subopt(self):
        self.scoring_function = self.scoring_function -self.esubopt_minus_mfe 
    
    def get_scoring_function_w_oligo(self):
        if self.oligomerization == True:
            self.scoring_function = self.scoring_function + 300*self.edesired_minus_mfe

            
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
