#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This Python module is designed for the analysis and visualization of RNA sequence simulation data.
It includes functionalities for plotting simulation metrics, calculating various statistical
properties of RNA sequences, and assessing their structural and energetic characteristics.

Key components of this module include:
- Functions for generating color-coded plots of simulation data across different replicas.
- A 'Stats' class for tracking statistical data during simulation processes.
- A 'ScoreSeq' class for representing and scoring RNA sequences based on various criteria.
- Utility functions for calculating suboptimal RNA structures and their energies.

The module makes use of the RNA package for RNA secondary structure prediction and energy calculations,
and matplotlib for generating plots. It's tailored for use in RNA sequence analysis and simulation environments.

Dependencies:
- matplotlib: For plotting simulation data.
- RNA: For RNA secondary structure prediction and energy calculations.
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
import RNA
from utils.pcofold_dimer_multichain_energy import *

matplotlib.use('Agg')


def get_replica_color(replica_num, num_replicas):
    """
    Determine the color for a given replica based on its number.

    Parameters:
    replica_num (int): The number of the replica.
    num_replicas (int): The total number of replicas.

    Returns:
    str: A hexadecimal color code corresponding to the replica number.
    """
    colors = ['#0000FF', '#00BFFF', '#00FF00', '#7CFC00', '#FFFF00',
              '#FFD700', '#FFA500', '#FF0000', '#8B0000', '#B22222']

    # Calculate the index of the color for the given replica
    if num_replicas == 1:
        color_index = 0
    else:
        color_index = int((replica_num - 1) * (len(colors) - 1) / (num_replicas - 1))

    return colors[color_index]


def plot_simulation_data_combined(simulation_data, alt, sim_options):
    """
    Plot combined simulation data for various metrics across different replicas.

    Parameters:
    simulation_data (list): A list of dictionaries containing simulation data.
    outname (str): The name of the output file for the plot.
    alt (bool): Flag to determine alternate plotting metrics.
    infile (str): The name of the input file used for simulation data.

    This function creates a plot for each metric in the simulation data, showing
    the performance of different replicas. The plot is saved to a file.
    """
    metrics = ['scoring_function', 'edesired_minus_Epf', 'Epf', 'edesired', 'mcc']

    num_plots = len(metrics)
    num_cols = 1
    num_rows = num_plots + 1  # + num_cols - 1) // num_cols

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
        row = i  # // num_cols
        ax = axs[row]  # [col]
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

        ax.yaxis.label.set_position(ylabel_pos)

        if metric == 'mcc':
            ax.set_ylabel('1 - MCC', fontsize=12)
            ax.set_xlabel('Simulation steps', fontsize=12)
        if metric == 'precision':
            ax.set_title('1 - Precision')
        if metric == 'recall':
            ax.set_title('1 - Recall')
        if metric == 'Epf':
            ax.set_ylabel('Epf', fontsize=12)
        if metric == 'temp_shelf':
            ax.set_title('Temperature shelf')
        if metric == 'edesired':
            ax.set_ylabel('E desired', fontsize=12)
        if metric == 'edesired_minus_Epf':
            ax.set_ylabel('E desired - Epf', fontsize=12)
        if (metric == 'esubopt_minus_Epf') and (alt == None):
            ax.set_title('dE Epf subopt')
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
    fig.suptitle(sim_options.infile.split(".")[0], fontsize=16, fontweight='bold')

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

    # Adjust the layout to make room for the legend
    fig.tight_layout(rect=[0, 0.0, 1, 0.95])

    # Save the plot to a file
    plt.savefig(sim_options.outname + '_replicas.png', dpi=400)


class Stats:
    """
    A class used to track various statistical metrics during a simulation process.

    Attributes:
    - global_step (int): The total number of steps processed in the simulation.
    - step (int): A counter for the number of steps.
    - acc_mc_step (int): The number of accepted Monte Carlo steps.
    - acc_mc_better_e (int): The count of accepted Monte Carlo steps with better energy.
    - rej_mc_step (int): The number of rejected Monte Carlo steps.
    - acc_re_step (int): The number of accepted Replica Exchange steps.
    - acc_re_better_e (int): The count of accepted Replica Exchange steps with better energy.
    - rej_re_step (int): The number of rejected Replica Exchange steps.
    """

    def __init__(self):
        """
        Initialize the Stats instance with all statistical counters set to zero.
        """
        self.global_step = 0
        self.step = 0
        self.acc_mc_step = 0
        self.acc_mc_better_e = 0
        self.rej_mc_step = 0
        self.acc_re_step = 0
        self.acc_re_better_e = 0
        self.rej_re_step = 0

    def update_global_step(self):
        """
        Increment the global step counter by one.
        """
        self.global_step += 1

    def update_step(self, steps):
        """
        Increment the step counter by a specified number.

        Parameters:
        steps (int): The number of steps to increment.
        """
        self.step += steps

    def update_acc_mc_step(self):
        """
        Increment the counter for accepted Monte Carlo steps by one.
        """
        self.acc_mc_step += 1

    def update_acc_mc_better_e(self):
        """
        Increment the counter for accepted Monte Carlo steps with better energy.
        """
        self.acc_mc_better_e += 1

    def update_rej_mc_step(self):
        """
        Increment the counter for rejected Monte Carlo steps by one.
        """
        self.rej_mc_step += 1

    def reset_mc_stats(self):
        """
        Reset the Monte Carlo statistical counters to zero.
        """
        self.acc_mc_step = 0
        self.acc_mc_better_e = 0
        self.rej_mc_step = 0

    def update_acc_re_step(self):
        """
        Increment the counter for accepted Replica Exchange steps by one.
        """
        self.acc_re_step += 1

    def update_acc_re_better_e(self):
        """
        Increment the counter for accepted Replica Exchange steps with better energy.
        """
        self.acc_re_better_e += 1

    def update_rej_re_step(self):
        """
        Increment the counter for rejected Replica Exchange steps by one.
        """
        self.rej_re_step += 1

    def get_tot_mc_steps(self):
        """
        Calculate and store the total number of Monte Carlo steps.
        """
        self.tot_mc_steps = self.acc_mc_step + self.rej_mc_step

    def get_mc_acc_ratio(self):
        """
        Calculate and store the acceptance ratio for Monte Carlo steps.
        """
        self.mc_acc_ratio = round(self.acc_mc_step / self.tot_mc_steps, 3)

    def get_mc_rand_acc_ratio(self):
        """
        Calculate and store the random acceptance ratio for Monte Carlo steps.
        """
        self.mc_rand_acc_ratio = round(self.acc_mc_better_e / self.acc_mc_step, 3)

    def get_tot_re_steps(self):
        """
        Calculate and store the total number of Replica Exchange steps.
        """
        self.tot_re_steps = self.acc_re_step + self.rej_re_step

    def get_re_acc_ratio(self):
        """
        Calculate and store the acceptance ratio for Replica Exchange steps.
        """
        self.re_acc_ratio = round(self.acc_re_step / self.tot_re_steps, 3)

    def get_re_rand_acc_ratio(self):
        """
        Calculate and store the random acceptance ratio for Replica Exchange steps.
        """
        self.re_rand_acc_ratio = round(self.acc_re_better_e / self.acc_re_step, 3)


class ScoreSeq:
    """
    A class for representing and scoring an RNA sequence. It tracks various metrics related to
    the sequence, such as energy scores, structural properties, and statistical measures.
    """

    def __init__(self, sequence):
        """
        Initialize a ScoreSeq instance with a specific RNA sequence and default values for various
        scoring parameters.

        Parameters:
        sequence (str): The RNA sequence to be scored and analyzed.
        """
        self.sequence = sequence
        self.scoring_function = 0
        self.replica_num = None
        self.temp_shelf = None
        self.sim_step = 0
        self.edesired_minus_Epf = 0
        self.Epf = 0
        self.edesired = 0
        self.mcc = 0
        self.mfe_ss = None
        self.subopt_e = 0
        self.esubopt_minus_Epf = 0
        self.sln_Epf = 0
#        self.score = 0
        self.MFE = 0
        self.edesired_minus_MFE = 0
        self.recall = 0
        self.precision = 0
        self.edesired2 = 0
        self.edesired2_minus_Epf = 0

    def get_replica_num(self, rep_num):
        """
        Set the replica number for this sequence.

        Parameters:
        rep_num (int): The replica number to be assigned.
        """
        self.replica_num = rep_num

    def get_temp_shelf(self, temp):
        """
        Set the temperature shelf value for this sequence.

        Parameters:
        temp (float): The temperature shelf value to be assigned.
        """
        self.temp_shelf = temp

    def get_sim_step(self, step):
        """
        Set the simulation step number for this sequence.

        Parameters:
        step (int): The simulation step number to be assigned.
        """
        self.sim_step = step

    def get_Epf(self, Epf):
        """
        Set the energy of the folded RNA structure (Epf, energy of partition function) for this sequence.

        Parameters:
        Epf (float): The Epf value to be assigned.
        """
        self.Epf = Epf

    def get_mfe_ss(self, ss):
        """
        Set the minimum free energy secondary structure (mfe_ss) for this sequence.

        Parameters:
        ss (str): The mfe_ss value to be assigned.
        """
        self.mfe_ss = ss

    def get_edesired(self, e_target):
        """
        Set the desired energy (edesired) for this sequence.

        Parameters:
        e_target (float): The edesired value to be assigned.
        """
        self.edesired = e_target

    def get_edesired_minus_Epf(self, Epf, e_target):
        """
        Calculate and set the difference between the desired energy and Epf for this sequence.

        Parameters:
        Epf (float): The Epf value.
        e_target (float): The desired energy value.
        """
        self.edesired_minus_Epf = e_target - Epf

    def get_edesired2(self, e_target):
        """
        Set an alternative desired energy (edesired2) for this sequence.

        Parameters:
        e_target (float): The edesired2 value to be assigned.
        """
        self.edesired2 = e_target

    def get_edesired2_minus_Epf(self, Epf, e_target):
        """
        Calculate and set the difference between the alternative desired energy and Epf for this sequence.

        Parameters:
        Epf (float): The Epf value.
        e_target (float): The alternative desired energy value.
        """
        self.edesired2_minus_Epf = e_target - Epf

    def get_subopt_e(self, e_subopt):
        """
        Set the suboptimal energy (subopt_e) for this sequence.

        Parameters:
        e_subopt (float): The suboptimal energy value to be assigned.
        """
        self.subopt_e = e_subopt

    def get_esubopt_minus_Epf(self, Epf, e_subopt):
        """
        Calculate and set the difference between the suboptimal energy and Epf for this sequence.

        Parameters:
        Epf (float): The Epf value.
        e_subopt (float): The suboptimal energy value.
        """
        self.esubopt_minus_Epf = e_subopt - Epf

    def get_precision(self, precision):
        """
        Set the precision metric for this sequence.

        Parameters:
        precision (float): The precision value to be assigned.
        """
        self.precision = 1 - precision

    def get_recall(self, recall):
        """
        Set the recall metric for this sequence.

        Parameters:
        recall (float): The recall value to be assigned.
        """
        self.recall = 1 - recall

    def get_mcc(self, mcc):
        """
        Set the Matthews correlation coefficient (MCC) for this sequence.

        Parameters:
        mcc (float): The MCC value to be assigned.
        """
        self.mcc = 1 - mcc

    def get_sln_Epf(self):
        """
        Calculate and set the Epf scaled by length value for this sequence.
        """
        self.sln_Epf = (self.Epf + 0.3759 * len(self.sequence) + 5.7534) / 10

    def get_MFE(self):
        """
        Calculate and set the minimum free energy (MFE) for the RNA folding of this sequence.
        """
        self.MFE = RNA.fold(self.sequence)[1]

    def get_edesired_minus_MFE(self):
        """
        Calculate and set the difference between the desired energy and MFE for this sequence.
        """
        self.edesired_minus_MFE = self.edesired - self.MFE

    def get_ensemble_defect(self, sec_struct):
        """
        Calculate and set the ensemble defect for a given secondary structure of this sequence.

        Parameters:
        sec_struct (str): The secondary structure for which to calculate the ensemble defect.
        """
        md = RNA.md()
        fc = RNA.fold_compound(self.sequence, md)
        (_, mfe) = fc.mfe()
        fc.exp_params_rescale(mfe)
        (_, pf) = fc.pf()
        self.ensemble_defect = fc.ensemble_defect(sec_struct)

    def get_scoring_function(self, scoring_f):
        """
        Calculate and set the scoring function for this sequence based on various metrics.

        Parameters:
        scoring_f (list of tuples): A list of tuples, where each tuple contains a scoring function name and its weight.
        """
        self.scoring_function = 0
        for function, weight in scoring_f:
            if function == 'Ed-Epf':
                self.scoring_function += self.edesired_minus_Epf * weight
            elif function == '1-MCC':
                self.scoring_function += self.mcc * 10 * weight
            elif function == 'sln_Epf':
                self.scoring_function += self.sln_Epf * weight
            elif function == 'Ed-MFE':
                self.scoring_function += self.edesired_minus_MFE * weight
            elif function == '1-precision':
                self.scoring_function += self.precision * 10 * weight
            elif function == '1-recall':
                self.scoring_function += self.recall * 10 * weight
            elif function == 'Edef':
                self.scoring_function += self.ensemble_defect * weight

    def get_scoring_function_w_alt_ss(self):
        """
        Update the scoring function for this sequence by including the effect of an alternative secondary structure.
        """
        self.scoring_function = self.scoring_function + self.edesired2_minus_Epf

    def get_scoring_function_w_subopt(self):
        """
        Update the scoring function for this sequence by including the effect of a suboptimal energy.
        """
        self.scoring_function = self.scoring_function - self.esubopt_minus_Epf

    def get_scoring_function_monomer(self):
        """
        Update the scoring function for this sequence by including the effect of monomer fraction.
        """
        fc = RNA.fold_compound(self.sequence + "&" + self.sequence)
        self.oligo_fraction = oligo_fraction(self.sequence + "&" + self.sequence, fc)
        self.monomer_bonus = kTlog_monomer_fraction(self.oligo_fraction)
        self.scoring_function = self.scoring_function + self.monomer_bonus

    def get_scoring_function_oligomer(self, fc):
        """
        Update the scoring function for this sequence by including the effect of oligomer fraction.

        Parameters:
        fc (RNA.fold_compound): The fold compound object for RNA secondary structure.
        """
        self.oligo_fraction = oligo_fraction(self.sequence, fc)
        self.oligomer_bonus = kTlog_oligo_fraction(self.oligo_fraction)
        self.scoring_function = self.scoring_function + self.oligomer_bonus

    def update_scoring_function_w_motifs(self, motif_bonus):
        """
        Update the scoring function for this sequence by adding a motif bonus.

        Parameters:
        motif_bonus (float): The motif bonus value to be added to the scoring function.
        """
        self.scoring_function += motif_bonus


def get_first_suboptimal_structure_and_energy(sequence, a):
    """
    Find the first suboptimal secondary structure and its energy for a given RNA sequence.

    Parameters:
    sequence (str): The RNA sequence.
    a (RNA.fold_compound): The RNA fold compound object for the given sequence.

    Returns:
    tuple: A tuple containing the first suboptimal structure and its corresponding energy.
    """
    RNA.cvar.uniq_ML = 1
    subopt_data = {'sequence': sequence}
    subopt_list = []

    def print_subopt_result(structure, energy, data):
        if structure != None:
            subopt_list.append([structure, energy])

    for i in range(100, 5000, 100):
        a.subopt_cb(i, print_subopt_result, subopt_data)
        if len(subopt_list) >= 2:
            break
        subopt_list = []

    subopt_list.sort(key=lambda x: x[1])

    if len(subopt_list) == 0:
        structure1 = "." * len(sequence)
        energy1 = 0
    else:
        structure1 = subopt_list[1][0]
        energy1 = subopt_list[1][1]

    return structure1, energy1
