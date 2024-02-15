#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DesiRNA.py

Usage:
Run the script with command-line arguments to specify input files, parameters for sequence design, and configuration options. The script supports various options for detailed control over the RNA design process.

Example:
python DesiRNA.py [options]

Dependencies:
This script uses several external libraries and modules, such as RNA for RNA secondary structure prediction, and pandas for data handling. Ensure these dependencies are installed and properly configured in your Python environment.

Note:
For detailed documentation on command-line arguments and options, refer to the accompanying README file or the help command (python DesiRNA.py -h).

Author: Tomasz Wirecki
Version: 11.2023
"""

import argparse as argparse
import sys
import random
import csv
import time
import re
import os

from datetime import datetime
from pathlib import Path
from shutil import copy

import RNA

from utils import sequence_utils as seq_utils
from utils import replica_exchange_monte_carlo as remc
from utils import stats_inputs_outputs as sio


def print_action_group_help(parser, action_group, include_all=False, print_usage=True):
    """
    Prints help information for a specific action group within an argparse parser.

    Parameters:
    - parser (ArgumentParser): The argparse parser object.
    - action_group (ArgumentGroup): The action group within the parser for which to print help.
    - include_all (bool, optional): If True, includes all actions in the usage. Defaults to False.
    - print_usage (bool, optional): If True, prints the usage information. Defaults to True.

    Returns:
    - None: This function does not return a value but prints help information to the console.
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
        """
        Executes the action to display advanced help options.

        This method is called when the action is triggered in the command-line interface. It prints the help information for both standard and advanced options, and then exits the parser.

        Parameters:
        - parser (ArgumentParser): The argparse parser instance.
        - namespace (Namespace): An object where the parser will add attribute data.
        - values: The associated command-line arguments; not used in this method.
        - option_string (str, optional): The option string that was used to trigger this action.

        """
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
    scoring_f = sio.parse_scoring_functions(args.scoring_f)
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

    sim_options = DesignOptions(infile, replicas, timlim, acgu_percs, t_max, t_min, oligo, param, exchange_rate, scoring_f, pm, tshelves,\
                                in_seed, subopt, diff_start_replicas, num_results, acgu_content, steps, tm_max, tm_min, motifs, dimer, sws)

    return sim_options


def initialize_simulation(input_file):
    """
    Initialize the simulation by setting up the necessary parameters and constraints.

    Args:
    input_file (InputFile): An object encapsulating input file parameters.

    Returns:
    list: List of nucleotides with applied constraints.
    """

    input_file.pairs = seq_utils.check_dot_bracket(input_file.sec_struct)  # check dotbracket correctness, assign as list of pairs
    input_file.set_target_pairs_tupl()

    if input_file.alt_sec_structs != None:
        print("Alternative structures are ok.")
        alt_pairs = seq_utils.get_pairs_for_graphs(input_file)
        input_file.graphs = seq_utils.generate_graphs(alt_pairs)
        seq_utils.update_graphs(input_file)


    seq_utils.check_seq_restr(input_file.seq_restr)
    seq_utils.check_length(input_file.sec_struct, input_file.seq_restr)
    nt_list = seq_utils.get_nt_list(input_file)
    seq_utils.check_input_logic(nt_list)
    
   # for i in range(len(nt_list)):
   #     print(vars(nt_list[i]))



    return nt_list


def generate_sequences(nt_list, input_file, sim_options):
    """
    Generates initial RNA sequences based on nucleotide constraints.

    Args:
    nt_list (list): A list of nucleotide constraints.
    input_file (InputFile): The input file object with sequence data.
    sim_options (DesignOptions): An object holding all the simulation options provided by the user.

    Returns:
    list: A list of initialized sequence objects.
    """

    sequence_score_list = seq_utils.generate_initial_list(nt_list, input_file, sim_options)
    seq_utils.generate_initial_list_random(nt_list, input_file, sim_options)

    return sequence_score_list


def handle_non_mutable_sequence(input_file, simulation_data, sim_options):
    """
    Handles the scenario where the sequence cannot be mutated due to sequence constraints.
    Scores the sequence and writes the results to a CSV file, then exits the program.

    Args:
    input_file (InputFile): The input file object containing sequence data.
    simulation_data (list): The simulation data to be processed.
    sim_options (DesignOptions): The object containing simulation options, including the output file name.
    """

    if all(char in "ACGU&" for char in input_file.seq_restr):
        print("Cannot mutate the sequence due to sequence constraints.\nScoring sequence.")
        sorted_results = sorted(seq_utils.round_floats(simulation_data), key=lambda d: (-d['mcc'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
        with open(sim_options.outname + '_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = sorted_results[0].keys()  # header from keys of the first dictionary
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for data in sorted_results:
                writer.writerow(data)
        sys.exit()


def run_functions(input_file, sim_options, now):
    """
    This function runs the design process, which includes mutation, scoring, and selection.

    Parameters:
    input_file (InputFile): An object that encapsulates the input parameters for the simulation.
    sim_options (DesignOptions): An object containing options and configurations for the simulation process.
    now (datetime): The current datetime, used for tracking the start of the simulation process.

    Returns:
    None: This function initiates the design process but does not return any value.
    """

    nt_list = initialize_simulation(input_file)
    seqence_score_list = generate_sequences(nt_list, input_file, sim_options)

    start_time = time.time()
    simulation_data = []

    for i in range(len(seqence_score_list)):
        simulation_data.append(vars(seqence_score_list[i]))

    handle_non_mutable_sequence(input_file, simulation_data, sim_options)

    stats = sio.Stats()

    while time.time() - start_time < sim_options.timlim:

        print('Designing sequences...\nETA', round((sim_options.timlim - (time.time() - start_time)), 0), 'seconds', end='\r')

        stats.update_global_step()
        seqence_score_list, stats = remc.mutate_sequence_re(seqence_score_list, nt_list, stats, sim_options, input_file)

        seqence_score_list, stats = remc.replica_exchange(seqence_score_list, stats, sim_options)

        for i in range(len(seqence_score_list)):
            seqence_score_list[i].get_sim_step(stats.step)
            simulation_data.append(vars(seqence_score_list[i]))

        if stats.global_step % 10 == 0:
            mcc_zero = sio.process_intermediate_results(simulation_data, sim_options)
            if mcc_zero == sim_options.num_results and sim_options.sws == "on":
                break

        if stats.global_step == sim_options.RE_steps:
            break
    finish_time = time.time() - start_time

    sio.parse_and_output_results(simulation_data, input_file, stats, finish_time, sim_options, now)


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


def update_options(sim_options, input_file_obj):
    """
    Updates the simulation options based on the input file and provided configurations.

    This function adjusts various simulation parameters, including time limits, ACGU content,
    secondary structure options, and temperature shelves, based on the input file and user-defined options.

    Parameters:
    sim_options (DesignOptions): An object containing simulation options and configurations.
    input_file_obj (InputFile): An object encapsulating the input file data, such as sequence and structure.

    Returns:
    tuple: A tuple containing the updated simulation options and the basename of the input file.
    """

    if sim_options.RE_steps != None:
        sim_options.update_timlim(100000000000000000)

    if sim_options.param == '1999':
        RNA.params_load(os.path.join(script_path, "rna_turner1999.par"))

    if sim_options.acgu_content == '':
        sim_options.add_nt_perc({"A": 15, "C": 30, "G": 30, "U": 15})
    else:
        acgu_l = [int(x) for x in sim_options.acgu_content.split(',')]
        if sum(acgu_l) != 100:
            print("The ACGU content should sum up to 100, check your command.")
            sys.exit()
        sim_options.add_nt_perc({"A": acgu_l[0], "C": acgu_l[1], "G": acgu_l[2], "U": acgu_l[3]})

    if input_file_obj.alt_sec_struct != None:
        sim_options.add_alt_ss("on")
    else:
        sim_options.add_alt_ss("off")

    oligo_opt = "none"

    if "&" in input_file_obj.sec_struct and sim_options.dimer == "off":
        too_much = input_file_obj.sec_struct.count("&")
        if too_much != 1:
            print("\nToo much structures in the input. Can only design RNA complexes of max two sequences.\nPlease correct your input file.\n")
            sys.exit()
        oligo_opt = "heterodimer"
    elif "&" in input_file_obj.sec_struct and sim_options.dimer == "on":
        oligo_opt = "homodimer"
    elif "&" not in input_file_obj.sec_struct and sim_options.oligo == "on":
        oligo_opt = "avoid"

    sim_options.add_oligo_state(oligo_opt)

    if set(input_file_obj.sec_struct).issubset('.()&'):
        sim_options.add_pks("off")
    else:
        sim_options.add_pks("on")

    filename = os.path.basename(sim_options.infile)

    sim_options.add_outname(sio.get_outname(filename, sim_options))

    if sim_options.tshelves == '':
        rep_temps_shelfs_opt = seq_utils.get_rep_temps(sim_options)
    else:
        rep_temps_shelfs_opt = [float(temp) for temp in sim_options.tshelves.split(",")]

    sim_options.add_rep_temps_shelfs(rep_temps_shelfs_opt)

    return sim_options, filename


class DesignOptions:
    """
    A class to store and manage design options for RNA sequence simulations.

    Attributes:
    infile (str): Input file path.
    replicas (int): Number of replicas in the simulation.
    timlim (int): Time limit for the simulation run.
    acgu_percentages (dict): Percentages of A, C, G, U nucleotides.
    T_max (float): Maximum temperature for replicas.
    T_min (float): Minimum temperature for replicas.
    oligo (str): Oligomerization option.
    param (str): Parameter set for RNA folding (e.g., '1999' for Turner 1999 parameters).
    RE_attempt (int): Frequency of replica exchange attempts.
    scoring_f (str): Scoring function for sequence evaluation.
    point_mutations (str): Type of point mutations allowed.
    tshelves (str): Temperature shelves.
    in_seed (int): Seed for random number generator.
    subopt (bool): Suboptimal design option.
    diff_start_replicas (bool): Start replicas with different sequences.
    num_results (int): Number of results to output.
    acgu_content (str): ACGU content configuration.
    RE_steps (int): Number of steps for replica exchange.
    tm_max (float): Maximum threshold for targeted mutation percentage.
    tm_min (float): Minimum threshold for targeted mutation percentage.
    motifs (str): Specific motifs to include in the design.
    dimer (str): Dimerization option.
    sws (str): Stop when solved option.
    L (float): Constant for Monte Carlo algorithm.

    Methods are designed to update various aspects of these attributes.
    """

    def __init__(self, infile, replicas, timlim, acgu_percentages, T_max, T_min, oligo, param, RE_attempt, scoring_f, point_mutations,
                 tshelves, in_seed, subopt, diff_start_replicas, num_results, acgu_content, RE_steps, tm_max, tm_min, motifs, dimer, sws):
        """
        Initializes the DesignOptions object with provided simulation parameters.
        Parameters are self-explanatory based on the attribute names.
        """
        self.infile = infile
        self.replicas = replicas
        self.timlim = timlim
        self.acgu_percentages = acgu_percentages
        self.T_max = T_max
        self.T_min = T_min
        self.oligo = oligo
        self.param = param
        self.RE_attempt = RE_attempt
        self.scoring_f = scoring_f
        self.point_mutations = point_mutations
        self.tshelves = tshelves
        self.in_seed = in_seed
        self.subopt = subopt
        self.diff_start_replicas = diff_start_replicas
        self.num_results = num_results
        self.acgu_content = acgu_content
        self.RE_steps = RE_steps
        self.tm_max = tm_max
        self.tm_min = tm_min
        self.motifs = motifs
        self.dimer = dimer
        self.sws = sws
        self.L = 504.12

    def update_timlim(self, timlim):
        """
        Updates the time limit for the simulation.
        Args:
        timlim (int): New time limit value.
        """
        self.timlim = timlim

    def add_nt_perc(self, nt_percentages):
        """
        Adds nucleotide percentages to the options.
        Args:
        nt_percentages (dict): Dictionary of nucleotide percentages.
        """
        self.nt_percentages = nt_percentages

    def add_oligo_state(self, oligo_state):
        """
        Updates the oligo state option.
        Args:
        oligo_state (str): New oligo state.
        """
        self.oligo_state = oligo_state

    def add_alt_ss(self, alt_ss):
        """
        Adds alternative secondary structure option.
        Args:
        alt_ss (str): Dot bracket of alternative structure.
        """
        self.alt_ss = alt_ss

    def add_pks(self, pks):
        """
        Adds pseudoknot inclusion state.
        Args:
        pks (str): State of pseudoknot inclusion ('on' or 'off').
        """
        self.pks = pks

    def add_rep_temps_shelfs(self, rep_temps_shelfs):
        """
        Updates the replica temperature shelves.
        Args:
        rep_temps_shelfs (list): List of temperature values.
        """
        self.rep_temps_shelfs = rep_temps_shelfs

    def add_outname(self, outname):
        """
        Sets the output filename.
        Args:
        outname (str): Output file name.
        """
        self.outname = os.path.basename(outname)


if __name__ == "__main__":

    #log_filename = "tmp_log.txt"
    #original_stderr, stderr_log = redirect_stderr_to_file(log_filename)

    #try:

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

    simulation_options = argument_parser()

    input_file_g = sio.read_input(simulation_options.infile)
    print(simulation_options.infile)

    simulation_options, filename = update_options(simulation_options, input_file_g)

    if simulation_options.in_seed != 0:
        original_seed = 2137 + simulation_options.in_seed
    else:
        original_seed = random.random()

    random.seed(original_seed)

    now = datetime.now()
    now = now.strftime("%Y%m%d.%H%M%S")

    WORK_DIR = simulation_options.outname + "_" + str(now)
    Path(WORK_DIR).mkdir(parents=True, exist_ok=True)
    copy(simulation_options.infile, WORK_DIR + "/" + filename)
    Path(WORK_DIR + "/trajectory_files").mkdir(parents=True, exist_ok=True)
    os.chdir(WORK_DIR)

    with open(simulation_options.outname + ".command", 'w', encoding='utf-8') as f:
        print(command, file=f)

    run_functions(input_file_g, simulation_options, now)
    '''    
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

    try:
        os.remove("../" + log_filename)
    except FileNotFoundError:
        pass
    '''
