"""
stats_inputs_outputs.py - Module for handling input, output, and statistics in a simulation process.

This module contains functions and classes for reading input data, processing simulation results,
and generating various output files. It also includes a Stats class for tracking statistical metrics
during a simulation.

Functions:
- get_replica_color(replica_num, num_replicas): Determine the color for a given replica based on its number.
- plot_simulation_data_combined(simulation_data, alt, sim_options): Plot combined simulation data for various metrics across different replicas.
- read_input(infile): Reads the input file and returns an InputFile object.
- parse_scoring_functions(scoring_f_str): Parses the scoring function string into a list of scoring function tuples.
- process_intermediate_results(simulation_data, sim_options): Processes intermediate results during the simulation.
- generate_replica_csv(simulation_data, output_name): Generates a CSV file containing simulation data for each replica.
- sort_trajectory(simulation_data): Sorts the simulation data based on simulation step and replica number.
- generate_trajectory_csv(sorted_data, output_name): Generates a CSV file from the sorted simulation data.
- generate_multifasta(sorted_data, sim_options, now): Generates a multi-FASTA file from the sorted simulation data.
- generate_best_fasta(simulation_data, sim_options, now): Generates a FASTA file from the best sorted simulation data.
- sort_and_filter_simulation_data(simulation_data, sim_options): Sorts and filters the simulation data based on specified criteria.
- generate_csv_from_data(sorted_data, output_name): Generates a CSV file from sorted simulation data.
- check_if_design_solved(sorted_results, input_file, sim_options): Checks the correctness of the results.
- write_best_str_file(correct_result_txt, output_name): Writes the correctness text to a file.
- generate_simulation_stats_text(stats, sorted_results, input_file, correct_bool, finish_time, sim_options): Generates text summarizing the statistics of the simulation.
- write_stats_to_file(stats_txt, output_name): Writes simulation statistics text to a file.
- move_results(file_list, directory, output_name): Moves a list of files to the specified directory.
- parse_and_output_results(simulation_data, input_file, stats, finish_time, sim_options, now): Parses simulation data and generates various outputs including CSV files, FASTA files, and statistics.
- get_outname(infile, sim_options): Generate an output name for a file based on the current parameters.

Classes:
- Stats: A class used to track various statistical metrics during a simulation process.
"""

import csv
import time
import sys
from shutil import move
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import ticker

from utils import sequence_utils as seq_utils

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
    alt (bool): Flag to determine alternate plotting metrics.
    sim_options (DesignOptions): A DesignOptions object with simulation settings, including motif definitions.

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


def process_intermediate_results(simulation_data, sim_options):
    """
    Processes intermediate results during the simulation by removing duplicates,
    sorting the results, and writing them to a CSV file.

    Args:
    simulation_data (list): The simulation data to be processed.
    sim_options (DesignOptions): A DesignOptions object with simulation settings, including sorting criteria.

    Returns:
    int: The number of results with MCC (Matthews Correlation Coefficient) equal to zero.
    """

    remove_duplicates = {item['sequence']: item for item in simulation_data}
    data_no_duplicates = list(remove_duplicates.values())

    if sim_options.oligo != "off":
        sorted_results = sorted(seq_utils.round_floats(data_no_duplicates), key=lambda d: (-d['mcc'], -d['oligo_fraction'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
    elif sim_options.dimer != "off":
        sorted_results = sorted(seq_utils.round_floats(data_no_duplicates), key=lambda d: (-d['mcc'], d['oligo_fraction'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
    else:
        sorted_results = sorted(seq_utils.round_floats(data_no_duplicates), key=lambda d: (-d['mcc'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)

    sorted_results = sorted_results[:sim_options.num_results]
    mcc_zero_count = sum(1 for item in sorted_results if item['mcc'] == 0.0)

    with open(sim_options.outname + '_mid_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
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
    return sorted(seq_utils.round_floats(simulation_data), key=lambda d: (d['sim_step'], d['replica_num']))


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


def generate_multifasta(sorted_data, sim_options, now):
    """
    Generates a multi-FASTA file from the sorted simulation data.

    Args:
    sorted_data (list): Sorted simulation data.
    sim_options (DesignOptions): A DesignOptions object with simulation settings, including input file name.
    now (str): Timestamp or unique identifier.

    This function creates a multi-FASTA file containing sequences from the sorted simulation data.
    Each sequence is annotated with information including the input file name, timestamp, replica number,
    simulation step, and scoring function.
    The multi-FASTA file is saved with a specified output name.
    """

    fasta_txt = ""
    for item in sorted_data:
        fasta_txt += f">{sim_options.infile}|{now}|{item['replica_num']}|{item['sim_step']}|{item['scoring_function']}\n{item['sequence']}\n"

    with open(sim_options.outname + '_multifasta.fas', 'w', encoding='utf-8') as fastafile:
        fastafile.write(fasta_txt)


def generate_best_fasta(simulation_data, sim_options, now):
    """
    Generates a FASTA file from the best sorted simulation data.

    Args:
    simulation_data (list): The simulation data to be processed.
    sim_options (DesignOptions): A DesignOptions object with simulation settings, including input file name.
    now (str): Timestamp or identifier for the sequence.

    This function generates a FASTA file containing sequences from the top sorted simulation data based on scoring function.
    The number of top results to include in the FASTA file is determined by 'num_results' in 'sim_options'.
    Each sequence is annotated with information including the input file name, timestamp, replica number,
    simulation step, and scoring function.
    The FASTA file is saved with a specified output name.
    """

    sorted_scores = sorted(seq_utils.round_floats(simulation_data), key=lambda d: (-d['scoring_function']), reverse=True)[:sim_options.num_results]

    fasta_txt = "".join(f">{sim_options.infile}|{now}|{score['replica_num']}|{score['sim_step']}|{score['scoring_function']}\n{score['sequence']}\n" for score in sorted_scores)

    with open(sim_options.outname + '_best_fasta.fas', 'w', encoding='utf-8') as fastafile:
        fastafile.write(fasta_txt)


def sort_and_filter_simulation_data(simulation_data, sim_options):
    """
    Sorts and filters the simulation data based on specified criteria.

    Args:
    simulation_data (list): The simulation data to be processed.
    sim_options (DesignOptions): A DesignOptions object with simulation settings, including sorting and filtering criteria.

    Returns:
    list: Sorted and filtered list of simulation data.

    This function takes a list of simulation data and sorts it based on the criteria specified in 'sim_options'.
    If 'oligo' in 'sim_options' is not "off", the data is sorted by multiple criteria including oligo_fraction, mcc, edesired_minus_Epf, and Epf.
    Otherwise, it's sorted by mcc, edesired_minus_Epf, and Epf.
    Duplicate sequences are removed, and the top 'num_results' results are retained.
    The sorted and filtered list of simulation data is returned.
    """

    remove_duplicates = {item['sequence']: item for item in simulation_data}
    data_no_duplicates = list(remove_duplicates.values())

    if sim_options.oligo != "off":
        sorted_results = sorted(seq_utils.round_floats(data_no_duplicates), key=lambda d: (-d['mcc'], -d['oligo_fraction'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
    elif sim_options.dimer != "off":
        sorted_results = sorted(seq_utils.round_floats(data_no_duplicates), key=lambda d: (-d['mcc'], d['oligo_fraction'], -d['edesired_minus_Epf'], -d['Epf']), reverse=True)
    else:
        sorted_results = sorted(seq_utils.round_floats(data_no_duplicates), key=lambda d: (-d['mcc'],  -d['edesired_minus_Epf'], -d['Epf'], -d['scoring_function']), reverse=True)

    return sorted_results[:sim_options.num_results]


def sort_and_filter_simulation_data_slternative(simulation_data, sim_options, input_file):

    remove_duplicates = {item['sequence']: item for item in simulation_data}
    data_no_duplicates = list(remove_duplicates.values())
    sorted_results = sorted(seq_utils.round_floats(data_no_duplicates), key=lambda d: (-d['mcc'],  -d['edesired_minus_Epf'], -d['Epf'], -d['scoring_function']), reverse=True)
    dat_alt_mcc = seq_utils.get_alt_mcc(sorted_results, input_file)

    mcc_list = []
    for i in range(len(input_file.alt_sec_structs)):
        mcc_list.append("mcc_"+str(i+1))

    generate_trajectory_csv(dat_alt_mcc, sim_options.outname)

    sorted_results = sorted(seq_utils.round_floats(dat_alt_mcc), 
                        key=lambda d: tuple([-d[mcc] for mcc in (['mcc'] + mcc_list)]+ [-d['mcc'], -d['edesired_minus_Epf'], -d['Epf']]),
                        reverse=True)



    return sorted_results[:sim_options.num_results]


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


def check_if_design_solved(sorted_results, input_file, sim_options):
    """
    Checks the correctness of the results and generates the corresponding text output.

    Args:
    sorted_results (list): The sorted results to check.
    input_file (InputFile): The input file object containing sequence data.
    sim_options (DesignOptions): A DesignOptions object with simulation settings, including the oligo option.

    Returns:
    Tuple[str, bool, int]: A tuple containing three elements:
        - str: Text indicating the correctness of the results.
        - bool: Boolean indicating whether at least one of the top 10 results has an MCC of 0.0.
        - int: The count of correct results among the top 10.

    This function checks the correctness of the results by examining the MCC (Matthews Correlation Coefficient) values
    of the top 10 sequences in 'sorted_results'. If at least one of the top 10 results has an MCC of 0.0, it is considered correct.
    The function generates text output indicating whether the design is correct, the count of correct results, and additional
    information such as oligo fraction if specified in 'sim_options'.
    """

    correct_count = sum(result['mcc'] == 0.0 for result in sorted_results[:10])
    correct_bool = correct_count > 0

    if sim_options.oligo != "off":
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


def generate_simulation_stats_text(stats, sorted_results, input_file, correct_bool, finish_time, sim_options):
    """
    Generates text summarizing the statistics of the simulation.

    Args:
    stats (Stats): Object containing simulation statistics.
    sorted_results (list): The sorted results for correctness check.
    input_file (InputFile): Object containing input file data.
    correct_bool (bool): Indicates if the design was solved successfully.
    finish_time (float): Finish time of the simulation.
    sim_options (DesignOptions): A DesignOptions object with simulation settings.

    Returns:
    str: Text summarizing the simulation statistics.

    This function generates text summarizing various statistics of the simulation, including acceptance ratios,
    iteration counts, replica exchange statistics, and the best solution (if the design was solved).
    The function returns the generated text as a string.
    """

    sum_mc = stats.acc_mc_step + stats.rej_mc_step
    acc_perc = round(stats.acc_mc_step / sum_mc, 3) if sum_mc else 0
    sum_mc_metro = sum_mc - stats.acc_mc_better_e
    acc_metro = stats.acc_mc_step - stats.acc_mc_better_e
    sum_replica_att = stats.acc_re_step + stats.rej_re_step if stats.acc_re_step + stats.rej_re_step else 1
    formatted_time = time.strftime("%H:%M:%S", time.gmtime(finish_time))

    best_solution_txt = "\nDesign solved succesfully!\n\nBest solution:\n" if correct_bool else "\nDesign not solved!\n\nTarget structure:\n"
    best_solution_txt += f"{sorted_results[0]['sequence']}\nMFE Secondary Structure: {sorted_results[0]['mfe_ss']}\nPartition Function Energy: \
                        {round(sorted_results[0]['Epf'], 3)}\n1-MCC: {round(sorted_results[0]['mcc'], 3)}\n"

    stats_txt = f">{sim_options.outname} time={sim_options.timlim}s\nAcc_ratio={acc_perc}, Iterations={stats.step}, Accepted={stats.acc_mc_step}/{sum_mc}, Rejected={stats.rej_mc_step}/{sum_mc}\n" \
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


def parse_and_output_results(simulation_data, input_file, stats, finish_time, sim_options, now):
    """
    Parses simulation data and generates various outputs including CSV files, FASTA files, and statistics.

    Args:
    simulation_data (list): The simulation data to be processed.
    input_file (InputFile): The input file object containing sequence data.
    stats (Stats): Object containing simulation statistics.
    finish_time (float): Finish time of the simulation.
    sim_options (DesignOptions): A DesignOptions object with simulation settings.
    now (str): Timestamp or unique identifier.

    This function processes simulation data and generates various output files and statistics. It sorts the simulation
    trajectory, generates CSV and FASTA files, checks for design correctness, and produces simulation statistics.
    The generated files and statistics are saved to the specified output directory.
    """

    # Generate necessary files and statistics
    sorted_simulation_trajectory = sort_trajectory(simulation_data)
    generate_trajectory_csv(sorted_simulation_trajectory, sim_options.outname)
    generate_multifasta(sorted_simulation_trajectory, sim_options, now)
    generate_replica_csv(simulation_data, sim_options.outname)
    generate_best_fasta(simulation_data, sim_options, now)

    if input_file.graphs == None:
        sorted_results = sort_and_filter_simulation_data(simulation_data, sim_options)
    else:
        sorted_results = sort_and_filter_simulation_data_slternative(simulation_data, sim_options, input_file)
    generate_csv_from_data(sorted_results, sim_options.outname)

    correct_result_txt, correct_bool, correct = check_if_design_solved(sorted_results[:10], input_file, sim_options)
    write_best_str_file(correct_result_txt, sim_options.outname)
    plot_simulation_data_combined(simulation_data, input_file.alt_sec_struct, sim_options)

    stats_txt = generate_simulation_stats_text(stats, sorted_results, input_file, correct_bool, finish_time, sim_options)
    print('\n' + stats_txt)
    write_stats_to_file(stats_txt, sim_options.outname)

    # Move results to the specified directory
    files_to_move = ['_best_str', '_traj.csv', '_replicas.csv', '_best_fasta.fas', '_multifasta.fas', '_random.csv', '.command']
    move_results(files_to_move, "trajectory_files/", sim_options.outname)


def get_outname(infile, sim_options):
    """
    Generate an output name for a file based on the current parameters.

    The output name includes various simulation parameters, such as the name of the input file, the number of replicas,
    energy settings, time limit, PKS setting, ACGU percentages, temperature range, parameter setting, scoring function,
    oligomerization setting, dimer setting, mutation setting, and point mutation setting.

    Args:
    infile (str): The name of the input file.
    sim_options (DesignOptions): A DesignOptions object with simulation settings.

    Returns:
    str: The generated output name.

    This function generates an output name for files based on the provided input file name and simulation options. The
    output name is constructed using a combination of these parameters and serves as a unique identifier for the simulation.
    """

    available_scoring_functions = ['Ed-Epf', '1-MCC', 'sln_Epf', 'Ed-MFE', '1-precision', '1-recall', 'Edef']

    for i in range(0, len(sim_options.scoring_f)):
        if sim_options.scoring_f[i][0] not in available_scoring_functions:
            print(sim_options.scoring_f[i][0], "is not an available option for scoring function. Check your command.")
            sys.exit()

    scoring_f_str = "_".join([f"{func}_{weight}" for func, weight in sim_options.scoring_f])

    outname = sim_options.infile.split(".")[0] + '_R' + str(sim_options.replicas) + "_e" + str(sim_options.RE_attempt) + "_t" +\
        str(sim_options.timlim) + "_pk" + str(sim_options.pks) + "_ACGU" + str(sim_options.acgu_percentages) +\
        "_Tmin" + str(sim_options.T_min) + "_Tmax" + str(sim_options.T_max) + "_p" + str(sim_options.param) + "_SF" +\
        scoring_f_str + "_O" + (str(sim_options.oligo)) + "_D" + (str(sim_options.dimer)) + "_PM" + (str(sim_options.point_mutations))
    return outname


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
        self.alt_pairs = None
        self.seed_seq = None
        self.alt_sec_struct = None
        self.alt_sec_structs = None
        self.target_pairs_tupl = {}
        self.graphs = None
        self.excluded_alt_pairs = None
        self.allsnakes = None

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
        """
        Converst pair list to tuples list.
        """
        self.target_pairs_tupl = {tuple(pair) for pair in self.pairs}
