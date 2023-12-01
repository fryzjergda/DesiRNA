"""
This module contains functions for performing replica exchange Monte Carlo simulations.

Replica exchange Monte Carlo (REMC) is a simulation technique used in computational biology to explore
conformational space and predict the folding behavior of RNA molecules. This module provides functions
for replica exchange, Metropolis criterion, and mutation processes necessary for REMC simulations.

Functions:
- mc_delta(deltaF_o, deltaF_m, T_replica, sim_options): Determine whether a mutation is accepted based on the Metropolis criterion.
- metropolis_score(temp, dE, sim_options): Calculate the Metropolis score.
- replica_exchange_attempt(T0, T1, dE0, dE1, sim_options): Attempt to exchange replicas.
- replica_exchange(seq_score_list, stats_obj, sim_options): Perform the replica exchange.
- single_replica_design(sequence_o, nt_list, worker_stats, sim_options, input_file): Perform replica design for a single replica.
- par_wrapper(args): Wrapper function for parallel processing.
- mutate_sequence_re(lst_seq_obj, nt_list, stats_obj, sim_options, input_file): Mutate a list of sequence objects in parallel.
"""

import random
import math

import multiprocess as mp

from utils import sequence_utils as seq_utils


def mc_delta(deltaF_o, deltaF_m, T_replica, sim_options):
    """
    Determine whether a mutation is accepted based on the Metropolis criterion.

    Parameters:
        T_replica (float): Temperature of the replica.
            The temperature of the replica at which the mutation is being considered.
        deltaF_o (float): Original free energy.
            The free energy of the original sequence.
        deltaF_m (float): Mutated free energy.
            The free energy of the mutated sequence.
        sim_options (object): DesignOptions object.
            An object containing simulation options and parameters.

    Returns:
        tuple: A tuple containing two boolean values indicating whether the mutation is accepted and whether the mutated scoring function is better.
            - accept (bool): True if the mutation is accepted, False otherwise.
            - accept_e (bool): True if the mutated scoring function is better, False otherwise.
    """

    accept_e = False
    if deltaF_m <= deltaF_o:
        accept = True
        accept_e = True
    else:
        diff = deltaF_m - deltaF_o
        p_dE = metropolis_score(T_replica, diff, sim_options)
        rand_num = random.random()

        accept = p_dE > rand_num

    return accept, accept_e


def metropolis_score(temp, dE, sim_options):
    """
    Calculate the Metropolis score.

    Args:
        temp (float): Temperature of the replica.
            The temperature at which the Metropolis score is calculated.
        dE (float): Energy difference.
            The energy difference between two states or configurations.
        sim_options (object): DesignOptions object.
            An object containing simulation options and parameters.

    Returns:
        tuple: A tuple containing two boolean values indicating whether the mutation is accepted and whether the mutated scoring function is better.
    """

    score = math.exp((-sim_options.L / temp) * (dE))

    return score


def replica_exchange_attempt(T0, T1, dE0, dE1, sim_options):
    """
    Attempt to exchange replicas.

    Args:
        T0 (float): Temperature of the first replica.
            The temperature of the first replica.
        T1 (float): Temperature of the second replica.
            The temperature of the second replica.
        dE0 (float): Energy of the first replica.
            The energy of the first replica.
        dE1 (float): Energy of the second replica.
            The energy of the second replica.
        sim_options (object): DesignOptions object.
            An object containing simulation options and parameters.

    Returns:
        tuple: A tuple containing two boolean values indicating whether the exchange is accepted and whether the exchange resulted in a better energy.
    """

    accept_e = False
    if dE1 <= dE0:
        accept = True
        accept_e = True
    else:
        rand_num = random.random()
        p = math.exp(sim_options.L * (1 / T0 - 1 / T1) * (dE0 - dE1))
        accept = p > rand_num

    return accept, accept_e


def replica_exchange(seq_score_list, stats_obj, sim_options):
    """
    Perform the replica exchange.

    Args:
        seq_score_list (list): List of ScoreSeq objects.
            A list of sequence objects representing the replicas.
        stats_obj (object): Statistics object.
            An object for storing statistics related to the replica exchange.
        sim_options (object): DesignOptions object.
            An object containing simulation options and parameters.

    Returns:
        tuple: A tuple containing the list of replicas after exchange (seq_score_list) and the updated statistics object (stats_obj).
    """

    re_pairs = []

    num_shelfs = [i for i in range(1, len(seq_score_list))]

    temps = []

    for i in range(0, len(seq_score_list)):
        temps.append(seq_score_list[i].temp_shelf)

    temps = sorted(temps)

    if stats_obj.global_step % 2 == 0:

        for i in range(0, len(num_shelfs) - 1, 2):
            re_pairs.append([i + 1, i + 2])
    else:
        for i in range(0, len(num_shelfs), 2):
            re_pairs.append([i, i + 1])

    seq_score_list = sorted(seq_score_list, key=lambda obj: obj.temp_shelf)

    for i in range(len(re_pairs)):
        rep_i = re_pairs[i][0]  # replica temporary number
        rep_j = re_pairs[i][1]  # replica temporary number

        T_i = seq_score_list[rep_i].temp_shelf
        T_j = seq_score_list[rep_j].temp_shelf

        dE_i = seq_score_list[rep_i].scoring_function
        dE_j = seq_score_list[rep_j].scoring_function

        accept, accept_e = replica_exchange_attempt(T_i, T_j, dE_i, dE_j, sim_options)

        if accept == True:
            seq_score_list[rep_i].get_temp_shelf(T_j)
            seq_score_list[rep_j].get_temp_shelf(T_i)
            stats_obj.update_acc_re_step()
            if accept_e == True:
                stats_obj.update_acc_re_better_e()
        else:
            stats_obj.update_rej_re_step()

    seq_score_list = sorted(seq_score_list, key=lambda obj: obj.replica_num)

    return seq_score_list, stats_obj


def single_replica_design(sequence_o, nt_list, worker_stats, sim_options, input_file):
    """
    Perform the replica design for a single replica.

    Args:
        sequence_o (ScoreSeq object): The original sequence to be mutated.
        nt_list (list): List of possible nucleotides for mutation.
        worker_stats (Stats object): Object to store the statistics of the worker.
        sim_options (DesignOptions object): Object containing simulation options and parameters.
        input_file (object): Object containing input file data.

    Returns:
        tuple: A tuple containing the best scored sequence after performing the replica design (sequence_o) and the updated statistics of the worker (worker_stats).
    """

    worker_stats.reset_mc_stats()

    for i in range(0, sim_options.RE_attempt):
        sequence_m = seq_utils.mutate_sequence(sequence_o, nt_list, sim_options, input_file)

        deltaF_o = sequence_o.scoring_function
        deltaF_m = sequence_m.scoring_function

        accept, accept_e = mc_delta(deltaF_o, deltaF_m, sequence_o.temp_shelf, sim_options)

        if accept == True:
            sequence_o = sequence_m
            worker_stats.update_acc_mc_step()
            if accept_e == True:
                worker_stats.update_acc_mc_better_e()

        if accept == False:
            worker_stats.update_rej_mc_step()

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

    return single_replica_design(*args)


def mutate_sequence_re(lst_seq_obj, nt_list, stats_obj, sim_options, input_file):
    """
    Mutate a list of sequence objects in parallel.

    Args:
        lst_seq_obj (list): List of sequence objects to mutate.
        nt_list (list): List of possible nucleotides for mutation.
        stats_obj (Stats object): Object to keep track of mutation statistics.
        sim_options (DesignOptions object): Object containing simulation options and parameters.
        input_file (object): Object containing input file data.

    Returns:
        tuple: A tuple containing the list of mutated sequence objects and the updated statistics object.
    """

    with mp.Pool(sim_options.replicas) as pool:
        # Generate a unique seed for each worker based on original_seed
        seeds = [i for i in range(len(lst_seq_obj))]

        # inputs = [(seq_obj, nt_list, stats_obj) for seq_obj in lst_seq_obj]
        inputs = [(seq_obj, nt_list, stats_obj, sim_options, input_file, seed) for seq_obj, seed in zip(lst_seq_obj, seeds)]
        lst_seq_obj_res_new = []
        workers_stats = []

        for result in pool.imap(par_wrapper, inputs):
            lst_seq_obj_res_new.append(result[0])
            workers_stats.append(result[1])

        attributes_to_update = ['acc_mc_step', 'acc_mc_better_e', 'rej_mc_step']

        stats_obj.update_step(sim_options.RE_attempt)

        for stats in workers_stats:
            for attr in attributes_to_update:
                current_value = getattr(stats_obj, attr)
                additional_value = getattr(stats, attr)
                setattr(stats_obj, attr, current_value + additional_value)

        return lst_seq_obj_res_new, stats_obj
