import random
import math

import multiprocess as mp

from utils import sequence_utils as seq_utils

def mc_delta(deltaF_o, deltaF_m, T_replica, sim_options):
    """
    Determine whether a mutation is accepted based on the Metropolis criterion.

    Args:
    T_replica (float): Temperature of the replica.
    deltaF_o (float): Original free energy.
    deltaF_m (float): Mutated free energy.

    Returns:
    tuple: A tuple containing two boolean values indicating whether the mutation is accepted and whether the mutated scoring function is better.
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
    temp (float): Temperature.
    dE (float): Energy difference.

    Returns:
    float: Metropolis score.
    """

    score = math.exp((-sim_options.L / temp) * (dE))

    return score


def replica_exchange_attempt(T0, T1, dE0, dE1, sim_options):
    """
    Attempt to exchange replicas.

    Args:
    T0 (float): Temperature of the first replica.
    T1 (float): Temperature of the second replica.
    dE0 (float): Energy of the first replica.
    dE1 (float): Energy of the second replica.

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
        replicas (list): List of replica objects.
        stats_obj (object): Statistics object.

    Returns:
        tuple: A tuple containing the list of replicas after exchange and the updated statistics object.
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

    Parameters
    ----------
    sequence_o : ScoreSeq object
        The original sequence to be mutated.
    nt_list : list
        List of possible nucleotides for mutation.
    worker_stats : Stats object
        Object to store the statistics of the worker.

    Returns
    -------
    best_score : ScoreSeq object
        The best scored sequence after performing the replica design.
    worker_stats : Stats object
        Updated statistics of the worker.
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
    This function mutates a list of sequence objects in parallel.

    Parameters:
    lst_seq_obj (list): A list of sequence objects to mutate.
    nt_list (list): A list of possible nucleotides to use in mutation.
    stats_obj (object): A statistical object to keep track of the mutation process.

    Returns:
    list: A list of mutated sequence objects.
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
