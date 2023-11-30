import RNA

from utils import dimer_multichain_energy as dme
from utils import sequence_utils as seq_utils

from utils.sim_score import SimScore

md = RNA.md()
md.compute_bpp = 0


def score_sequence(seq, input_file, sim_options):
    """
    Calculate various scoring metrics for a given RNA sequence.

    Parameters
    ----------
    seq : str
        RNA sequence to score.

    Returns
    -------
    scored_sequence : ScoreSeq object
        Object with various scoring metrics, including mfe, target energy,
        energy difference, precision, recall, mcc, and overall score.
        If the 'oligo' option is set, the oligomerization of the sequence
        is also calculated.

    Notes
    -----
    This function uses the ViennaRNA package to calculate the mfe
    and the energy of the target structure. The precision, recall,
    and MCC are calculated using the SimScore class. The overall score
    is a combination of these metrics.
    """

    if sim_options.oligo_state != "homodimer":
        scored_sequence = ScoreSeq(sequence=seq)
    else:
        seq = seq.split("&")[0]
        seq = seq + "&" + seq
        scored_sequence = ScoreSeq(sequence=seq)
    pf_energy, mfe_structure, fold_comp = get_mfe_e_ss(seq, sim_options)

    scored_sequence.get_Epf(pf_energy)
    scored_sequence.get_mfe_ss(mfe_structure)

    scored_sequence.get_edesired(fold_comp.eval_structure(input_file.sec_struct.replace("&", "")))

    scored_sequence.get_edesired_minus_Epf(scored_sequence.Epf, scored_sequence.edesired)

    ssc = SimScore(input_file.sec_struct.replace("&", "Ee"), scored_sequence.mfe_ss.replace("&", "Ee"))
    ssc.find_basepairs()
    ssc.cofusion_matrix()


    scored_sequence.get_precision(ssc.precision())
    scored_sequence.get_recall(ssc.recall())
    scored_sequence.get_mcc(ssc.mcc())

    for function, weight in sim_options.scoring_f:
        if function == 'sln_Epf':
            scored_sequence.get_sln_Epf()
        if function == 'Ed-MFE':
            scored_sequence.get_MFE()
            scored_sequence.get_edesired_minus_MFE()
        if function == 'Edef':
            scored_sequence.get_ensemble_defect(input_file.sec_struct)

    scored_sequence.get_scoring_function(sim_options.scoring_f)

    if input_file.alt_sec_struct != None:
        energies = [fold_comp.eval_structure(alt_dbn) for alt_dbn in input_file.alt_sec_structs]
        scored_sequence.get_edesired2(sum(energies) / len(energies))
        scored_sequence.get_edesired2_minus_Epf(scored_sequence.Epf, scored_sequence.edesired2)
        scored_sequence.get_scoring_function_w_alt_ss()

    if (sim_options.subopt == "on") and (scored_sequence.mcc == 0):
        scored_sequence.get_subopt_e(func.get_first_suboptimal_structure_and_energy(seq, fold_comp)[1])
        scored_sequence.get_esubopt_minus_Epf(scored_sequence.Epf, scored_sequence.subopt_e)
        scored_sequence.get_scoring_function_w_subopt()

    if sim_options.oligo_state == "homodimer" or sim_options.oligo_state == "heterodimer":
        scored_sequence.get_scoring_function_oligomer(fold_comp)

    if sim_options.oligo_state == "avoid":
        scored_sequence.get_scoring_function_monomer()

    if sim_options.motifs:
        motif_score = score_motifs(seq, sim_options)
        scored_sequence.update_scoring_function_w_motifs(motif_score)

    return scored_sequence


def get_mfe_e_ss(seq, sim_options):
    """
    This function calculates the minimum free energy (MFE) and ensemble free energy (EFE)
    of a given RNA sequence.
    Parameters:
    seq (str): The RNA sequence to analyze.
    Returns:
    (float, float): A tuple containing the MFE and EFE of the sequence.
    """
    fc = RNA.fold_compound(seq, md)

    if sim_options.oligo_state in {"none", "avoid"}:
        structure, energy = fc.pf()
        structure = fc.mfe()[0]
        if sim_options.pks == "on":
            structure = seq_utils.get_pk_struct(seq, structure, fc)
    elif sim_options.oligo_state in {"homodimer", "heterodimer"}:
        seqa_len = len(seq.split("&")[0])
        structure_dim, energy = fc.mfe_dimer()
        structure = structure_dim[:seqa_len] + "&" + structure_dim[seqa_len:]
    return energy, structure, fc


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
        self.oligo_fraction = dme.oligo_fraction(self.sequence + "&" + self.sequence, fc)
        self.monomer_bonus = dme.kTlog_monomer_fraction(self.oligo_fraction)
        self.scoring_function = self.scoring_function + self.monomer_bonus

    def get_scoring_function_oligomer(self, fc):
        """
        Update the scoring function for this sequence by including the effect of oligomer fraction.

        Parameters:
        fc (RNA.fold_compound): The fold compound object for RNA secondary structure.
        """
        self.oligo_fraction = dme.oligo_fraction(self.sequence, fc)
        self.oligomer_bonus = dme.kTlog_oligo_fraction(self.oligo_fraction)
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
