#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module provides functionalities for analyzing RNA sequences, specifically for calculating various thermodynamic properties and minimum free energies (MFEs) of RNA oligomers and dimers. It utilizes the RNA package for RNA secondary structure prediction and energy calculation.

The module includes functions to calculate the oligomer fraction of a dimer sequence, the kT logarithm of both oligomer and monomer fractions, the MFE of an RNA oligomer, and the MFE and individual component energies of an RNA dimer.

Dependencies:
    - RNA: For RNA secondary structure prediction and energy calculations.
    - numpy: For numerical operations.

Constants:
    - KB: Boltzmann constant in kcal/(mol*K).
    - RHO: Concentration of H2O in mol/L.
    - TEMP: Physical temperature in Kelvin.
    - CONC: RNA concentration in mol/L.

Functions:
    - oligo_fraction: Calculates the oligomer fraction of a dimer sequence.
    - kTlog_oligo_fraction: Computes the kT logarithm of the oligomer fraction.
    - kTlog_monomer_fraction: Computes the kT logarithm of the monomer fraction.
    - energy_of_oligomer: Calculates the MFE of an RNA oligomer.
    - mfe_e_dimer: Calculates the MFE and individual component energies of an RNA dimer.
"""

import RNA
import numpy as np

KB = 0.001987204259
RHO = 55.14       # H2O CONCentration in mol/L
TEMP = 273.15 + 37  # Physical Temperature
CONC = 1e-3       # RNA CONCentration im mol/L


def oligo_fraction(seq_dimer, fc):
    """
    Calculate the oligomer fraction of a dimer sequence.

    Parameters:
    seq_dimer (str): The RNA dimer sequence.
    fc (RNA.fold_compound): The fold compound object for RNA secondary structure.

    Returns:
    float: Oligomer fraction of the dimer.
    """
    dimer_ss, pfa, pfb, pfab, dimer_pf = fc.pf_dimer()  # RNA.co_pf_fold(seq_dimer)
    dF = pfab - pfa - pfb
    rhs = CONC / RHO * np.exp(-dF / (KB * TEMP))
    return 1 - (np.sqrt(1 + 4 * rhs) - 1) / (2 * rhs)


def kTlog_oligo_fraction(oligo_frac):
    """
    Calculate the kT logarithm of the oligomer fraction.

    Parameters:
    oligo_frac (float): The oligomer fraction.

    Returns:
    float: The kT logarithm of the oligomer fraction.
    """
    return -KB * TEMP * np.log(oligo_frac)


def kTlog_monomer_fraction(oligo_frac):
    """
    Calculate the kT logarithm of the monomer fraction.

    Parameters:
    oligo_frac (float): The oligomer fraction.

    Returns:
    float: The kT logarithm of the monomer fraction.
    """
    return -KB * TEMP * np.log(1 - oligo_frac)


def energy_of_oligomer(seq):
    """
    Calculate the minimum free energy (MFE) of an RNA oligomer.

    Parameters:
    seq (str): The RNA sequence.

    Returns:
    float: The MFE of the RNA oligomer.
    """
    fc = RNA.fold_compound(seq + "&" + seq)
    oligo_structure, mfe_oligo = fc.mfe_dimer()
    return mfe_oligo


def mfe_e_dimer(seq_dimer):
    """
    Calculate the MFE and individual component energies of an RNA dimer.

    Parameters:
    seq_dimer (str): The RNA dimer sequence, separated by '&'.

    Returns:
    tuple: Contains MFE of the dimer, MFE of component A, MFE of component B, and the secondary structure.
    """
    fc = RNA.fold_compound(seq_dimer)
    mfe_ss_dimer_joint = fc.mfe_dimer()[0]

    seqa, seqb = seq_dimer.split("&")[0], seq_dimer.split("&")[1]

    cofolded = RNA.co_pf_fold(seqa + "&" + seqb)

    mfe_e_dimer = cofolded[-1]

    mfe_e_a = RNA.fold(seqa)[1]
    mfe_e_b = RNA.fold(seqb)[1]

    mfe_ss_dimer = mfe_ss_dimer_joint[:len(seqa)] + "&" + mfe_ss_dimer_joint[len(seqa):]

    return mfe_e_dimer, mfe_e_a, mfe_e_b, mfe_ss_dimer
