#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import RNA
import numpy as np

kB = 0.001987204259
rho = 55.14       # H2O concentration in mol/L
temp = 273.15+37  # Physical Temperature
conc = 1e-3       # RNA concentration im mol/L
    
def oligo_fraction(seq_dimer):
    dimer_ss, pfa, pfb, pfab, dimer_pf = RNA.co_pf_fold(seq_dimer)
    dF = pfab - pfa - pfb
    rhs = conc/rho*np.exp(-dF/(kB*temp))
    return 1-(np.sqrt(1+4*rhs)-1)/(2*rhs)

def kTlog_oligo_fraction(oligo_frac):
    return -kB*temp*np.log(oligo_frac)

def kTlog_monomer_fraction(oligo_frac):
    return -kB*temp*np.log(1-oligo_frac)



def energy_of_oligomer(seq):
    fc = RNA.fold_compound(seq + "&" + seq)
    oligo_structure, mfe_oligo =fc.mfe_dimer()
    return mfe_oligo

def mfe_e_dimer(seq_dimer):
    fc = RNA.fold_compound(seq_dimer)
    mfe_ss_dimer_joint =fc.mfe_dimer()[0]
    
    seqa, seqb = seq_dimer.split("&")[0], seq_dimer.split("&")[1]

    cofolded = RNA.co_pf_fold(seqa+"&"+seqb)
    
    mfe_e_dimer = cofolded[-1]
    
    mfe_e_a = RNA.fold(seqa)[1]
    mfe_e_b = RNA.fold(seqb)[1]

    mfe_ss_dimer = mfe_ss_dimer_joint[:len(seqa)]+"&"+mfe_ss_dimer_joint[len(seqa):]

    return mfe_e_dimer, mfe_e_a, mfe_e_b, mfe_ss_dimer


def if_oligomer(mfe_mono, mfe_oligo):
    
    
    if 2*mfe_mono < mfe_oligo:
        oligo = False
    else:
        oligo = True   

    return oligo

    
def if_dimer(mfe_dimer, mfe_a, mfe_b):
    

    if (mfe_a+mfe_b) > mfe_dimer:
        dimer = True
    else:
        dimer = False   

    return dimer


if __name__ == "__main__":
    main()

