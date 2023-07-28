#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import RNA


def energy_of_oligomer(seq):
    fc = RNA.fold_compound(seq + "&" + seq)
    oligo_structure, mfe_oligo =fc.mfe_dimer()
    return mfe_oligo

def mfe_e_dimer(seq_dimer):
    fc = RNA.fold_compound(seq_dimer)
#    mfe_ss_dimer_joint, mfe_e_dimer =fc.mfe_dimer()
    mfe_ss_dimer_joint =fc.mfe_dimer()[0]
    
#    mfe_ss_dimer_joint, mfe_e_dimer = RNA.pf_co_fold()
    
    seqa, seqb = seq_dimer.split("&")[0], seq_dimer.split("&")[1]

    print(seqa, seqb)
    cofolded = RNA.co_pf_fold(seqa+"&"+seqb)
#    print(cofolded,"cofo")
#    print(fc.mfe_dimer(),"fc")
    
    mfe_e_dimer = cofolded[-1]
    
    mfe_e_a = RNA.fold(seqa)[1]
    mfe_e_b = RNA.fold(seqb)[1]

    mfe_ss_dimer = mfe_ss_dimer_joint[:len(seqa)]+"&"+mfe_ss_dimer_joint[len(seqa):]

    print(mfe_e_dimer, mfe_e_a, mfe_e_b, mfe_ss_dimer)
#    quit()
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

