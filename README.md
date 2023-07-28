# 2DesignRNA
Sure, I apologize for misunderstanding. Here's a full detailed explanation of each argument in Markdown format:

```markdown
## Usage

You can execute the `2DesignRNA.py` program using the following command:

```bash
python 2DesignRNA.py [arguments]
```

Here is a detailed explanation of each argument:

- `-f`, `--filename`: **(Required)** The name of the input file containing secondary structures and constraints. This file should be in the correct format for the program to interpret the secondary structures and constraints correctly.

- `-R`, `--replicas`: **(Optional, default = 10)** The number of replicas to create for the simulation. Increasing this number will likely result in a more diverse set of RNA sequences, but may increase computation time.

- `-e`, `--exchange`: **(Optional, default = 10)** The number of accepted exchanges for Metropolis Monte Carlo simulation. Increasing this number will likely result in a more diverse set of RNA sequences, but may increase computation time.

- `-t`, `--timelimit`: **(Optional, default = 60)** The time limit for running the program, in seconds. Increase this to allow the program more time to find optimal solutions.

- `-n`, `--number_of_sequences`: **(Optional, default = 10)** The number of best sequences that the program should return. Increasing this number will yield a greater number of resultant sequences, but may increase computation time.

- `-acgu`, `--ACGU`: **(Optional, default = 'off')** This option enables or disables the ACGU constraints in the sequence design. 

- `-pk`, `--PK`: **(Optional, default = 'off')** This option enables or disables pseudoknots in the sequence design.

- `-l`, `--lambda`: **(Optional, default = 1000)** This option sets the temperature factor for the Metropolis criterion used in the Monte Carlo simulation. 

- `-o`, `--oligomerization`: **(Optional, default = 'off')** This option enables or disables oligomerization in the sequence design.

- `-d`, `--dimer`: **(Optional, default = 'off')** This option enables or disables dimerization in the sequence design.

- `-p`, `--param`: **(Optional, default = '2004')** This option allows you to choose the energy parameter model. The default model is the Turner 2004 model. The '1999' option allows you to use the Turner 1999 model.

- `-sf`, `--scoring_function`: **(Optional, default = 'dmt')** This option allows you to choose the scoring function for the sequence design. The default is 'dmt'. Other options include 'mcc', 'mfe', and 'mix'.

To use the program with specific options, include the corresponding arguments in your command. For example, to execute the program with 20 replicas, your command would look like this:

```bash
python 2DesignRNA.py -f yourfile -R 20
```

Replace `yourfile` with the path to your input file and `20` with the number of replicas you want to use. If your file path contains spaces, please enclose it in quotation marks.

Remember, the more precise the input parameters, the more accurate the resulting RNA sequence designs will be.
```

Just replace `yourfile` with the path to your file when running the program. Also, remember to replace 'TBA' with the appropriate descriptions once they become available.









DESCRIPTION FROM PARNASSUS README BELOW



## Generation of NOISE families

We need to generate "fake" families in order to sort out noise while processing BLAST matrices. We will create 1000 RNA families, basing on a representative seqence of length 50 - 1000 nts. 
First we need to generate those sequences, and the design a family basing on them. The family will be gapped sequences, designed by DesiRNA with sec-struct restraint being prediction of RNAfold for the representative sequence.
The consensus sec struct for the fake family will be the same as the restraint used for generating the family. Then the whole procedure of mapping-predicting-enriching-translating will be performed as on teh regular Rfam families.

### 1. Generation of random sequences - representatives of fake families

**generate_random_sequences.py**

To generate 1000 fake families, we need 1000 inputs for DesiRNA. We want to generate 1000 inputs of sequences ranging 50-1000 nts, and the population of the family would be random number between 10-1000.
The proportions of AUGC in the representative sequences is also randomised, each nucleotide will have probability in range 10-35, e.g. A = 25, G = 30, C = 35, U = 30 (**note** teh probabiliteis may sum to > or < than 100, but this eventually is rescaled to sum probablities = 100%).

There is no input for the script, only flags. 

Command:

`generate_random_sequences.py -n 1000 -min 50 -max 1000`

where:
* `-n` is the number of generated fake families representatives
* `-min` and `-max` are the minimal and maximal length of randomly generated sequences

the output of the scipt is a series of files for the modified DesiRNA (`design_fake_families.py`):
`RF_rand_XXXX.input`


Procedure:
```
generate_random_sequences.py -n 250 -min 50 -max 300  -s 1
generate_random_sequences.py -n 250 -min 50 -max 500  -s 251
generate_random_sequences.py -n 250 -min 50 -max 800  -s 501
generate_random_sequences.py -n 250 -min 50 -max 1000  -s 751
```


### 2. Design of the fake families

**design_fake_families.py**

No we need to run the design of families. Each represenattive is used as an input for the modified DesiRNA, and the number of generated sequences for each family is randomized in range 10-1000 sequences.
Additionally to each sequence there is a randomized gap introduction ragnging 0-40% of sequence length. Each of the fake families will also have a "fake" consensus structure, which is the secondary structure restraint used in the design.

Example input for the script:
```
>name
RF_rand_1
>sec_struct
..((((...(((.................))).))))............................
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>number
14
```

Command:

`design_fake_families.py -f RF_rand_1.input`

The output is a stockholm file with gapped sequences and "fake" consensus structure:

`RF_rand_1.output`

```
RF_rand_1.seq.1  AGUC-UG-UCCCGUCCGAUG--CAGCAGUGGGCAUGACAAAGCGCUUUGGUAGGUAU-UAAAUAC
RF_rand_1.seq.2  UGC--CGAG-UUACCAC-AC-GA-CUAGAAACUGGAGAU-UUGU-GUAUUCUUUG-GGC--UAGA
RF_rand_1.seq.3  C-AUAUGC--CC-CGAU--AUC-AAAGGUG--A-UAU--G-AG-U-CGCUG--CCGCG-C-AAGA
RF_rand_1.seq.4  GCUGUUGCAGUUCCAGGAGAGACUCCGGUAACUAACAUUUGGCGCUAAUCAAGCACGAUUUAUGA
RF_rand_1.seq.5  CAAA-AGCAUUCG-AC-A--A-AA-G-G-GAAGUCU-U-UA-UUC-AU-A-UUGAGUAUCUU--A
RF_rand_1.seq.6  CGU-AA-U-CGG-UCAUC--A-AUAGA-GCC-GU-AAACUU--C-UA-GU--CUGCGG-C-A--A
RF_rand_1.seq.7  GAUCUCAAGCGCU-AGUAAACAGCAACCCGCGCGAGAACUUUUGUACAUGGGGGUGUGUGCUACU
RF_rand_1.seq.8  -CAAGU--GUUGUA-U-UGAUAU--UG-UC---A-UUAAGCAA--U--A--U-GAGGAA--A-UG
RF_rand_1.seq.9  G--UGCGA-A-GUUGA-CC-UC--CACUUCGU-GCA-UA-G--A--C---GGC-CCG-CUGG--C
RF_rand_1.seq.10 CAGGAUAGGCUGAUGC-GUUCGGAGGCGACAGGAUCCAACGAUGGCAGUUCCAGACACUUUACAC
RF_rand_1.seq.11 GGG-AGUGA-GGA--A-A-CACUCGUC-GCCCU-U-CGAU---CG--GC-GACG-GGU--UGG-U
RF_rand_1.seq.12 CG-GCAGC-CGCCG-U---CU-GC-GG-A-CG-U-CGACG-GUCGCGCCUGCC--UAACUCGACU
RF_rand_1.seq.13 -U-CA--UGGGUAAG---AUU-GUCU-AAA-CGUUGCG--GAACGU-UGCCAAC--C-GGGCCGC
RF_rand_1.seq.14 UUACUUCACAUCCGACGCUUGCC-AUAGGGA-UAAGUGAGUGACUGA-UUCCACUCU-GAUAUUG
#=GC_SS_cons     ..((((...(((.................))).))))............................
```

