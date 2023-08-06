# DesiRNA

DesiRNA is a state-of-the-art RNA sequence design tool, that stands out for its speed, lightweight nature, ease of installation, and user-friendly interface.

### Features
- **Super Fast:** Utilizes Replica Exchange Monte Carlo for efficient computations.
- **Installation:** Simple to install and requires only Python3.
- **Versatile Design Capabilities:**
  - Single chain RNAs
  - Two-strand RNA complexes
  - RNAs with alternative structures
  - RNAs with pseudoknots
  - Checks if the designed RNA tends to oligomerize
  - Designs RNAs with natural ACGU content
- **Positive and Negative Design:** Follows both positive and negative design principles.
- **Exceptional Performance:** Proven extraordinary performance in the EteRNA benchmark, solving 99 out of 100 design problems in under 24 hours.


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Files](#files)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)


## Installation

The installation of DesiRNA can be accomplished through several methods. Below are the instructions for different scenarios:

### Without Conda

Execute the following commands:

```bash
python3 -m pip install -r requirements.txt
python3 -m pip install viennarna
```

### With Conda - Creating a Custom Environment

Execute the following commands:

```bash
conda create -n desirna python=3 numpy matplotlib pandas shutil pathlib
conda activate desirna
pip install viennarna
```

### With Conda - Utilizing Conda Environment YAML File

Execute the following commands:

```bash
conda env create -f DesiRNA-env.yml
conda activate DesiRNA-env
```

#### Adding DesiRNA to PATH

To facilitate access to DesiRNA from the command line, the path to `DesiRNA.py` may be added to the `.bashrc` file:

```bash
echo 'export PATH=$PATH:/path/to/DesiRNA.py' >> ~/.bashrc
source ~/.bashrc
```

Replace `/path/to/DesiRNA.py` with the correct path to the `DesiRNA.py` file.

#### Additional Installation Guidance

- **Conda Installation:** Should Conda be required, the [official Conda installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) provides comprehensive instructions.
- **ViennaRNA Installation:** ViennaRNA may be installed through various methods. Detailed instructions are available in the [ViennaRNA installation guide](https://www.tbi.univie.ac.at/RNA/#download).




## Usage

DesiRNA is a command-line tool, and its functionality can be accessed through various options and arguments. Below is a detailed explanation of how to use DesiRNA:

### Basic Command

The basic command to run DesiRNA requires specifying the filename that contains secondary structures and constraints:

```bash
DesiRNA.py -f NAME
```

### Optional Arguments

- `-R, --replicas REPLICAS`: Number of replicas (default: 10).
- `-e, --exchange EXCHANGE`: Frequency of replica exchange attempt (default: 100).
- `-t, --timelimit TIMLIM`: Time limit for running the program in seconds (default: 60).
- `-s, --steps STEPS`: Number of Replica Exchange steps, overwrites the `-t` option (default: None).
- `-r, --results_number NUM_RESULTS`: Number of best results to be reported (default: 10).
- `-p, --param {2004,1999}`: Turner energy parameter for calculating MFE (default: 1999).
- `-tmin, --tmin T_MIN`: Minimal Replica Temperature (default: 10).
- `-tmax, --tmax T_MAX`: Maximal Replica Temperature (default: 150).
- `-ts, --tshelves TSHELVES`: Custom temperature shelves for replicas, provide comma-separated values.
- `-sf, --scoring_function {dmt,mcc}`: Scoring function used to guide the design process (default: dmt).
- `-nd, --negative_design {off,on}`: Use negative design approach (default: off).
- `-acgu, --ACGU {off,on}`: Keep 'natural' ACGU content, with default content A:15%, C:30%, G:30%, U:15% (default: off).
- `-acgu_content, --ACGU_content ACGU_CONTENT`: Provide user-defined ACGU content, comma-separated values (e.g., 10,40,40,10).
- `-o, --oligomerization {off,on}`: Check if the designed sequence tends to oligomerize (default: off).
- `-tm, --target_mutations {off,on}`: Targeted mutations (default: on).
- `-tm_perc_max, --target_mutations_percentage_max TM_MAX`: Highest percentage of targeted mutations for the lowest temperature replica (default: 0.7).
- `-tm_perc_min, --target_mutations_percentage_min TM_MIN`: Lowest percentage of targeted mutations for the highest temperature replica (default: 0.0).
- `-seed, --seed IN_SEED`: User-defined seed number for simulation (default: 0).
- `-re_seq, --replicas_sequences {different,same}`: Choose whether replicas will start from the same or different random sequence (default: same).


### Example Usage

To run DesiRNA with a specific file and custom parameters:

```bash
DesiRNA.py -f example.txt -R 20 -e 50 -t 120
```

This command will run DesiRNA with the file `example.txt`, 20 replicas, an exchange frequency of 50, and a time limit of 120 seconds.

For further assistance with the command-line options, you can use the help command:

```bash
DesiRNA.py -h
```

## Files

### Input Files

The input file for DesiRNA must contain specific information related to the RNA sequence design. Here's an overview of the expected format and additional options:

#### Basic Single Chain RNA Design

- **Name Line:** `>name` followed by a unique identifier.
- **Sequence Restraints Line:** `>seq_restr` followed by sequence restraints using IUPAC nomenclature.
- **Secondary Structure Line:** `>sec_struct` followed by the secondary structure.

##### Example:

```
>name
Ete_11
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((.((((((((....))))).)).).))))))
```

#### Design with Pseudoknots

Include different brackets for pseudoknot representation.

**Note:** Available brackets for pseudoknots include `()`, `[]`, `<`, `>`, `{}`, `Aa`, `Bb`, `Cc`, `Dd`. Up to three levels of pseudoknots are accepted.

##### Example:

```
>name
Pseudoknot_Example
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((....[[[[..))))))......]]]]....
```

#### Two-Chain RNA Complex Design

Include an `&` symbol in both the sequence restraints and structure lines.

**Note:** A maximum of two-chain RNAs can be designed.


##### Example:

```
>name
Two_Chain_Example
>seq_restr
NNNNNNNNNNNNNNNNN&NNNNNNNNNNNNNNNNNN
>sec_struct
(((.(((((....))..&(((....)))..))))))
```

#### Design with Alternative Structures

Include an additional line for alternative structures.

**Note:** One can design a sequence folding into two structures maximum.

##### Example:

```
>name
Alt_Struct_Example
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((.((((((((....))))).)).).))))))
>alt_sec_struct
(((((((((((((....)))..)).)).).))))).
```

#### Design with Seed Sequence

Include a seed sequence for the starting point of the design simulation.

##### Example:

```
>name
Seed_Seq_Example
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((.((((((((....))))).)).).))))))
>seed_seq
GCCCCGGCCCCCGGCGAAAGCCGGUGGAGGCGGGGC
```

#### Sequence Restraints Dictionary:

The sequence restraints line uses IUPAC symbols such as `N` for any nucleotide, `W` for A or U, `S` for C or G, etc.

For more details on the IUPAC nomenclature, please refer to the [Wikipedia page on nucleic acid notation](https://en.wikipedia.org/wiki/Nucleic_acid_notation).




## Citation

If you use DesiRNA in your research, please cite our paper:

[DesiRNA: structure-based design of RNA sequences with a Monte Carlo approach](https://www.biorxiv.org/content/10.1101/2023.06.04.543636v1.abstract)

You can find the full citation details at the link above. Your citation helps support our work and ensures that it reaches a wider audience.

For any inquiries related to the paper or the use of DesiRNA, please feel fr
