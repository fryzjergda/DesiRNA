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
- `-p, --param {2004,1999}`: Turner energy parameter for calculating MFE (default: 1999).
- `-tmin, --tmin T_MIN`: Minimal Replica Temperature (default: 10).
- `-tmax, --tmax T_MAX`: Maximal Replica Temperature (default: 150).
- `-ts, --tshelves TSHELVES`: Custom temperature shelves for replicas, provide comma-separated values.
- `-sf, --scoring_function {dmt,mcc,mfe,mix,mix2,alt}`: Scoring function used to guide the design process (default: dmt).
- `-nd, --negative_design {off,on}`: Use negative design approach (default: off).
- `-acgu, --ACGU {off,on}`: Keep 'natural' ACGU content (default: off).
- `-pk, --PK {off,on}`: Design of pseudoknotted structures (default: off).
- `-o, --oligomerization {off,on}`: Check if the designed sequence tends to oligomerize (default: off).
- `-d, --dimer {off,on}`: Design of an RNA complex, of two strands (default: off).
- `-tm, --target_mutations {off,on}`: Targeted mutations (default: on).
- `-a, --alt_ss {off,on}`: Design of sequences folding into two structures (default: off).
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

The input file for DesiRNA must contain information related to the RNA sequence design, specifically for single chain RNA. Here's an overview of the expected format:

- **Name Line:** Starts with `>name` on one line, followed by a unique identifier for the RNA design on the next line (e.g., `Ete_11`).
- **Sequence Restraints Line:** Starts with `>seq_restr` on one line, followed by the sequence restraints on the next line, where `N` represents any nucleotide. The sequence restraints define the allowed nucleotides at each position.
- **Secondary Structure Line:** Starts with `>sec_struct` on one line, followed by the secondary structure in dot-bracket notation on the next line.

#### Example:

A basic single chain RNA design file may look like this:


```
>name
Ete_11
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((.((((((((....))))).)).).))))))
```

