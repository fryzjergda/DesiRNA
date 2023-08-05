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

### Without Conda

```
python3 -m pip install -r requirements.txt
python3 -m pip install viennarna
```

### With Conda - creating your own environment

```
conda create -n desirna python=3 numpy matplotlib pandas shutil pathlib
conda activate desirna
pip install viennarna
```

### With Conda - using Conda environment YAML file

```
conda env create -f DesiRNA-env.yml
conda activate DesiRNA-env
```


Certainly! Here's a more formal version of the "Installation" section:

```markdown
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

### Adding DesiRNA to PATH

To facilitate access to DesiRNA from the command line, the path to `DesiRNA.py` may be added to the `.bashrc` file:

```bash
echo 'export PATH=$PATH:/path/to/DesiRNA.py' >> ~/.bashrc
source ~/.bashrc
```

Replace `/path/to/DesiRNA.py` with the correct path to the `DesiRNA.py` file.

### Additional Installation Guidance

- **Conda Installation:** Should Conda be required, the [official Conda installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) provides comprehensive instructions.
- **ViennaRNA Installation:** ViennaRNA may be installed through various methods. Detailed instructions are available in the [ViennaRNA installation guide](https://www.tbi.univie.ac.at/RNA/#download).

It is imperative to ensure that all dependencies are properly installed prior to utilizing DesiRNA.
```







## Usage

You can execute the `DesiRNA.py` program using the following command:

```bash
python DesiRNA.py [arguments]
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
python DesiRNA.py -f yourfile -R 20
```

Replace `yourfile` with the path to your input file and `20` with the number of replicas you want to use. If your file path contains spaces, please enclose it in quotation marks.


