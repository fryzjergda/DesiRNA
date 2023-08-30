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
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Citation](#citation)
- [Benchmark](#benchmark)
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

- `-R, --replicas`: Number of replicas (default: 10).
- `-e, --exchange`: Frequency of replica exchange attempt (default: 100).
- `-t, --timelimit`: Time limit for running the program in seconds (default: 60).
- `-s, --steps`: Number of Replica Exchange steps, overwrites the `-t` option (default: None).
- `-r, --results_number`: Number of best results to be reported (default: 10).
- `-p, --param {2004,1999}`: Turner energy parameter for calculating MFE (default: 1999).
- `-tmin, --tmin`: Minimal Replica Temperature (default: 10).
- `-tmax, --tmax`: Maximal Replica Temperature (default: 150).
- `-ts, --tshelves`: Custom temperature shelves for replicas, provide comma-separated values.
- `-sf, --scoring_function`: Scoring functions and weights used to guide the design process. Multiple scoring functions can be selected. Please provide desired scoring functions and their weigths e.g., ```-sf Ed-Epf:0.5,1-MCC:0.5```. (default: `Ed-Epf:1.0`). Available options:
  - `Ed-Epf`: Energy of the desired structure (Ed) minus free energy of the thermodynamic ensemble (Epf).
  - `Ed-MFE`: Energy of the desired structure (Ed) minus Minimum Free Energy (MFE)
  - `1-MCC`: One minus Matthews Correlation Coefficient (MCC).
  - `sln_Epf`: Sequence Length Normalized Epf.
  - `1-precision`: One minus precision (TP/(TP+FP)).
  - `1-recall`: One minus recall (TP/(TP+FN)).

- `-nd, --negative_design {off,on}`: Use negative design approach (default: off).
- `-acgu, --ACGU {off,on}`: Keep 'natural' ACGU content, with default content A:15%, C:30%, G:30%, U:15% (default: off).
- `-acgu_content, --ACGU_content`: Provide user-defined ACGU content, comma-separated values e.g., ```-acgu_content 10,40,40,10```.
- `-o, --oligomerization {off,enforce,avoid}`: Check if the designed sequence tends to oligomerize. User may enforce or avoid oligomerization. Slows down the simulation (default: `off`).
- `-d, --dimer {off,on}`: Design of a homodimer complex, of two strands. Requires input file complying with RNA-RNA complex format (default: `off`).
- `-tm, --target_mutations {off,on}`: Targeted mutations (default: on).
- `-tm_perc_max, --target_mutations_percentage_max`: Highest percentage of targeted mutations for the lowest temperature replica (default: 0.7).
- `-tm_perc_min, --target_mutations_percentage_min`: Lowest percentage of targeted mutations for the highest temperature replica (default: 0.0).
- `-motifs, --motif_sequences`: Prevent or enforce specific sequence moitif. Provide sequence motifs along with their bonuses(-)/penalties(+), e.g., ```-motifs GNRA,-1,CCCC,2```.
- `-seed, --seed_number`: User-defined seed number for simulation (default: 0).
- `-re_seq, --replicas_sequences {different,same}`: Choose whether replicas will start from the same or different random sequence (default: same).


### Example Usage

To run DesiRNA with a specific file and custom parameters:

```bash
DesiRNA.py -f Standard_design_input.txt -R 20 -e 50 -t 120
```

This command will run DesiRNA with the file `Standard_design_input.txt`, 20 replicas, an exchange frequency of 50, and a time limit of 120 seconds.

For further assistance with the command-line options, you can use the help command:

```bash
DesiRNA.py -h
```


## Input Files

The input file for DesiRNA must contain specific information related to the RNA sequence design. Here's an overview of the expected format and additional options:

### Basic Single Chain RNA Design

- **Name Line:** `>name` followed by a unique identifier.
- **Sequence Restraints Line:** `>seq_restr` followed by sequence restraints using IUPAC nomenclature.
- **Secondary Structure Line:** `>sec_struct` followed by the secondary structure.

#### Example:

```
>name
Design
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((.((((((((....))))).)).).))))))
```

####

### Design with Pseudoknots

Include different brackets for pseudoknot representation.

**Note:** Available brackets for pseudoknots include `()`, `[]`, `<`, `>`, `{}`, `Aa`, `Bb`, `Cc`, `Dd`. Up to three levels of pseudoknots are accepted.

#### Example:

```
>name
Pseudoknot_Example
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((....[[[[..))))))......]]]]....
```

### Two-Chain RNA Complex Design

Include an `&` symbol in both the sequence restraints and structure lines.

**Note:** A maximum of two-chain RNAs can be designed.


#### Example:

```
>name
RNA-RNA Complex_Example
>seq_restr
NNNNNNNNNNNNNNNNN&NNNNNNNNNNNNNNNNNN
>sec_struct
(((.(((((....))..&(((....)))..))))))
```

### Two-Chain Homodimer Design

Include an `&` symbol in both the sequence restraints and structure lines. The two sequences, must be of the same length.

**Note:** Turn on the homodimer option `-d on`.


#### Example:

```
>name
Homodimer_Example
>seq_restr
NNNNNNNNNNNNNNNNN&NNNNNNNNNNNNNNNNN
>sec_struct
((((....((((.....&))))....)))).....

```



### Design with Alternative Structures

Include additional lines for alternative structures.

#### Example:

```
>name
Alt_Struct_Example
>seq_restr
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>sec_struct
((((((.((((((((....))))).)).).))))))
>alt_sec_struct
(((((((((((((....)))..)).)).).))))).
(((((((((((((....)))))...)).).))))).
```

### Design with Seed Sequence

Include a seed sequence for the starting point of the design simulation.

#### Example:

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


## Output Files

The results of a simulation are organized into a specific directory, named according to the parameters and options used in the simulation. The directory includes the following files and subdirectory:

### Main Directory

- `*_mid_results.csv`: A CSV file containing intermediate results of the simulation.
- `*_replicas.png`: A PNG image showing the design simulation process.
- `*_results.csv`: A CSV file containing the final results of the simulation.
- `*_stats`: A file containing statistical information about the simulation.
- `*.txt`: The original input file used for the simulation.

#### Trajectory Files Subdirectory

Within the main directory, there is a `trajectory_files` folder containing additional files related to the simulation trajectory:

- `*_best_fasta.fas`: A FASTA file containing the best sequence(s) from the simulation.
- `*_best_str`: A file containing the best solution.
- `*_command`: A file containing the command used to run the simulation.
- `*_multifasta.fas`: A FASTA file containing all sequences from the simulation.
- `*_random.csv`: A CSV file containing random sequences fitting compying to the desired structure.
- `*_replicas.csv`: A CSV file containing information about the varius stats for each replica used in the simulation.
- `*_traj.csv`: A CSV file containing the trajectory of the simulation.


## Benchmark

The benchmark files are organized into two main directories: `Eterna100V1_benchmark_results` and `Eterna100V1_inputs`.

### Eterna100V1_benchmark_results

This directory contains the results of the benchmark tests. The files include:

- `Eterna100V1_1h_results.txt`: Results for the 1-hour benchmark test.
- `Eterna100V1_1min_results.txt`: Results for the 1-minute benchmark test.
- `Eterna100V1_24h_results.txt`: Results for the 24-hour benchmark test.
- `Eterna100V1_all_results.txt`: Consolidated results for all benchmark tests.

### Eterna100V1_inputs

This directory contains the input files used for the benchmark tests. There are 100 individual text files (`eteV1_01.txt`, `eteV1_02.txt`, ..., `eteV1_100.txt`) that contain the structures of 100 Eterna benchmark puzzles.

These files are organized in the same format as described in the [Input Files](#input-files) section.



## Citation

If you use DesiRNA in your research, please cite our paper:

[DesiRNA: structure-based design of RNA sequences with a Monte Carlo approach](https://www.biorxiv.org/content/10.1101/2023.06.04.543636v1.abstract)

You can find the full citation details at the link above. Your citation helps support our work and ensures that it reaches a wider audience.

For any inquiries related to the paper or the use of DesiRNA, please feel fr


## License

DesiRNA is licensed under the Apache License, Version 2.0. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

For the full license text, please see the [LICENSE](LICENSE) file in the repository.

## Contact

If you have any questions, need support, or would like to collaborate, please feel free to reach out to us:

- **Dr. Tomasz. K Wirecki**: [twirecki@iimcb.gov.pl](mailto:twirecki@iimcb.gov.pl)
- **Prof. Janusz M. Bujnicki**: [janusz@iimcb.gov.pl](mailto:janusz@iimcb.gov.pl)

You can also learn more about our research and other projects on our lab's website: [GeneSilico Lab](https://genesilico.pl/)

We welcome feedback and contributions to DesiRNA, and we look forward to hearing from you!
