# 2DesignRNA

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

