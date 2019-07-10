# IPANEMAP
### Integrative Probing Analysis of Nucleic Acids Empowered by Multiple Accessibility Profiles.

IPANEMAP is a software for predicting stable RNA secondary structures compatible from multiple chemical probing (SHAPE, DMS...)  reactivities profiles. From one or several input sequences, along with several reactivity profiles, it computes and outputs one or several secondary structures, corresponding to the conformers best supported by experimental data and thermodynamics.

## Installing IPANEMAP

IPANEMAP consists in a set of Python 2.7+ scripts, and requires the prior installation, and accessibility from the command line, of the following dependencies:
1. `ViennaRNA` package, version posterior to 2.0, [downloadable from the TBI](https://www.tbi.univie.ac.at/RNA/#download "Download the Vienna package")
2. `scikit-learn`, [version 0.2 for Python 2.7+](https://scikit-learn.org/stable/install.html "Download scikit-learn")

## Executing IPANEMAP

Once all dependencies are satisfied, IPANEMAP can be invoked through: 

      python IPANEMAP.py [--conditions c1 c2 ...]
   
The method will run with a configuration specified within `IPANEMAP.cfg`, possibly overriding the list of considered conditions with the one provided through the optional command-line option `--conditions`.

## Configuration
Most configuration options are set by modifying the content of the `IPANEMAP.cfg` file.

### Main options
 - `RNA`: Specifies a path (relative to the working directory) to a FASTA file where the nucleotides sequence of the main RNA of interest can be found. Note that the filename is important, as it will be used as a base names for the other input files. Example: `RNA: fasta_files/didymium.fa` will process the sequence found in the file, and `didymium` will be used as the *base name* of reactivities/hard contraints files (see `Conditions` option).
 - `HardConstraintsDir` and `SoftConstraintsDir`: Those options specify the directories where (optional) hard and soft (reactivities) constraints, associated with the various conditions, can be found.
 - `Conditions`: Can be used to specify the list of probing conditions used for the prediction. Should be set to a comma-separated list of conditions, i.e. the names of reactivity profiles/experiments to be considered for the structure prediction. For a condition name `{Cond}`, the method will lead the method to  attempt to locate files named `{RNA}{Cond}.txt` in the `{HardConstraintsDir}` and `{SoftConstraintsDir}` folders, where `{RNA}` is the base name of the processed RNA.
Example: `Conditions: 1M7,1M7MG,NMIA` will jointly consider the three conditions. For `RNA: ` is set to  `RNA1M7.txt,`

 
### Input folders
 - `Constraints_probing`: contains constraint files with reactivity scores. For each probing condition, two files are required: `RnaCondition.fa` and `RnaConditionProbing.txt` where `RnaCondition.fa` is a `Fasta` format file (sequence) with the identifier `>Rnacondition(reagent name)`
 - `fasta_files`: Fasta-format file of the studied RNA `Rna.fa`
 - `Constraints_Hard`: Contains files with hard constraints

### Output folders
 - `Multiprobing` contains predicted structures (= optimal centroid structure) in dot-bracket format
 - `output` contains all the clustering properties
 - `logfile` lists all informations pertaining the sampling, clustering and Pareto selection processes

