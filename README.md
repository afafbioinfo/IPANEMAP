# IPANEMAP: an Integrative Probing Analysis of Nucleic Acids Empowered by Multiple Accessibility Profiles.

A software for predicting stable RNA secondary structures compatibles from multiple chemical probing (SHAPE, DMS...)  reactivities profiles.
The output is an optimal centroid structure.

## Installing IPANEMAP

IPANEMAP consists in a set of Python 2.7+ scripts, and requires the prior installation, and accessibility from the command line, of the following dependencies:
1. `ViennaRNA` package, version posterior to 2.0
2. `scikit-learn`, version posterior to x.y
3. `VARNAv3-93.jar` should be found in the working repository

## Executing IPANEMAP

Once all dependencies are satisfied, IPANEMAP can be invoked through: 
  `python2.7 IPANEMAP.py`

## Configuration
All configuration options are set by modifying the content of the `IPANEMAP.Config` file

Parameters that should be set:
 - `constraint`: the names of reactivity profiles/experiments to be integrated in the analysis
 Example: `constraints: 1M7,1M7MG,NMIA`
 - `numberofsruct`: number of structures generated for each condition, once the value is fixed a folder `OutputSamples{numberofsruct}` is created

## Input folders
 - `Constraints_probing`: contains constraint files with reactivity scores. For each probing condition, two files are required: `RnaCondition.fa` and `RnaConditionProbing.txt` where `RnaCondition.fa` is a `Fasta` format file (sequence) with the identifier `>Rnacondition(reagent name)`
 - `fasta_files`: Fasta-format file of the studied RNA `Rna.fa`
 - `Constraints_Hard`: Contains files with hard constraints

## Output folders:

 - Multiprobing contains predicted structures (= optimal centroid structure) in dot-bracket format
 - output with all the clustering properties
 - logfile: Sampling, clustering and Pareto selection results

