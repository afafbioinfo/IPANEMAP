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
 - `SoftConstraintsDir` and `HardConstraintsDir`: Those options specify the directories where soft (reactivities) and hard constraints files can be found (if available).
 - `Conditions`: Can be used to specify the list of probing conditions used for the prediction. Should be set to a comma-separated list of conditions, i.e. the names of reactivity profiles/experiments to be considered for the structure prediction. 
 
For an RNA having base name `{RNA}`, and a condition name `{Cond}`, IPANEMAP will attempt to locate files named `{SoftConstraintsDir}/{RNA}{Cond}.txt` and `{HardConstraintsDir}/{RNA}{Cond}.txt`. If none of these files is found, the method will rely on a purely thermodynamic sampling.

Example: Given a configuration
 
      RNA: fasta_files/5sRNA.fa
      SoftConstraintsDir: soft
      HardConstraintsDir: hard
      Conditions: DMSMG,NMIA
   
the method will attempt to locate, and use for the sampling phase of the method, two files `5sRNADMSMG.txt` and `5sRNANMIA.txt` in each of the `soft` and `hard` directories.

### Paths-related options
 - `WorkingDir`: Main output directory for temp files and final results
 - `LogFile`: Name of file gathering the accumulated log

### Output folders
 - `Multiprobing` contains predicted structures (= optimal centroid structure) in dot-bracket format
 - `output` contains all the clustering properties
 - `logfile` lists all informations pertaining the sampling, clustering and Pareto selection processes

## How to
 - How do I perform a *pure thermodynamic*/constraints-free prediction? 
 Simply make sure that no constraint file named `{RNA}{Cond}.txt` is found in either `{SoftConstraintsDir}` or `{HardConstraintsDir}`, and IPANEMAP will default to a purely thermodynamic sampling (you may safely ignore the warning).
 - How do I specify a different sequence for some specific condition? 
If a FASTA file named `{RNA}{Cond}.fa` is found in either of the condition directories, then its sequence will be used instead of the main FASTA file. This can be useful when minor variants of the original sequence have been probed (eg Mutate-and-Map protocols).



