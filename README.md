# universal_chromatin_hmm

Universal chromatin HMMs (now superseeded by ChroMATIC)

# About

This is a repository of analysis scripts for the Hidden Markov Modelling of chromatin states, as used in Josserand et al 2023.  Chromatin state definitions and calls are based on unpublished data and research from the Marshall lab.

> **Warning** These scripts are forerunners of **ChroMATIC**, a complete gaussian Chromatin HMM analysis package that provides fully-automated analysis, machine-learning assignment of chromatin states, multi-modal inference and more.
>
> Please see our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.30.514435v1.abstract) for more details.
>
> The full source code of ChroMATIC will be released upon publication, but if you're interested in using ChroMATIC before this, please email the Marshall lab.

# Installation

*chromatin.universal.hmm.R* requires the following R packages:

-   tools

-   parallel

-   RHmm

-   ggplot2

-   ggcorrplot

-   stringr

-   ComplexHeatmap

# Usage

For Viterbi path fitting used for chromatin state predictions, R needs to be run with additional protected memory via the command-line argument `--max-ppsize=500000`.

```         
Rscript --max-ppsize=500000 ~/bin/chromatin.universal.hmms.R
```

### Required run-time parameters

A number of additional data files are optional to provide extra or different functionality, and can set via the following command-line options:

| Option            | Description                                                       |
|----------------|--------------------------------------------------------|
| `--states.matrix` | A matrix file for chromatin state prediction                      |
| `--genes.file`    | A GFF of all genes from the relevant genome release               |
| `--tfs`           | A list of transcription factors to determine enrichment by state. |

### Input data

The script will search for all files in bedGraph format present within the working directory, and generate models from these files.

### Other options

All command-line options and a description are listed below. Parameters are set with `--option=value` (no spaces around the equals or in parameter values).

To see these at runtime, along with any set default values, run the script with the `--help` option.
