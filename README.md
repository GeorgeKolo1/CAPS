# CAPS

[![CI](https://github.com/GeorgeKolo1/CAPS/actions/workflows/python-package.yml/badge.svg)](https://github.com/GeorgeKolo1/CAPS/actions/workflows/python-package.yml)

## Overview

CAPS is a command-line tool that enables subtype-phenotype association testing. CAPS calculates uses one vs rest comparisons against subtypes and phenotypes to quantify the strength of association between subtypes and a given phenotype. CAPS does this in two primary ways:

1. Using an odds ratio that is calculated from a 2x2 contingency table
2.  If data is perfectly seperated, Firth logistic regression is used to avoid infinite estimates

CAPS also provides statistical testing to compare subtyping methods, including the discrimination index (Hunter & Gaston, 1988) and the adjusted Wallace coefficient (Severiano et al., 2011).

We acknowledge comparingpartitions.com as a predeccesor which is a platform that enables extensive comparison of subtyping methods using statistical tests (including the discrimination index and adjusted wallace coefficient)

## Quickstart

```caps -i path/to/your/subtyping_data -o path/to/output_directory -c comparator_column -p phenotype```

## Instructions

### Input requirements
#### Input

´´´-i or --i´´´

The input is required to be in either csv or tsv format. It should contain subtyping methods as columns and subtypes as values within the columns.
There is no limit on the number of subtypes (columns) you can include.

#### Comparator column

```-e or --e```

The comparator column should be the name of one of the subtyping columns within the input file. The comparator column is used to compare all subtyping methods to a single subtyping method in the adjusted wallace coefficient test. 

#### Phenotype

```-p or --p```

The phenotype file should be in comma seperated (.csv) or tab-seperated (.tab) format and include phenotypes 

### Output

The output file is generated in the user-specified directory. If the user provides a filename in the output for example ```-o path/to/your_files/results.tsv (or .csv)``` then the file will be named after the user defined output. If the user specifies only a directory, for example ```-o path/to/your_files/``` then the results file will be saved as results.csv in the user-specified directory.
