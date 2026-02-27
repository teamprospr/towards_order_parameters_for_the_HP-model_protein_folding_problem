# Towards Order Parameters for the HP-Model Protein Folding Problem

This is an anonymized repository for the paper: Towards Order Parameters for the HP-Model Protein Folding Problem.
It offers the bare minimum code required for reproducing the required results and generating the figures used in the paper.

The code uses a protein folding framework called PROSPR.
We wrote this framework and made it open source available on GitHub under one of our public users.
This is an explicit note to **NOT** look up the framework for maintaining the anonymity.

## Directory structure

The directory tree looks as follows:

```console
.
├── code
│   ├── datasets
│   ├── experiments
│   ├── figures
│   ├── prospr
│   └── TC5b
├── data
│   ├── [REDACTED]_hratio
│   └── [REDACTED]_random
├── figures
└── jobs
    ├── hratio_l10_results
    ├── hratio_l15_results
    ├── hratio_l20_results
    ├── hratio_l25_results
    ├── hratio_l30_results
    ├── random_l10_results
    ├── random_l15_results
    ├── random_l20_results
    ├── random_l25_results
    ├── random_l30_results
    └── slurm
```

Within the `code/` directory you can find all the code required for:

1) Acquiring stats on the dataset and solved proteins (`code/datasets/`)
2) The protein folding experiments (`code/experiments/`)
3) Generating the figures used in the paper (`code/figures/`)
4) The PROSPR framework used for folding the proteins (`code/prospr/`)
5) Folding and visualizing the TC5b protein (`code/TC5b/`)

The datasets can be found under `/data/<dataset>/`, where the subdirectory is either the H-ratio or Random dataset.
We have redacted the name of the datasets here as it contains the sirname of one of the authors.
This is an explicit note to **NOT** look up the subdirectory names for maintaining the anonymity.

The `figures/` folder is a placeholder where all generated figures will be stored.

The `jobs/` folder contains all scripts for folding the proteins from the datasets within a multi-node SLURM-based environment.
There are some MPI distribution bugs within our main experiment `dfs_bnb_mpi`, and thus you might see that not all proteins from a dataset are folded.
That is why there is the `solve_remaining` code which takes a CSV as input for solving some specific proteins of a dataset.
The subdirectories (e.g. `jobs/hratio_l10_results/`) contain the results from our experiments, which are used by the code in the other directories.
The `jobs/slurm/` directory contains the slurm .out and .err files of our experiments.
These also contain some test executions, and thus do not fully correspond with soley folding the proteins.
Moreover, they contain references to our users and datasets, and thus you should **NOT** open these files to maintain the anonymity.

## Installing requirements

The required Python packages can be simply installed through pip:

```shell
pip install -r requirements.txt
```

or

```shell
python -m pip install -r requirements.txt
```

## Reproducing figures

Install the requirements and go into the `code/figures/` folder.
Here you'll find Python scripts for generating the figures, which will have names starting with `fig<num>_`.
Simply execute the script with Python and the figure should be generated in the `figures/` folder.
