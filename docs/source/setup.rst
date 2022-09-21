Setup a project environment
===========================

.. _cookiecutter:

Project structure
-----------------

The first step when starting a new project is to create a new working directory. Snakemake developers recommend a folder structure, which appropriately separates the workflow from data, results and configurations files.
Here is an example for this tutorial.

.. code-block:: none

    ├── .gitignore
    ├── workflow
    │   ├── rules
    |   │   ├── QC.smk
    │   ├── envs
    |   │   ├── pyQC.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   └── Snakefile
    ├── config
    │   ├── config.yaml
    ├── results
    │   ├── plots
    └── resources




Setting up snakemake
====================

Install the latest snakemake package as described in the `setup guide <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba>`_

For M1 MacOS
----------------
For many packages there are no versions available for the M1 chips. This occurs frequently, mainly with Bioconductor packages. If you encounter this situation, you can install packages that were compiled for Intel chips instead. For this you need to create the snakemake environment by setting appropriate flags. Snakemake will then also install packages that were designed for Intel.

.. code-block:: console
    
    CONDA_SUBDIR=osx-64 conda create -n snk64   # create a new environment called snk64
    conda activate snk64
    conda env config vars set CONDA_SUBDIR=osx-64  # subsequent commands use intel packages

You can then follow the same guide as above to install the snakemake package