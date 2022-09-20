Setup a project environment
===========================

.. _cookiecutter:

Project structure
-----------------

The first step when starting a new project is to create a new working directory. Snakemake developers recommend the folder structure shown below, which appropriately separates the workflow from data, results and configurations files.

.. code-block:: none

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── workflow
    │   ├── rules
    |   │   ├── module1.smk
    |   │   └── module2.smk
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── notebooks
    |   │   ├── notebook1.py.ipynb
    |   │   └── notebook2.r.ipynb
    │   ├── report
    |   │   ├── plot1.rst
    |   │   └── plot2.rst
    |   └── Snakefile
    ├── config
    │   ├── config.yaml
    │   └── some-sheet.tsv
    ├── results
    └── resources

While one can just set it up manually, we can also use cookiecutter in order to make our life easier. First, check if you have cookiecutter installed.

.. code-block:: console

    which cookiecutter

If you do not have it installed you can do so using conda

.. code-block:: console
    
    conda install cookiecutter

Once this is done, navigate to a folder where you would like to create your new working directory of your project. The command below will create a new directory containing most files necessary for a new project. Cookiecutter will prompt you for some 

.. code-block:: console

    cookiecutter gh:snakemake-workflows/cookiecutter-snakemake-workflow
    full_name [Johannes Köster]: Cameron Doe
    email [johannes.koester@protonmail.com]: cameron.doe@uni-heidelberg.de
    username [johanneskoester]: cameron.d
    project_name [RNA-Seq]: tutorial
    repo_name [rna-seq]: snk-tutorial
    min_snakemake_version [5.7.0]: 6.0

Please be aware that you can already find some workflows for various bioinformatic pipelines from the `snakemake workflows project <https://github.com/snakemake-workflows/docs>`_

.. note::
    This template creates several files which are used for standardised workflows, including validation steps as well as dockers singularities. You can just ignore this for now, and check up on them once you feel you need these particular features. I usually do not go that much into details, but for bigger projects - that you might publish - it might be useful to start including these validations from the start.

Setting up snakemake
--------------------

Install the latest snakemake package as described in the `setup guide https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba`_

**For M1 MacOS**
For many packages there are no versions available for the M1 chips. This occurs frequently, mainly with Bioconductor packages. If you encounter this situation, you can install packages that were compiled for Intel chips instead. For this you need to create the snakemake environment by setting appropriate flags. Snakemake will then also install packages that were designed for Intel.

.. code-block::
    CONDA_SUBDIR=osx-64 conda create -n snk64   # create a new environment called snk64
    conda activate snk64
    conda env config vars set CONDA_SUBDIR=osx-64  # subsequent commands use intel packages

You can then follow the same guide as above to install the snakemake package