Setup project environment
=========================

Setting up snakemake
====================
Of course, any decent snakemake tutorial will start by installing snakemake. Install the latest snakemake package as described in the `setup guide <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba>`_

For M1 MacOS
----------------
For many packages there are no versions available for the M1 chips. If you encounter this situation, you can install packages that were compiled for Intel chips instead. For this, you need to create the snakemake environment by setting appropriate flags. Snakemake will then automatically install packages that were designed for Intel. This occurs frequently, mainly with Bioconductor packages.

.. note:: 
    Be aware that the first time you load/execute packages from this conda environment it might seems like the execution is hanging. This is normal as the code is translated for your processor and subsequent executions will be faster.

.. code-block:: console
    
    CONDA_SUBDIR=osx-64 conda create -n snk64   # create a new environment called snk64
    conda activate snk64
    conda env config vars set CONDA_SUBDIR=osx-64  # subsequent commands use intel packages

Conda might prompt you to activate the base environment again to save changes. It will also display warning messages every time you (de-)activate this environment.

You can now follow the same guide as above to install the snakemake package.


Project structure
=================
Now that you have a conda environment with basically only snakemake installed, let's clone the tutorial git repository. Navigate to a folder where you would like to clone it into.

.. code-block:: console
    git clone https://github.com/saezlab/snk-tutorial

The project repository structure is shown below. Snakemake developers recommend a folder structure, which separates the workflow from data, results and configurations files. 

.. code-block:: none

    ├── .gitignore
    ├── workflow
    │   ├── rules
    |   │   ├── download_data.smk
    │   ├── envs
    |   │   ├── scanpy.yaml
    │   ├── scripts
    |   │   ├── fake_samples.py
    |   │   ├── QC_samples.py
    |   └── Snakefile
    └── config
        └── config.yaml

This structure is a more lightweight version of the cookiecutter version provided by the snakemake workflows project. You can find how to setup your own project directory :ref:`here <cookiecutter_setup>`