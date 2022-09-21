Cookiecutter project setup
==========================

.. note::

   This project is under active development. This might be updated later with an appropriate template for the saezlab.

In order to have a well-organised and clean workflow, it is recommended to use a structured working directory. Cookiecutter is very useful to setup a github repo structure from scratch in a few lines. First, check if it is installed

.. code-block:: console

    which cookiecutter

If you do not have it installed you can do so using conda

.. code-block:: console
    
    conda install cookiecutter

Once this is done, navigate to a folder where you would like to create your new working directory of your project. The command below will create a new directory containing most files necessary for a new project. Cookiecutter will prompt you for some 

.. _cookiecutter:
.. code-block:: console

    cookiecutter gh:snakemake-workflows/cookiecutter-snakemake-workflow

    full_name [Johannes Köster]: Cameron Doe
    email [johannes.koester@protonmail.com]: cameron.doe@uni-heidelberg.de
    username [johanneskoester]: cameron.d
    project_name [RNA-Seq]: tutorial
    repo_name [rna-seq]: snk-tutorial
    min_snakemake_version [5.7.0]: 6.0

This template creates several files which are used for standardised workflows, including validation steps as well as dockers singularities. You can just ignore this for now, and check up on them once you feel you need these particular features. I usually do not go that much into details, but for bigger projects - that you might publish - it might be useful to start including these validations from the start.

.. note::
    Please be aware that you can already find some workflows for various bioinformatic pipelines from the `snakemake workflows project <https://github.com/snakemake-workflows/docs>`_