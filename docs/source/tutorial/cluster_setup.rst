SLURM cluster setup
===================
Snakemake has the great functionality of being able to run just as well on your local computer as on a cluster. For the latter, it can also take care of submitting your jobs to the Slurm scheduler.

Profile
-------
A profile, is a folder which contains information on the default options you would like to execute snakemake with, as well as how to submit jobs. Luckily, there already exists a premade cookiecutter template for SLURM.

Firstly, you will need to connect to your cluster of choice and go to your working directory, install cookiecutter to your conda environment and then load the template. 

.. code-block:: console

    #connect to bioquant server
    ssh user@login1.bioquant.uni-heidelberg.de
    cd working-dir

    conda activate base
    conda install cookiecutter # in order to install profiles
    cookiecutter https://github.com/Snakemake-Profiles/slurm.git
    profile_name [slurm]: slurm
    sbatch_defaults []:
    cluster_config []:
    Select advanced_argument_conversion:
    1 - no
    2 - yes
    Choose from 1, 2 [1]: 1
    cluster_name []:

The cookiecutter command will interactively let you choose what kind of profile you would like to create. The profile_name will name a folder with your profile stored inside. I recommend storing this profile (folder) inside the ``config/`` directory.

The profile contains ``config.yaml`` that defines default options for snakemake when you run it with the profile. Here is what I am currently using on the BioQuant cluster.

.. code-block:: yaml

    restart-times: 0
    jobscript: "slurm-jobscript.sh"
    cluster: "slurm-submit.py"
    cluster-status: "slurm-status.py"
    cluster-cancel: "scancel"
    max-jobs-per-second: 1
    max-status-checks-per-second: 10
    local-cores: 16
    latency-wait: 60
    use-conda: True
    use-envmodules: True
    jobs: 1

You can also add default resources by adding it in this file, e.g.:

.. code-block:: yaml

    default-resources:
        - runtime=14400
        - mem_mb=16000

Renaming logs
-------------
When you run jobs on the cluster, all console output from snakemake will be redirected automatically to the slurm .out and .err files. You can change their naming by modifying ``settings.json`` to the following:

.. code-block:: json

    {
        "SBATCH_DEFAULTS": "job-name=smk-{rule}-{wildcards} output=logs/{rule}/{rule}-{wildcards}-%j.out",
        "CLUSTER_NAME": "",
        "CLUSTER_CONFIG": "",
        "ADVANCED_ARGUMENT_CONVERSION": "no"
    }

Change the naming scheme as is most convenient to you, or to a dedicated slurm logs folder


Rules with resources
--------------------
Rules can have ``resources`` and ``threads`` directives. For example:

.. code-block:: python

    rule run_misty_views:
        input:
            expr = "results/exp2/{sandwich}/Misty/view_p{radius}_expr.rds",
            metab = "results/exp2/{sandwich}/Misty/view_p{radius}_metab.rds"
        output: 
            directory("results/exp2/{sandwich}/Misty/model_p{radius}_int-{dtype}")
        threads: 4
        resources:
            mem_mb=25000,
            disk_mb=1000,
            time='12:00:00'
        script: "../scripts/R/misty_exp2_run.R"

These overwrite the defaults from the profile and ask for the appropriate resource allocation.

Using tmux
----------
When you execute snakemake on the cluster, it runs for the whole time your jobs are running as well and submits new jobs whenever necessary. However, if you just log in normally, any task running when you exit the ssh connection will be stopped.

tmux allows you to have a virtual session running, that you can connect to and disconnet from without cancelling any tasks. It is already preinstalled on the BioQuant cluster. The following command opens a new session called snakes:

.. code-block:: console

    tmux new -s snakes

You detach from it with ``Ctrl-B + D``. To re-attach to the previous session:

.. code-block:: console

    tmux a

Running snakemake in interactive mode
-------------------------------------
`Installing snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba>`_ and your conda environments is no different than installing it on your local computer.

Snakemake monitors your jobs while they are running. It should therefore **not be run on the login nodes** but in a separate interactive job.

.. code-block:: console

    tmux a #attach to your previous tmux session

    srun -t 5:00:00 --mem=5G --pty bash #slurm interactive job

    conda activate snk #activate your snakemake environment

    #check dry-run your rule all
    snakemake -n --profile ./path_profile_dir 

    #-j N specifies the max number of simultaneous jobs submitted at the same time
    #launch snakemake with max 10 parallel jobs
    snakemake -j 10 --profile ./path_profile_dir 

.. note::
    Be aware that when your interactive job ends, snakemake will not be able to track correct job completions etc. The interactive job should therefore have a longer max runtime than any of your jobs
