Advanced rule design
====================

Global wildcards
----------------
If you look back at the ``merge`` rule in the ``QC.smk`` module, you can see that we specifically input three samples in the ``expand()`` function.

.. code-block:: python

    rule merge:
        input:
            expand('results/QC/filtered_{sample}.h5ad', sample = ['sample1', 'sample2', 'sample3'])
        output:
            data = 'results/merged.h5ad',
            plot = 'plots/umap.pdf'
        conda:
            '../envs/scanpy.yaml'
        script:
            '../scripts/merge.py'

However, if we added additional samples to the directory and run the pipeline by calling for the creation of ``merged.h5ad``, it would only include data for these three samples.

To remedy this, we can dynamically create a list of samples from the contents of the ``data/`` folder using the ``glob_wildcards()`` function. Add the following lines before ``rule all`` in the ``Snakefile``

.. code:: python

    SAMPLES = glob_wildcards('data/{sample}.h5ad')
    print(SAMPLES)

.. note:: 
    Global wildcards are often written in upper case

If you now do a dry-run, you will see which files correspond to this naming scheme.

.. code-block:: console

    snakemake -n

.. code-block:: console

    Wildcards(sample=['sample1', 'sample2', 'sample3'])

You can see that the function therefore returns a dictionary, and our sample names are under the ``sample`` key. Let's change the Snakefile to the following:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("7.0")

    configfile: "config/config.yaml"

    SAMPLES = glob_wildcards('data/{sample}.h5ad').sample
    print(SAMPLES)

    # define which output files you would like to generate
    rule all:
        input:
            'data/sample1.h5ad',
            expand("results/QC/filtered_{sample}.h5ad", sample = SAMPLES)


    module download_data:
        snakefile: "rules/download_data.smk"
        config: config

    use rule * from download_data as dwn_*

    module QC:
        snakefile: "rules/QC.smk"
        config: config

    use rule * from QC as QC_*

You can see that we are calling for the creation of the filtered objects. Thus, if you just type ``snakemake -n`` you will see which files would be created.

.. note:: 
    Remember that ``rule all`` specifies which files to create, when no files are given to the snakemake command.


Checkpoints and input functions
-------------------------------

In the above example, it assumes that the samples are already present in the ``data/`` folder when you start with the execution. If you were to delete them, of course snakemake would throw an error.

In this section, we will take a look at how to handle steps in pipelines where it is not possible to know the names of output files before they are executed. This is often the case for scripts or command line commands that return folders.

Take a look at the download module (``workflow/rules/download.smk``):

.. code-block:: python

    import os

    checkpoint download:
        output:
            directory('data/filtered_gene_bc_matrices/hg19')
        shell:
            '(test -d data || mkdir data)'
            '&& wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz '
            '&& cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz'

    def get_downloaded_file(wildcards):
        checkpoint_output = checkpoints.download.output[0]
        return os.path.join(checkpoint_output, "matrix.mtx")

    rule make_samples:
        input:
            get_downloaded_file
            # 'data/filtered_gene_bc_matrices/hg19/matrix.mtx'
        output:
            'data/sample1.h5ad',
            'data/sample2.h5ad',
            'data/sample3.h5ad'
        conda:
            '../envs/scanpy.yaml'
        script:
            '../scripts/fake_samples.py'

The first block is a special rule, called a checkpoint, that has a folder as output flagged with ``directory()``. It downloads and extracts the single cell data. 

The second block is a so-called input function that only has ``wildcard`` as argument. This is necessary in order to propagate the wildcards from output to input (remember that snakemake is a backwards-looking worflow manager), and tell snakemake which rules' inputs/outputs are connected.

`Input functions<https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#input-functions>`_ have to return a list of filenames, or alternatively a dict (e.g. with named keys), that you need to `unpack. <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#input-functions-and-unpack>`_

The last block is a rule that just splits the data into three, and creates toy samples. The input directive is only the input function defined above, but other files could have been added as well.

In this example, if we are asking snakemake to create a sample file (e.g. ``data/sample3.h5ad``), snakemake knows that the download checkpoint needs to be run, but does not yet know the filenames given as input to ``rule make_samples``. When it is done downloading and extracting the data, it will run the input function and then pass its output to the last rule.

You can see this in the rule blocks in the output of the following command:

.. code-block:: console

    snakemake -n --forceall data/sample1.h5ad

.. code-block:: console

    Building DAG of jobs...
    Job stats:
    job                 count    min threads    max threads
    ----------------  -------  -------------  -------------
    dwn_download            1              1              1
    dwn_make_samples        1              1              1
    total                   2              1              1


    [Mon Sep 26 11:34:38 2022]
    checkpoint dwn_download:
        output: data/filtered_gene_bc_matrices/hg19
        jobid: 1
        resources: tmpdir=/var/folders/vl/1y1qg3c911x2hvqbsl7zfpz40000gn/T
    Downstream jobs will be updated after completion.


    [Mon Sep 26 11:34:38 2022]
    rule dwn_make_samples:
        input: <TBD>
        output: data/sample1.h5ad, data/sample2.h5ad, data/sample3.h5ad
        jobid: 0
        resources: tmpdir=/var/folders/vl/1y1qg3c911x2hvqbsl7zfpz40000gn/T

    Job stats:
    job                 count    min threads    max threads
    ----------------  -------  -------------  -------------
    dwn_download            1              1              1
    dwn_make_samples        1              1              1
    total                   2              1              1

    This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

lambda expresssions
-------------------

In case you want to feed parameters (or file names) based on wildcards present in rules, you can also do so using lambda expressions. 

Let's look at an example where you set QC parameters dynamically based on the sample name. 

Before you start, just execute this command to make sure your files are up to date with your current workflow:

.. code-block:: console

    snakemake results/merged.h5ad --forceall --use-conda -c3

You will first change the ``config.yaml`` file and make entries for each sample, e.g.:

.. code-block:: yaml

    QC:
        'min_gene': 200
        'min_cells': 3
        'max_pct_mt': 5

    sample1:
        'min_gene': 100
        'min_cells': 3
        'max_pct_mt': 10

    sample2:
        'min_gene': 300
        'min_cells': 5
        'max_pct_mt': 5

    sample3:
        'min_gene': 200
        'min_cells': 10
        'max_pct_mt': 5

.. note:: 
    To be clear, it is not recommende to use different QC parameters for different samples, but this is just showing off some of the capabilities of snakemake based on this example.

You now need to change your ``QC_sample`` rule, and pass these parameters instead:

.. code-block:: python
    
    rule QC_sample:
        input:
            'data/{sample}.h5ad'
        output:
            temp('results/QC/filtered_{sample}.h5ad')
        params:
            lambda w: config[w.sample]
        conda:
            '../envs/scanpy.yaml'
        script:
            '../scripts/QC_samples.py'

You can see that we are passing only the part of the config file that is corresponding to the sample wildcard. If this entry is missing for one sample, it will throw an error, but not if any of the subkeys is missing. 

.. note:: 
    Instead of using the config file, you can also load the parameters from a panda dataframe. However, changes to this file will not be tracked.

Check what happens when you run the following:

.. code-block:: console

    snakemake results/merged.h5ad --use-conda -n

Snakemake should tell you that the file already exists, but that some parameters have changed:

.. code-block:: console

    Building DAG of jobs...
    Updating job dwn_make_samples.
    The params used to generate one or several output files has changed:
        To inspect which output files have changes, run 'snakemake --list-params-changes'.
        To trigger a re-run, use 'snakemake -R $(snakemake --list-params-changes)'.
    Nothing to be done (all requested files are present and up to date).

It can also tell you which parameters have changed, with this command:

.. code-block:: console

    snakemake --list-params-changes results/merged.h5ad

You will need to force execution in order to get it to run:


.. code-block:: console
    
    snakemake results/merged.h5ad --use-conda --force -c3