.. _qc tutorial:
QC samples
==========

In order to illustrate the basic functionalities in snakemake, this tutorial uses the pbmc dataset from 10x Genomics. After running the :ref:`setup tutorial <setup for tutorial>`, you should have three samples in compressed Anndata format in the ``data`` folder.

You will first learn how to implement a snakemake ``rule`` that does basic quality control filtering of a sample, and learn how to generalise it for any sample using so-called ``wildcards``. Then you will learn how to aggregate several files/samples using the ``expand()`` function.

.. note:: 
    This tutorial is supposed to teach how to use snakemake in order to develop reproducible pipelines. Even though it uses single cell data analysis as a way to illustrate snakemake's functionalities, it should not be taken as a state-of-the-art template of single cell expression analysis, for that we recommend looking at the `scanpy tutorials <https://scanpy.readthedocs.io/en/latest/tutorials.html>`_.

Simple rule
-----------
In Snakemake, a ``rule`` is an instruction - written in YAML - that ties a (set of) input file(s) to output files/directories, and an associated operation, which can either be an inline chunk of python code, a script (python/R), a notebook or any set of operation that can be written out in a shell command.

The rule below specifies that in order to create a filtered version of the data from sample1, it requires one file as input, and uses a python script, that should be run in a given conda environment. You can take a look at the script in order to see what it does exactly, but it does some basic filtering and normalises the count data.

.. note:: 
    The paths for input and output files are relative to the working directory, whereas the paths for the conda environment and the script are relative to the location of the ``Snakefile`` (from the file location where the rule is defined).

.. code-block:: python

    rule QC_sample:
        input:
            'data/sample1.h5ad'
        output:
            'results/QC/filtered_sample1.h5ad'
        conda:
            'envs/scanpy.yaml'
        script:
            'scripts/QC_samples.py'


.. note:: 
    Each directive (input, output, etc) is available in python as a list: e.g. the input path as ``snakemake.input[0]``

You can add it at the bottom of the Snakefile and run the command below in order to produce the filtered data. The command specifies which file to create, that it should load a specific environment, and that you want up to use up to one core.

.. code-block:: console

    snakemake results/QC/filtered_sample1.h5ad --use-conda -c1

The output shows you some snakemake information about the jobs that will be run, the activation of the conda environment and then the some script output. Notice that your working directory now contains a ``results`` folder with the data, as snakemake takes care of folder and file management after they are encoded in the rules. It also takes care of cleaning up files produced by failed jobs, since those might be corrupted in some way.

Wildcards
---------
The rule shown above takes care of the quality control filtering for one sample, but we have three. So-called ``wildcards`` are a kind of placeholder string that is only determined at execution time and snakemake finds out which operations it needs to do.

Modify the previous rule to the following:

.. code-block:: python

    rule QC_sample:
        input:
            'data/{sample}.h5ad'
        output:
            'results/QC/filtered_{sample}.h5ad'
        conda:
            'envs/scanpy.yaml'
        script:
            'scripts/QC_samples.py'

This rule will now be able to QC any files that have the right filepath pattern. You can rerun the previous command for the second sample:

.. code-block:: console

    snakemake results/QC/filtered_sample2.h5ad --use-conda -c1

Add parameters
--------------
The ``'scripts/QC_samples.py'`` script is throwing a warning because it expects filtering parameters and is currently using hard-coded defaults. Modify your rule so that you can define your own parameters for your samples:

.. code-block:: python

    rule QC_sample:
        input:
            'data/{sample}.h5ad'
        output:
            'results/QC/filtered_{sample}.h5ad'
        params:
            min_gene = 200,
            min_cells = 3,
            max_pct_mt = 5
        conda:
            'envs/scanpy.yaml'
        script:
            'scripts/QC_samples.py'

Notice that the lines are separated with a comma and have a name (mainly for readability).

.. note:: 
    In python, the names of the parameters do not matter and will be available in a list, i.e. ``snakemake.params[0]``, ``snakemake.params[1]``, etc.

    In R, named elements are duplicated and available either by index or by name, i.e. here ``snakemake@params`` has length **6** and elements are accessible e.g. as ``snakemake@params$min_gene`` and ``snakemake@params[[0]]``

For rules/pipelines with many parameters, it can be quite a hassle to parse all of these parameters and keep track of where you need to change them. Instead of adding each individually, you can pass specific keys of the ``config`` file that contains these parameters:

.. code-block:: yaml

    #contents of the yaml file
    project: 'snk-tutorial'

    QC:
        'min_gene': 200
        'min_cells': 3
        'max_pct_mt': 5


.. code-block:: python

    rule QC_sample:
        input:
            'data/{sample}.h5ad'
        output:
            'results/QC/filtered_{sample}.h5ad'
        params:
            config['QC']
        conda:
            'envs/scanpy.yaml'
        script:
            'scripts/QC_samples.py'

There are two advantages of using this approach: firstly, it simplifies and centralises parameter management to one single file, and secondly, the changes of parameters are also tracked by snakemake and you will be prompted to rerun the pipeline if they do.

.. note:: 
    These two rule examples are equivalent in what they do, but the parsing of the parameters is different. The latter actually passes a whole python ``dict`` object in ``snakemake.params[0]``. Check out ``'scripts/QC_samples.py'`` if you want to know more.

Merge
=====
Now that we have QCed and normalised all the files, we can proceed with combining them to do clustering and create some nice-looking UMAPs. 

As with parameters, a rule can take more than one file as input and output. You could write them out one-by-one, however that can be inconvenient when pooling many samples at a time. You can simplify this by using the ``expand()`` function as shown in the example below.

.. code-block:: python

    rule merge:
        input:
            expand('results/QC/filtered_{sample}.h5ad', sample = ['sample1', 'sample2', 'sample3'])
        output:
            data = 'results/merged.h5ad',
            plot = 'plots/umap.pdf'
        conda:
            'envs/scanpy.yaml'
        script:
            'scripts/merge.py'

You can add this new rule at the bottom of the ``Snakefile`` and check what would be executed by doing a dry-run again. Note that you need to specify only one of the output files, the others will be produced anyway.

.. code-block:: console

    snakemake results/merged.h5ad -n

The output of snakemake should show you a table similar to the following:

.. code-block:: console

    Job stats:
    job          count    min threads    max threads
    ---------  -------  -------------  -------------
    QC_sample        1              1              1
    merge            1              1              1
    total            2              1              1

You can see that it will QC the remaining sample, and then do the merge. Since you have executed the QC for the two first samples, and neither the input data nor the script having been changed in the meantime, these don't need to be run again.

Temporary files
---------------
If we run the pipeline as is, we will retain the raw data *and* the filtered data of each sample, as well as the merged data. If you have many more samples, this might represent a sizeable amount of duplicate data. Therefore, we can flag the filtered data as being only temporary using ``temp()``. This way, snakemake will delete the files once all jobs requiring them have been executed successfully.

.. note:: 
    Be aware that it might not be wise to flag data which takes complex computations as temporary, at least during pipeline development.

The two rules in your ``Snakefile`` should now look this way:

.. code-block:: python

    rule QC_sample:
        input:
            'data/{sample}.h5ad'
        output:
            temp('results/QC/filtered_{sample}.h5ad')
        params:
            config['QC']
        conda:
            'envs/scanpy.yaml'
        script:
            'scripts/QC_samples.py'

    rule merge:
        input:
            expand('results/QC/filtered_{sample}.h5ad', sample = ['sample1', 'sample2', 'sample3'])
        output:
            data = 'results/merged.h5ad',
            plot = 'plots/umap.pdf'
        conda:
            'envs/scanpy.yaml'
        script:
            'scripts/merge.py'

You can now execute the last step in the pipeline by running the command:

.. code-block:: console

    snakemake results/merged.h5ad --use-conda -c1

Snakmeake should tell you that the temporary files were deleted upon successfull job completion

.. code-block:: console

    Removing temporary output results/QC/filtered_sample1.h5ad.
    Removing temporary output results/QC/filtered_sample2.h5ad.
    Removing temporary output results/QC/filtered_sample3.h5ad.

.. note:: 
    If a job depending on temporary files fails, the temp files are not deleted.

    For debugging your workflow, it might still be usefull to use the ``--notemp`` option, which prevents temporary file deletion, especially if they a long time to compute.

In your ``results`` directory, you should now indeed have only the ``merged.h5ad`` file. You can also look at your nice UMAPs in the ``plots`` folder.

.. note:: 
    In general, I have found it usefull to clearly separate rules that do data cleaning, transformation, etc from plotting, as the latter can change quite often in order to make plots more readable. It is therefore not very practical if the plots are tied to expensive computations.

Force execution
---------------
If you now try to execute the previous command again, snakemake will tell you that there is nothing to be done:

.. code-block:: console

    snakemake results/merged.h5ad --use-conda -c1

.. code-block:: console

    Building DAG of jobs...
    Updating job dwn_make_samples.
    Nothing to be done (all requested files are present and up to date).
    Complete log: .snakemake/log/2022-09-22T111259.106356.snakemake.log

This is exactly the functionality that makes snakemake so useful: only do what is necessary. However, while working you might want to make sure that you are using the latest script/config/env/rules in your pipleine as, generally, changes in rules are not tracked!)

You can either force the last rule:

.. code-block:: console

    snakemake results/merged.h5ad --use-conda --force -n

.. code-block:: console

    ...
    Would remove temporary output results/QC/filtered_sample1.h5ad
    Would remove temporary output results/QC/filtered_sample2.h5ad
    Would remove temporary output results/QC/filtered_sample3.h5ad
    Job stats:
    job          count    min threads    max threads
    ---------  -------  -------------  -------------
    QC_sample        3              1              1
    merge            1              1              1
    total            4              1              1

    This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

Or the whole pipeline:

.. code-block:: console

    snakemake results/merged.h5ad --use-conda --forceall -n

.. code-block:: console

    ...
    Would remove temporary output results/QC/filtered_sample1.h5ad
    Would remove temporary output results/QC/filtered_sample2.h5ad
    Would remove temporary output results/QC/filtered_sample3.h5ad
    Job stats:
    job                 count    min threads    max threads
    ----------------  -------  -------------  -------------
    QC_sample               3              1              1
    dwn_download            1              1              1
    dwn_make_samples        1              1              1
    merge                   1              1              1
    total                   6              1              1

    This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

In both these cases, it would be advantageous to use more cores, since the QC can be run on the three samples in parallel.

With this we are done with the main functionalities of snakemake. In the next tutorial, we will change your two QC rules into a QC module.