Modules
=======

:ref:`In the previous tutorial <qc tutorial>`, you have made two rules that take care of single-cell quality control and clustering. We will now turn it into a module.

For longer workflows, having all the rules in the ``Snakefile`` would make it completely unreadable and unmanageable. Furthermore, modules are also useful in order to compartimentalise tasks that could be reused in a different project.

.smk file
---------
Make a new file inside ``workflow/rules`` and name it ``QC.smk``. Copy the rules we have implemented previously there, and remove them from the ``Snakefile``. The latter should now look like this:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("7.0")

    configfile: "config/config.yaml"

    # define which output files you would like to generate
    rule all:
        input:
            'data/sample1.h5ad'


    module download_data:
        snakefile: "rules/download_data.smk"
        config: config

    use rule * from download_data as dwn_*

and the ``workflow/rules/QC.smk``:

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

Now that you have made a new module, you need to give its location to the Snakefile. Try to do this by taking a look at the ``download_data`` module.

.. note:: 
    The rules from module ``download_data`` are loaded and renamed with a prefix. This is not strictly neceassary, but makes it more easy to locate the rules across modules

The implementation should look something like the following. Notice that the path is relative to the ``Snakefile`` location.

.. code-block:: python

    module QC:
        snakefile: "rules/QC.smk"
        config: config

    use rule * from QC as QC_*

Now try executing the creation of the merged object again:

.. code-block:: console

    snakemake results/merged.h5ad --use-conda --force -n

This should now throw an error since it cannot locate the environment definition:

.. code-block:: console

    Building DAG of jobs...
    Updating job dwn_make_samples.
    WorkflowError:
    Failed to open source file /Users/user/Documents/Projects/snk-tutorial/workflow/rules/envs/scanpy.yaml
    FileNotFoundError: [Errno 2] No such file or directory: '/Users/user/Documents/Projects/snk-tutorial/workflow/rules/envs/scanpy.yaml'

Indeed, while the input and output paths are relative to the working directory, the ``conda`` and ``script`` directives take paths relative to the file location where the rule is defined. These should therefore be changed to:

.. code-block:: python

    rule QC_sample:
        input:
            'data/{sample}.h5ad'
        output:
            temp('results/QC/filtered_{sample}.h5ad')
        params:
            config['QC']
        conda:
            '../envs/scanpy.yaml'
        script:
            '../scripts/QC_samples.py'

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

After these changes are made, rerunning the previous command should not throw an error anymore. While this might seem a bit counterintuitive, this is neceassary when modules are located in completely different working directories, together with their scripts and environment definitions.

Modularisation
--------------

In general, it is recommended to structure a project as such that scripts, results and plots used/created by a module are put in a subfolder, e.g. ``scripts/QC``, ``results/QC`` and ``plots/QC``.

This would therefore look something like:

.. code-block:: python
    
    rule QC_sample:
        input:
            'data/{sample}.h5ad'
        output:
            temp('results/QC/filtered_{sample}.h5ad')
        params:
            config['QC']
        conda:
            '../envs/scanpy.yaml'
        script:
            '../scripts/QC/QC_samples.py'

    rule merge:
        input:
            expand('results/QC/filtered_{sample}.h5ad', sample = ['sample1', 'sample2', 'sample3'])
        output:
            data = 'results/QC/merged.h5ad',
            plot = 'plots/QC/umap.pdf'
        conda:
            '../envs/scanpy.yaml'
        script:
            '../scripts/QC/merge.py'

and a project structure:

.. code-block:: console

    ├── .gitignore
    ├── LICENSE
    ├── README.md
    ├── workflow
    │   ├── rules
    |   │   ├── download_data.smk
    |   │   ├── QC.smk
    │   ├── envs
    |   │   ├── scanpy.yaml
    │   ├── scripts
    |   │   ├── fake_samples.py
    │   │   └── QC
    |   │       ├── merge.py
    |   │       └── QC_samples.py
    |   └── Snakefile
    └── config
        └── config.yaml

In the next section, you will learn about global wildcards, and so-called checkpoints.