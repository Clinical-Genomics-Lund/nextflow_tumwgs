## Nextflow tumor WGS pipeline 

The main script for running the pipeline is "main.nf". You have to change the config file accordingly to the server requirements and panel configuration. Currently, our worker server is hopper and single command line for running the pipeline is given below  

**Running the "hema" pipeline**


```console
~$ ml singularity Java nextflow/21.10.6
~$ nextflow run main.nf -profile hopper \
    --panel hema \
    -csv sample.csv \ 
    --queue grace-normal \
    -with-singularity ~/path/to/container/tumwgs_2022-02-08.sif \
    -with-report /path/to/workdir/nextflow/reports/sample.tumwgs-all.report.html \
    -with-trace /path/to/workdir/nextflow/reports/sample.tumwgs-all.trace.txt \
    -with-timeline /path/to/workdir/nextflow//reports/sample.tumwgs-all.timeline.html \
    -work-dir /path/to/workdir/nextflow/nextflow//val-6883-19.tumwgs-all
```

**Running the "solid" pipeline**

```console
~$ ml singularity Java nextflow/21.10.6
~$ nextflow run main.nf -profile hopper \
    --panel solid \
    -csv sample.csv \ 
    --queue grace-normal \
    -with-singularity ~/path/to/container/tumwgs_2022-02-08.sif \
    -with-report /path/to/workdir/nextflow/reports/sample.tumwgs-all.report.html \
    -with-trace /path/to/workdir/nextflow/reports/sample.tumwgs-all.trace.txt \
    -with-timeline /path/to/workdir/nextflow//reports/sample.tumwgs-all.timeline.html \
    -work-dir /path/to/workdir/nextflow/nextflow//val-6883-19.tumwgs-all
```

But ofcourse if you have already defined these pipelines in `pipeline_files.config` defined as follows:

```
[tumwgs-hema]
pipeline = /fs1/saile/prj/nextflow_tumwgs_dev_v2.0/main.nf -profile hopper --panel hema
container = /fs1/saile/prj/nextflow_tumwgs_dev_v2.0/container/tumwgs_2022-02-08.sif
singularity_version = 3.8.0
nextflow_version = 21.04.2
executor = slurm
cluster = grace
queue = normal

[tumwgs-solid]
pipeline = /fs1/saile/prj/nextflow_tumwgs_dev_v2.0/main.nf -profile hopper --panel solid
container = /fs1/saile/prj/nextflow_tumwgs_dev_v2.0/container/tumwgs_2022-02-08.sif
singularity_version = 3.8.0
nextflow_version = 21.04.2
executor = slurm
cluster = grace
queue = normal
```

Then running the pipeline from our own `start_nextflow_analysis.pl`starter with the bash script is as follows:

```console
~$ /path/to/your/start_nextflow_analysis.pl sample.csv

```
 In our production server the bash script for the individual file and the csv files are automatically generated thereby nextflow script is activated upon availability of the csv and shell scripts.
