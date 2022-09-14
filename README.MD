## Center for molecular diagnostics (CMD) tumor whole genome sequencing pipeline 

The main script for running the pipeline is "main.nf". You have to change the config file accordingly to the server requirements and panel accordingly to research question - i.e either hematolgy and solid tumor. Currently, our computing server is "hopper and grace" hence a profile $hopper is created with all the server descriptions (see nextflow.config for more information).Running the pipeline based on biologial question types in our CMD server (hopper/grace) is given below:  

**Running the "Hematolgy (ALL/AML) Tumor" pipeline**


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

**Running the "Solid Tumor" pipeline**

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

But ofcourse if you have already predefined these configurations in `pipeline_files.config` defined as follows the `start_nextflow_analysis.pl` ´ starter will triggers the correct pipelines

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

Under the hood `start_nextflow_analysis.pl`starter and necessary `*.csv' (usually gerated from earlier Bjorn systen)" will run in the bash as follows:

```console 
## Hematolgy
~$ /path/to/your/start_nextflow_analysis.pl sample_hema.csv

## Solid 
~$ /path/to/your/start_nextflow_analysis.pl sample_solid.csv

```
In our production server the bash script for the individual file and the csv files are automatically generated thereby nextflow script is activated upon availability of the csv and shell scripts.
