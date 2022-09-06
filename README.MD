## Nextflow tumor WGS pipeline 

The main script for running the pipeline is "main.nf". You have to change the config file accordingly to the server requirements and panel configuration. Currently, our worker server is hopper and single command line for running the pipeline is given below  

**Running the "ALL" pipeline**

```
nextflow run main.nf -profile hopper --panel all -csv sample.csv --queue grace-normal -with-singularity ~/path/to/container/tumwgs_2022-02-08.sif -with-report /path/to/workdir/nextflow/reports/sample.tumwgs-all.report.html -with-trace /path/to/workdir/nextflow/reports/sample.tumwgs-all.trace.txt -with-timeline /path/to/workdir/nextflow//reports/sample.tumwgs-all.timeline.html -work-dir /path/to/workdir/nextflow/nextflow//val-6883-19.tumwgs-all
```


**Running the "solid" pipeline*

```
nextflow run main.nf -profile hopper --panel solid -csv sample.csv --queue grace-normal -with-singularity ~/path/to/container/tumwgs_2022-02-08.sif -with-report /path/to/workdir/nextflow/reports/sample.tumwgs-all.report.html -with-trace /path/to/workdir/nextflow/reports/sample.tumwgs-all.trace.txt -with-timeline /path/to/workdir/nextflow//reports/sample.tumwgs-all.timeline.html -work-dir /path/to/workdir/nextflow/nextflow//val-6883-19.tumwgs-all
```