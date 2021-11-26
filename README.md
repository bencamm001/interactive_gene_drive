# interactive_gene_drive
Local Interactive Gene Drive Model

Needs at least Julia v1.5 and Plots v1.19

# ReadMe

This script creates an interactive gene drive model at localhost:8000.

## Quick Install
```
/path/to/julia -i /path/to/julia_shiny_model18.jl
```


## Quick Guide
The first time running the script will take some time to install necessary Julia packages. 

The solid blue line is the allele frequency of the drive, the green dashed line is the wild-type allele frequency and the red dashed line is the resistance allele frequency.

![](gene_drive_model.png)
