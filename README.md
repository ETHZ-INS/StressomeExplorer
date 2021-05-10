# StressomeExplorer

This app enables the inspection of raw data associated with the publication (REF) on a per gene basis.
The user can input the name of a gene of interest and data across datasets will be presented in a intuitive way.
There are three different types of datasets that can be inspected, transcriptomic, proteomic and phospho-proteomic.

To install the app on your local computer, simply downlad the ziped repository and extract its contents. Ensure R is installed on your computer. Then, install all required packages (see package versions below).

To launch the app, use:
```{r}
shiny::runApp('/path/app')
```

A new window will open with a pre-selected gene and associated plots. on the top right, a selection box will enable to change the gene. Initialization of the app may take a few seconds (<20), changing the gene may take a few seconds (<10).

a runing web based version (no installation required) can be found under https://bohaceklab.hest.ethz.ch/StressomeExplorer

Further, scripts and results related to the publication (REF) can be found in the /Results and /Scripts Folders

## System Info (used for testing and compatible with the code)

R version 4.0.3 (2020-10-10)

Platform: x86_64-pc-linux-gnu (64-bit)

Running under: Ubuntu 18.04 LTS

## Package versions

shiny_1.6.0

SummarizedExperiment_1.20.0

SEtools_1.4.9

ggplot2_3.3.3

cowplot_1.1.1

## Pseudocode

This software enables visual inspection of raw data. Data in the form of Rds files containing summarized experiment files for proteomic and transcriptomic data are loaded first. the software then initializes to one gene (Fos) and plots all data associated with this gene across multiple tabs. Further, it determines all available genes in any datasets and enables the user to select any of these. Then, the software will plot any data available for this gene on multiple tabs. Further, data is shown as log transformed or raw abundace based on the users input.
