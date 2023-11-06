# Phyfum

[![PyPI - Version](https://img.shields.io/pypi/v/phyfum.svg)](https://pypi.org/project/phyfum)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/phyfum.svg)](https://pypi.org/project/phyfum)

-----

**Table of Contents**

- [Quick start](#quick-start)
- [Installation](#installation)
- [Preparing the sample sheet / metadata](#preparing-the-sample-sheet--metadata)
- [License](#license)

## Quick start

Phyfum allows two different workflows. If you're working with raw data (IDAT files), you can try to run phyfum in __complete__ mode. It will preprocess the files with [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html). If needed and if normal samples are available, it will also run a copy number pipeline based on [rascal](https://github.com/crukci-bioinformatics/rascal). This will allow to blacklist fCpGs that won't fluctuate as the model expects. 

An example run for this workflow would look like this:

```{bash}
phyfum run --input /path/to/epic_array --output experiment1 --workdir experiment1 --patientinfo /path/to/epic_array/sample_sheet.csv --patient-col patient --age-col age --patient-col Patient --sample-col Sample --stemcells 3-10-3 
```

In case you already have the beta values, you can run phyfum in __trees__ mode. The pipeline will simply deploy the XMLcreator tool to format the input data as expected by [BEAST](https://beast.community/) and run the inference.

```{bash}
phyfum run --input /path/to/exampleBeta.csv --output onlybetas --workdir onlybetas --patientinfo /path/to/meta.csv --patient-col patient --age-col age --stemcells 3-10-3
```

Phyfum detect what kind of input are provided and select automatically what workflow is needed.

**Important**. This a preliminary version of the pipeline. Right now, only samples of same age should be provided (per individual). The pipeline won't work properly otherwise.

## Installation

Make sure you have installed python3 and R (>4.0.0). Then simply run:

```console
pip install phyfum
```

The preprocessing of IDAT files relies on Minfi, an R package. Missing dependencies will be installed during the first run, so it may take longer than usual. You can also install them yourself with:

```{r}
if (!require("pacman")) install.packages("pacman")
p_load(optparse, cli, conumee, minfi, parallel, tibble, tidyr, dplyr, data.table, gtools)
p_load_gh("crukci-bioinformatics/rascal")
```

## Preparing the sample sheet / metadata

Phyfum relies on the Array Sample sheet for the complete workflow and a custom metadata file for the short workflow. In any case, the file must be a comma-separated file. 

- __Sample sheet__. Please, make sure you have a column called 'Basename', that leads to the idat files of each sample. For instance, suppose the folder tree of your array looks like this:


```
├── 200360145522
│   ├── 200360145522_R15C01_Grn.idat --> Sample: A
│   ├── 2003601455222_R15C01_Red.idat --> Sample: A
├── 200531490071
│   ├── 200531490071_R02C01_Grn.idat
│   ├── 200531490071_R02C01_Red.idat
└── allSampleSheet.csv
```

The Basename for Sample A would be: _200360145522/200360145522_R15C01_. This is required by Minfi when it performs the matching of the IDAT files with the sample sheet metadata.

Additionally, make sure that a column called 'Group' exists, to differentiate tumor and normal samples. The pipeline will try to find how many "normal" or "control" samples exist to use them as controls for the CNV pipeline. If you don't have any normals, this part of the pipeline will be skipped.

> We'll soon parametrize the column name for the group, so it can be named differently.

- __Custom metadata__. It doesn't require anything special as long as it is in CSV format. In order to identify what the columns are, you can use the arguments `--patient-col`, `--sample-col` and `age-col`, if the column names in your file are different from the defaults.

> Both the custom metadata and the sample sheet are passed through `--patientinfo`.



## License

`phyfum` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

