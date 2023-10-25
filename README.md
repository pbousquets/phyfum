# Phyfum

[![PyPI - Version](https://img.shields.io/pypi/v/phyfum.svg)](https://pypi.org/project/phyfum)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/phyfum.svg)](https://pypi.org/project/phyfum)

-----

**Table of Contents**

- [Quick start](#quick-start)
- [Installation](#installation)
- [License](#license)

## Quick start

Phyfum allows two different workflows. If you're working with raw data (IDAT files), you can try to run phyfum in __complete__ mode. It will preprocess the files with [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html). If needed and if normal samples are available, it will also run a copy number pipeline based on [rascal](https://github.com/crukci-bioinformatics/rascal). This will allow to blacklist fCpGs that won't fluctuate as the model expects. 

In case you already have the beta values, you can run phyfum in __trees__ mode. The pipeline will simply deploy the XMLcreator tool to format the input data as expected by [BEAST](https://beast.community/) and run the inference.


## Installation

```console
pip install phyfum
```

## License

`phyfum` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

