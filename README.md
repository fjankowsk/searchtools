# SearchTools: Tools for manipulating search mode data #

[![GitHub issues](https://img.shields.io/badge/issue_tracking-GitHub-blue.svg)](https://github.com/fjankowsk/searchtools/issues/)
[![License - MIT](https://img.shields.io/pypi/l/fitpdf.svg)](https://github.com/fjankowsk/searchtools/blob/master/LICENSE)

This repository contains tools to manipulate time domain search mode data. This is relevant for pulsar or fast radio burst (FRB) data analysis.

## Author ##

The software is primarily developed and maintained by Fabian Jankowski. For more information, feel free to contact me via: fabian.jankowski at cnrs-orleans.fr.

## Installation ##

The easiest and recommended way to install the software is via the Python command `pip` directly from the `searchtools` GitHub software repository. For instance, to install the master branch of the code, use the following command:  
`pip install git+https://github.com/fjankowsk/searchtools.git@master`

This will automatically install all dependencies. Depending on your Python installation, you might want to replace `pip` with `pip3` in the above command.

## Usage ##

```console
$ stools-align-lanes -h
usage: stools-align-lanes [-h] files [files ...]

Time align and frequency splice search mode data.

positional arguments:
  files       Names of search mode data files to process.

optional arguments:
  -h, --help  show this help message and exit
```

```console
$ stools-fix-chunks -h
usage: stools-fix-chunks [-h] files [files ...]

Fix multi-chunk PSRFITS search mode data.

positional arguments:
  files       Names of search mode data files to process.

optional arguments:
  -h, --help  show this help message and exit
```

```console
$ stools-truncate-fits -h
usage: stools-truncate-fits [-h] --nrows value files [files ...]

Truncate PSRFITS search mode data.

positional arguments:
  files          Names of search mode data files to process.

optional arguments:
  -h, --help     show this help message and exit
  --nrows value  The number of the data rows to retain in the output files. (default: None)
```
