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

This program time-aligns and splices frequency sub-band data together. It works for SIGPROC filterbank files, where each sub-band or lane contains a fraction of the total bandwidth. The resulting data are written to the merged.fil file. The time alignment is performed with better than sampling-time precision using Fourier shifting.
```

```console
$ stools-fix-chunks -h
usage: stools-fix-chunks [-h] [--mode {default,parkes}] files [files ...]

Fix multi-chunk contiguous PSRFITS search mode data.

positional arguments:
  files                 Names of search mode data files to process.

optional arguments:
  -h, --help            show this help message and exit
  --mode {default,parkes}
                        Operating mode. In 'default' mode, it recomputes the starting time by advancing it by the given NSUBOFFS offsets and resets NSUBOFFS to zero. In 'parkes' mode, it only resets the NSUBOFFS values to zero. (default:
                        default)

This program fixes time-tagging issues in multi-chunk, contiguous PSRFITS search mode data. Contiguous search mode data can be split into several files (data chunks), where the second and further files are identified by non-zero values of the NSUBOFFS keyword in the SUBINT HDU, which specifies the subint offset from the start of the first file. Later files carry the starting time of the first file. Processing software is expected to advance the starting time of the later files by the NSUBOFFS amount. However, in practice, this does not work in DSPSR (May 2026).

This program fixes the issue by advancing the starting times of the later files by their NSUBOFFS amounts and resetting their NSUBOFFS values to zero. Each file in the fixed multi-chunk dataset has NSUBOFFS = 0, but now has fully consecutive starting times as given by STT_IMJD, STT_SMJD, and STT_OFFS.

Select the operating mode depending on your data.
* Use the 'default' mode if the telescope backend advances the NSUBOFFS value, but keeps the high-resolution time stamp (STT_IMJD, STT_SMJD, and STT_OFFS) from the first file. This is the case for NRT data.

* Use the 'parkes' mode if the telescope backend advances BOTH the high-resolution multi-chunk time stamp AND the NSUBOFFS value. This is the case for Parkes UWL / Medusa data.
```

```console
$ stools-truncate-rows -h
usage: stools-truncate-rows [-h] --nrows value files [files ...]

Truncate PSRFITS search mode data.

positional arguments:
  files          Names of search mode data files to process.

optional arguments:
  -h, --help     show this help message and exit
  --nrows value  The number of the data rows to retain in the output files. (default: None)

This program truncates PSRFITS search mode data to a given number of sub:nrows of samples, specified on the command line.

This is useful for time-aligned data that differ in row count due to subtle offsets in backend process or thread synchronisation when a stop command is received. The program truncates all input files to the same number of rows so that they can be processed. A common task is to truncate all frequency lanes to the same nrows to splice them together.
```
