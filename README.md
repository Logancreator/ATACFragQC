# ATACFragQC

The Python toolkits designed to control the fragment quality of Bulk/SingCell ATAC-seq.

## Installation
~~~
python3 -m pip install --upgrade ATACFragQC
~~~

## Usage
~~~
# Basic usage
ATACFragQC [options] -i <input.bam> -r <reference.gtf>

# For more information
ATACFragQC -h
~~~

## Features
* The distrubution of fragments in chromosomes
* The distrubution of fragment lengths
* The distrubution of fragments around transcription start sites (TSSs)
* Other feature would be supported in the future ...

## Overview
![Overview of ATACFragQC](https://raw.githubusercontent.com/0CBH0/ATACFragQC/main/Images/MCBULK_qc.png)
