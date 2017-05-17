# Metagenomic Taxonomic Inference
 Metagenomic Taxonomic Inference is an open-source end-to-end pipeline that infers and visualizes relative abundances of microbial genomes present in a metagenomic sample. The pipeline consists of three stages: read alignment, statistical inference of relative abundance and visualization. 

## Read Alignment
 [Burrows-Wheeler Aligner](https://github.com/lh3/bwa), an open-source read alignment tool, performs read mapping of the user provided single-end and paired-end reads with the provided reference genome. 
 General information about the software package and the BWT algorithm is available at [Fast and accurate short read alignment with Burrows–Wheeler transform](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705234/)

## Statistical Inference of Relative Abundance
 GRAMMy (Genome Relative Abundance using Mixture Model Theory) uses the alignment results to compute the relative abundance. GRAMMy is a C++ implementation of the Expectation Maximization (EM) algorithm described in [Accurate Genome Relative Abundance Estimation Based on Shotgun Metagenomic Reads](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027992)

## Visualization Module
 The visualization module is implemented in python using Pandas dataframes, ETE toolkit and Seaborn. Sample visualizations are shown below.

### Heatmap
![Alt text](/Visualization/mock_images/heatmap.png?raw=true "Heatmap")

### Violin Plot
![Alt text](/Visualization/mock_images/violin.png?raw=true "Violin Plot")

### Bar Graph
![Alt text](/Visualization/mock_images/bar.png?raw=true "Bar Graph")


# Usage

## Visualization

> Python 2.7

    $ python visualize.py <<GRA file>>.gra


## Docker

Build an image using `Dockerfile`:

    $ docker build -t mti-dev relative/path/to/mti/Visualization

 Then (note the `//c/` rather than `c:/`):

    $ docker run -it --rm -v //c/path/to/mti:/src/mti-dev mti-dev bash

The data in this folder will be synced to the folder in the container on the right-hand side of the colon.

If running `visualize.py` in a headless environment, prepend `xvfb-run` to the command to avoid problems with the X server. For example:

    $ xvfb-run python visualize.py

## Troubleshooting

If you encounter issues with docker, try running the commands below and rebuilding the image.

> **WARNING**: This will stop all running containers and remove all container images

    $ docker stop $(docker ps -a -q)
    $ docker rm $(docker ps -a -q)

