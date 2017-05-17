# Taxonomic Inference
> Metagenomic Taxonomy Inference Project
> | University of Central Florida


![Alt text](/Visualization/mock_images/heatmap.png?raw=true "Heatmap")
![Alt text](/Visualization/mock_images/violin.png?raw=true "Violin Plot")
![Alt text](/Visualization/mock_images/bar.png?raw=true "Bar Graph")


## Usage

### Visualization

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

### Troubleshooting

If you encounter issues with docker, try running the commands below and rebuilding the image.

> **WARNING**: This will stop all running containers and remove all container images

    $ docker stop $(docker ps -a -q)
    $ docker rm $(docker ps -a -q)

