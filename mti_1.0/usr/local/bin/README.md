# Taxonomic Inference
> Metagenomic Taxonomy Inference Project
> | University of Central Florida


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


## GitHub Rules
1. Be clear in your commits.
  * Title should be clear and commit must have description
  * If the commit description goes over the "view limit" (it can't be contained in a single line) then make a short single sentence and explain the rest after two 'returns'
2. When deciding to edit the master branch of any file, you must fork. Avoid editing the master branch directly.
3. If you've forked from the master, you must finish that fork and merge back to master before making a new fork or editing the master.


## Programming Rules
1. Commenting Format:
  * All comments must be as clear as possible. No shorthand, etc. Pretend Dr. Yooseph is presenting your code to his colleagues (just might happen).
  * All major pieces of code (methods, classes, etc.) must be preceded with a block comment that includes a summary, used inputs, parameters, and outputs.
  * Each logical block (loops, etc.) must be preceded with inline comments
2. Naming Conventions:
  * Constants == ALL CAPS. "const int CONSTANT = 100"
  * Camel Casing for variables. "varOne"
  * Classes begin with a capital first letter. "ClassOne"
3. Indenting:
  * C++ - whatever is default
  * Python - four spaces
4. After you finish a major block of code, before merging, team performs a code review. Then you can merge.


## Resources
> References to concepts and technical details
