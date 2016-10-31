FROM ubuntu:16.04


# Install apt-get packages
RUN apt-get update && apt-get install -y \
	# Python & ETE dependencies
	python \
	python-pip \
	python-qt4 \
	python-lxml \
	python-six \
	python-numpy \
	# Misc.
	git-all \
	xvfb \
	nano

# Install pip packages
RUN pip install --upgrade \
	pip \
	biopython \
	ete3


# Create a new project folder
RUN mkdir -p /src && cd /src
WORKDIR /src/vis-dev