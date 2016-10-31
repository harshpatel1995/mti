FROM ubuntu:16.04


# Install ubuntu packages
RUN apt-get update && apt-get install -y \
	# BWA dependencies
	bwa \
	zlib1g-dev \
	# GRAMMy dependencies
	libstdc++6 \
	# Python & ETE dependencies
	python \
	python-pip \
	python-qt4 \
	python-lxml \
	python-six \
	python-numpy \
	python-scipy \
	python-setuptools \
	python-dev \
	python-biopython \
	xvfb \
	# Misc.
	git-all \
	curl \
	nano 

# Install pip packages
RUN pip install --upgrade \
	pip \
	ete3


# Download and compile wgsim
RUN git clone https://github.com/lh3/wgsim.git
RUN cd wgsim && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm

# Create a new project folder
RUN mkdir -p /src && cd /src
WORKDIR /src/mti-dev