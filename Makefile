# Metagenomic Taxonomic Inference
# Authors: Trevor Ballard, Harsh Patel, Felix Sosa, Austin Vo
# University of Central Florida, Fall 2016 - Spring 2017

all:
	# Install ubuntu packages
	apt-get update && apt-get install -y \
		# BWA dependencies
		bwa \
		zlib1g-dev \
		# GRAMMy dependencies
		libstdc++6 \
		# Python & ETE dependencies
		python3 \
		python3-pip \
		python3-pyqt4 \
		python3-scipy \
		python3-biopython \
		python3-tk \
		xvfb \
		# Misc.
		git-all \
		curl \
		nano

	# Install pip packages
	pip3 install --upgrade \
		pip \
		ete3 \
		docopt \
		six \
		seaborn \
		pandas

	# Download and compile wgsim
	git clone https://github.com/lh3/wgsim.git
	cd wgsim && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm

grammy:
	# Compile GRAMMy.
	g++ -std=c++11 -o GRAMMy/grammy GRAMMy/grammy.cpp
