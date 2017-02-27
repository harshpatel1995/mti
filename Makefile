# Metagenomic Taxonomic Inference
# Authors: Trevor Ballard, Harsh Patel, Felix Sosa, Austin Vo
# University of Central Florida, Fall 2016 - Spring 2017

SYSTEM ?= generic
ENVIRONMENT ?= generic

all:

environment:
	# Install ubuntu packages
	sudo apt-get update && apt-get install -y \
		bwa \
		zlib1g-dev \
		libstdc++6 \
		python3 \
		python3-pip \
		python3-pyqt4 \
		python3-scipy \
		python3-biopython \
		python3-tk \
		xvfb \
		git-all \
		curl \
		nano

	# Install pip packages
	sudo pip3 install --upgrade \
		pip \
		ete3 \
		docopt \
		six \
		seaborn \
		pandas
	
	# Clean up.
	sudo apt-get autoremove

wgsim:
	# Download and compile wgsim
	# git clone https://github.com/lh3/wgsim.git
	cd wgsim-master && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm

grammy:
	# Compile GRAMMy.
	g++ -std=c++11 -o GRAMMy/grammy GRAMMy/grammy.cpp

mainfile:
	cd MainFile && chmod 777 main.py && chmod +x main.py
