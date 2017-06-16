#!/bin/bash

declare -a genomes_of_interest=("archaea" "bacteria" "viral")

function download {
	for i in "${genomes_of_interest[@]}"
	do
		rm -r ../../genomes/"$i"
		mkdir ../../genomes/"$i"
		cd ../../genomes/"$i"
		curl "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/${i}/assembly_summary.txt" | \
		awk '{FS="\t"} !/^#/ {print $20} ' | \
		sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)|\1\2/\2_genomic.fna.gz|' > genomic_file
		mkdir fasta
		cd fasta/
		wget --input ../genomic_file
		find . -type f -exec gunzip {} +  	
		cd ../../../mti/Scripts/
	done
}

function index {
	cd ../../genomes/
	rm reference.fna
	find . -name '*.fna' -exec cat {} + > reference.fna | tee reference.fna
}

function bwa {
	cd ../../genomes/
	bwa index reference.fna
}

download
index
bwa
