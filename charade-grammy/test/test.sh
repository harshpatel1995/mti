#!/usr/bin/env bash

###the script is an test example of grammy with step to step explanation, they are also echoed###
echo "###get the sample genome references file"""
if [ ! -d grefs ]
then
	echo "wget http://meta.usc.edu/softs/grammy/grefs.tgz"
	wget "http://meta.usc.edu/softs/grammy/grefs.tgz" -O grefs.tgz
	echo "tar -zxvf grefs.tgz"
	tar -zxvf grefs.tgz
fi

echo "###the script is an test example of grammy with step to step explanation###"
echo "###put grefs under this test directory###"

echo "###[first] generate the genome data needed by grammy"
echo "#let us give the output name 'test' and supply the taxids we want to include in the 'test' set"
echo "grammy_gdt test 228908,290400,292414,314225,335992,190192,272557,272844,70601,326442,312309,216895,243274,243090,167879,316279,177439,267377,178306,283942,223926,243232,196600,188937,224325,74547,84588,298386,167539,246200,186497,235909,59919,221109"

grammy_gdt test 228908,290400,292414,314225,335992,190192,272557,272844,70601,326442,312309,216895,243274,243090,167879,316279,177439,267377,178306,283942,223926,243232,196600,188937,224325,74547,84588,298386,167539,246200,186497,235909,59919,221109

echo "#now we have test.gdt and test.fna.1 and test.fna.2 three new files"
echo "#test.gdt is genome data file needed by grammy"
echo "#test.fna.1 and test.fna.2 is genome data"

echo "###[second] we generate the read data needed by grammy"
echo "grammy_rdt . ."

grammy_rdt . .

echo "#now we have two additional files"
echo "#c0_L100_N1000_S1.rdt is the read data file need by grammy"
echo "#c0_L100_N1000_S1.fasta.gz is the zipped reads fasta file"

echo "###[third] we prepare the probablity matrix file will be used by grammy"
echo "#c0_L100_N1000_S1.tblat.1 and c0_L100_N1000_S1.tblat.2 are two mapping result files"
echo "grammy_pre c0_L100_N1000_S1 test"

grammy_pre c0_L100_N1000_S1 test

echo "#now we parsed the two result files and get the c0_L100_N1000_S1.mtx file need by grammy"

echo "###[fourth] we carry out the EM esimation"
echo "grammy_em c0_L100_N1000_S1.mtx"

grammy_em c0_L100_N1000_S1.mtx

echo "#we have three new files"
echo "#c0_L100_N1000_S1.est, where the estimation is, not normalized"
echo "#c0_L100_N1000_S1.btp, where the extra bootstrap estimations are"
echo "#c0_L100_N1000_S1.lld, final log likelihood"

echo "###[fifth] we normalize the result and get the genome length and relative abundance estimates"
echo "grammy_post c0_L100_N1000_S1.est test c0_L100_N1000_S1.btp"

grammy_post c0_L100_N1000_S1.est test c0_L100_N1000_S1.btp

echo "#c0_L100_N1000_S1.avl, average genome length estimates"
echo "#c0_L100_N1000_S1.gra, genome relative abundance, first line is taxon id, second line is relative abundance, last line is error bound"
