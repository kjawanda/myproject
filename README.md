
CSC8311 Advanced Programming for Biologists - Coursework Version 1.4. 11/03/2018

Getting Started

The program contained in the CSC8311 folder focuses on aligning 2 or more sequences, that are in fasta format.
There are two types of alignment, pairwise and multiple sequence alignment. Pairwise alignment is when two sequences are being aligned using either global or local alignment. For global alignment, a specific algorithm is used called the Needleman_Wunsch algorithm. For local alignment the algorithm used is Smith-Waterman algorithm. If there are more than two sequences that need to be aligned, tools that can be be used include Clustalw and MUSCLE. These are the tools that are incorporated within this program.

General Usage Notes

This program was produced on Ubuntu, Linux Version 17.10.1 using Oracle's VirtualBox. Programs that were initially installed onto the virtual machine include Python3 and Biopython package using the MATE Terminal:
sudo apt-get update
sudo apt-get install python3.6
sudo pip install biopython (install pip might be necessary)
pip install biopython --upgrade

Other operating systems can be used. Please refer to the companies' website to find
out how to do so. These websites can be found in the references document in the CSC8311 directory.
Softwares that need to be installed are MUSCLE, Clustalw:
sudo apt-get install clustalw
sudo apt-get install muscle

For pairwise alignment, packages are installed using python.
These are found at the top of the page of the code.
The program is found within the project called myproject, under the file name alignment.py

Running the Tests

In order to use any of the alignment functions that have been created, the files must be in fasta format.
Therefore, the first step is to make the format of the file is in fasta. If not, then the function 'format_conversion' can be used.
It takes in 2 arguments, the file name and the file format.
If the file(s)/string have 2 sequences, then pairwise alignment needs to take place. If the file/string has more than 2 sequences, multiple alignment needs to be used.

For pairwise alignment, choose whether you want global or local alignment.
For global alignment use pairwise2_global, which has two arguments, sequence1 and sequence 2 where the sequences are strings or use pairwise2_global_file where the arguments are the two$
For local alignment use pairwise2_local, which has two arguments, sequence1 and sequence2(strings) or use pairwise2_local_file where the arguments are in files.

For multiple alignment, the user has a choice between using clustalw or muscle. With clustalw, a phylogenetic tree is also created.
The function for clustalw is clustalw. For muscle, there are two functions, one for when the input is small (muscle_smallinput) and when the input is large (muscle_largeinput) .
The input is considered small if it is less than 10 and large if it is equal to or greater than 10.
The approach used for muscle_small input is simple but if you are dealomg with very large output text, all of stdout and stderr is loaded into memory as a string.
This is a potential drawback. Therefore, by using the 'subprocess' mmodule, it is possible t directly work with handles instead.

For each function included in the program, an example has been given (it has been commented out). For pairwise alignment, we use files alpha.faa, beta.faa, and haemoglobin.faa.
For MSA, opuntia.fasta was used.
To see if the format_conversion function works, you can use 'cor6_6.gb'.

Alignment program can be reached at:

Name: Kirendeep Kaur Jawanda

Website:n/a

Email: Kjawanda@hotmail.co.uk

Copyright 1994-2018 Example Corporation. All rights reserved.
This program and its uses are subject to a license agreement and are also
subject to copyright, trademark, patent and/or other laws.
All other brand and product names are trademarks or registered.
