CSC8311 Advanced Programming for Biologists - Coursework Version 1.4. 11/03/2018

Getting Started

The program contained in the CSC8311 folder focuses on aligning 2 or more sequences, that are in fasta format.
There are two types of alignment, pairwise and multiple sequence alignment. Pairwise alignment is when two sequences are being aligned 
using either global or local alignment. For global alignment, a specific algorithm is used called the Needleman_Wunsch algorithm. For local
alignment the algorithm used is Smith-Waterman algorithm. If there are more than two sequences that need to be aligned, tools that can be 
be used include Clustalw and MUSCLE. These are the tools that are incorporated within this program. 

General Usage Notes
This program was produced on Ubuntu, Linux Version 17.10.1 using Oracle's VirtualBox. Programs that were initially installed onto the 
virtual machine include Python3 and Biopython package using the MATE Terminal:
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

Running the Tests
Questions are asked throughout the program. Answer these in lower case ALWAYS!
The first question is about making sure the file(s) are in the correct format. 
The next questions asks if you need pairwise or multiple alignment. 
If pairwise is chosen, the question asks if you want global or local alignment.
If multiple alignment is chosen, the question asks which software would you like to use.
finally, it tells you which function to use. 

Alignment program can be reached at:
Name: Kirendeep Kaur Jawanda
Website: 
Email: Kjawanda@hotmail.co.uk

Copyright 1994-2018 Example Corporation. All rights reserved. 
This program and its uses are subject to a license agreement and are also
subject to copyright, trademark, patent and/or other laws. 
All other brand and product names are trademarks or registered. 
