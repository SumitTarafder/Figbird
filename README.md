# Figbird
Figbird(FIlling Gaps Based on Insert Range Distribution) is a software designed for filling gaps in draft genome assemblies using
second generation Illumina sequencing reads.

#Dependencies:

The software can run on Linux and Mac systems with a few dependencies listed below:

1. Bowtie2: Used for mapping read pairs to gapped scaffolds. The default bowtie2 version used and given with the software is 2.2.3(Linux). But user can download their preferred version from:
 https://github.com/BenLangmead/bowtie2
- If you want to use the given version inside the software, then  unzip the folder and compile it with command “make”.

2. The software is developed in C++ and requires GNU g++(version 4.8 or greater) to compile the codes and the driver script is written in bash and requires GNU bash (version 4.3 or greater)

3. The software also needs GNU uitlity 'bc'[basic calculator]. If you don't have bc in your system, run the following command:
	 - sudo apt install bc

4. A command line JSON processor library 'jq'. You can install jq from the following github page:
	 https://stedolan.github.io/jq/

5. [Optional]Python is required only if you want to assess the quality of filled gaps using QUAST software. The exact version of QUAST along with necessary correction files as depicted in paper is already attached with the software. Unzip the folder before using it. There is no need to install QUAST.
