# Figbird
Filling gaps based on insert range distribution
- A software developed in C++ for filling gaps in draft genome assemblies using second generation Illumina sequencing reads.
- Supports read pairs of both smaller inserts(~200 bp) and larger inserts (~3500 bp).
- Utilizes probabilistic methods instead of graph based methods based on insert size information of read pairs.
- Makes maximum use of available sequence information by using both partially aligned reads as well as unampped reads.

## Dependencies
The software can run on Linux and Mac systems with a few dependencies listed below:

- <strong>Bowtie2</strong>: Used for mapping read pairs to gapped scaffolds. The default bowtie2 version used and provided within the software is 2.2.3(Linux). But users can download their preferred version from https://github.com/BenLangmead/bowtie2. If you want to use the version of bowtie inside the software, then  give the following commands:
  ```
  unzip bowtie2-2.2.3-source.zip
  cd bowtie2-2.2.3/
  make
  ```
- The software is developed in C++ and requires <strong>GNU g++</strong>(version 4.8 or greater) to compile the codes and the driver script is written in bash and requires GNU bash (version 4.3 or greater)
- <strong>GNU uitlity 'bc'</strong> [basic calculator]. If you don't have bc in your system, run the following command:
```
sudo apt install bc
```
- A command line JSON processor library <strong>'jq'</strong>. You can install jq from the following github page https://stedolan.github.io/jq/
- [Optional] <strong>Python</strong> is required only if you want to assess the quality of filled gaps using <strong>QUAST</strong> software. The exact version of QUAST along with necessary correction files as depicted in paper is already attached with the software. Unzip the folder before using it. There is no need to install QUAST.

## Input configuration 
Figbird uses a configuration file in JSON format to take scaffold and read library information. A sample configuration file named "Config.json" is provided in the installation folder which must be editted accordingly. Users can also use this website (http://jsonlint.com/) to check the validity of their input config file. Following is the list of required information specified in the JSON file with explanations:

- <strong>Draft_genome</strong>: Path to the gapped draft genome to fill.
- <strong>Bowtie2</strong>: Path to the bowtie2 executables. If you use the bowtie2 version inside the folder, then put "bowtie2-2.2.3", otherwise if it is in system path, then put “” in the path. In case you want to use your preferred version, put that path of installed directory in input.
- <strong>Output_Folder</strong>: Path to the directory where all the outputs will be stored. 
- <strong>Reference_Genome[optional]</strong>: Only needed if you want to evaluate the quality of the filled assembly using QUAST.
- <strong>Read_Pairs</strong>: Input your paired read libraries one by one along with all the following information:
	
	1. <strong>path_1</strong>: Path to first of the read pair files.
 	2. <strong>path_2</strong>: Path to second of the read pair files.
	3. <strong>avg_insert_size</strong>: Average insert size of the read pair library.
	4. <strong>is_reverse</strong>: If your read pair files are already in forward-reverse(FR) orientation then put 0, otherwise put 1. In case a 1 is given, we will reverse complement both the files of the the input read pair.
	5. <strong>max_read_len</strong>: Maximum read length of the library.
	6. <strong>serial_num</strong>: The order of reads usage for filling gaps.
	7. <strong>num_itr_partial</strong>: We will use both one end partially aligned and one end unmapped reads for each read pair for gap filling purpose. Enter the itration count for partial approach here.
	8. <strong>num_itr_unmapped</strong>: Enter the itration count for unmapped approach here.
	9. <strong>order</strong>: Put the order for Which one between partial and unmapped method will be applied first.
	
	* [Users must input atleast one library of read pair files and all 9 required information per library to start gap filling] 
	
- <strong>Parameters</strong>: 
	1. <strong>numthreads</strong>: Number of threads used during bowtie2 alignment and gap filling procedure.[Default:4]
	2. <strong>evaluation</strong>: Put 1 if you want to assess with QUAST or 0 otherwise.[Default:0]
	3. <strong>gaplen_negative_overlap</strong>: We have allowed negative gap lengths in our method i.e a gap can be diminished if the corresponding left and right flank has an overlap with supporting verification of aligned reads. Enter the maximum length of the gaps for which this method will be applicable.[Default: 30]
	4. <strong>default</strong>: If you want to manually input the order of the reads usage along with their number of iterations, put 0. Otherwise, put 1 for default approach. If you put 1, then information [6-9] for read pairs won't be needed to specify.[Default:1]
	5. <strong>trim_len</strong>: The amount of nucleotides being chopped off from either side of the gapped regions as this is the stopping point for the assemblers and highly likely to contain erroneous sequence.[Deafult:10]
	
## Running Figbird
Download the folder https://github.com/SumitTarafder/Figbird. Users can directly run the tool if all the dependencies are installed beforehand.
To run Figbird
```
tar xzf Figbird.tar.gz
cd Figbird
chmod a+x RunFigbird.sh && ./RunFigbird.sh Config.json
```

## Citations
