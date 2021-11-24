
"""
	OVERVIEW

	This script extends and adjusts the metrics reported by QUAST by parsing the contig_reports*.stdout
	file of alignments. This module does:
	1 In addition to QUASTs report of total mismatches and indels, this script calculates the total length
		of local misassemblies (summing the lengths of the "local misassembly" in contig_reports*.stdout).
	2 Corrects extensive misassemblies caused by small (<N bp) sequences aligning eleswhere and report them
	as local. They are added to the total length of local assemblies and the number of extensive misassemblies
		reported by QUAST is adjusted. More information and examples about these cases below. 
	3 Sums up the total erronious sequence 'e' calculated from all the lengths in 1 and 2 (i.e. total length of
		local assemblies).
	4 Output: Total sum of erronoius sequence and total number of extensive misassemblies after "correction".

	ABOUT 3:

	Calculating total sum:

	Total 'local' erronius sequence length as reported by QUAST is either mismatches, indels or local 
	missassemblies. Total erronius sequence length 'e' is calculated as:

	e = nr_mismatches + sum_i [small_indel_i * len(small_indel_i) ] 
		+ sum_i [large_indel_i * len(large_indel_i) ] 
		+ sum_{i in local misassembly} [local_misassembly_i * len(local_misassembly_i)]
 

	ABOUT 2:

	QUAST reporta a break point between two adjacent alignments as a local misassembly if the 
	alignmnet coordinates between left and right flanking sequence of adjacent alignments
	differs with < Nbp (N=1000 in QUAST). If they differ with more than N bp, it is classified as a
	misassembly. see section 1.1  in 
	http://bioinformatics.oxfordjournals.org/content/suppl/2013/02/18/btt086.DC1/quast_supplement_updated.pdf.

	We report a break point between two alignments (not neccessarily adjacent) as a local misassembly if the 
	alignmnet coordinates between left and right flanking sequence of two alignments, a_i and a_(i+k)
	differs with < Nbp (N=1000 in QUAST) and alignments lenght[a_(i+1)+..+ a_(i+k-1)] < N . 
	
	Example: Assume scaffold ABC. We classify B to be a local assembly error if 
	1 B is misplaced and the length of the piece B is smaller than N bp.
	2 Flanking contigs A and C are placed correctly with respect to each other (ie the distance on the genome
	compared to the distance in the scaffold is less than Nbp)  

	This opens up for adjusting small pieces aligned at distant locations on the genome to still 
	be classified as local errors. The difference becomes clear when studying Example 1 and 2 below 
	(both with real data in real scenario). Note also that B can be made up of several small 
	alignemnents (as in Example 2 below).

	METHOD:

	More specifically, we parse all regions where "Extensive misassemblies occur" and call this a 
	missassembly region. Note that a "misassembly region" can be made up by several 
	alignments a_err = a_i,..a_j. For each such misassembly region, we check the 
	distance discrepancy "d" between the flanking left and right sequence to a_err, namely
	alignment a_(i-1) and a_(r+1). We correct the extensive missassembly(ies) if 
	max( len(a_err), d) < N. That is, the total length of the small fragmemts and the distance discrepancy
	between the flanking sequence of that region needs to be less than N bp. Otherwise we let the region
	be, i.e. quasts default reporting is left alone.

	Main algorithm of this "correction" is found in the function "get_sum_large_misassemblies(infile, N)".
	There are three corner cases:
	1 First alignment in a contig/scaffold is incorrect. 
	2 Last alignment in a contig/scaffold is incorrect. 
	3 Scaffold is only made up by 2 alignments with extensive misassembly between them.


	EXAMPLES

	Theoretical example:

	Consider a scaffold ABCDEF with B,D,F < N bp
	if true alignement are ACE and BDF with BDF being "far away from" ACE then
	QUAST would give 5 extensive misassemblies.
	With this script, we would give 3 local misassembies, namely the "small" pieces
	B,D,F are all wrongly placed.

	Biological data examples:

    Example1:
    From GAGE SGA assembly on S aureus, ABCDEFG (Real Alignment 4-10) with B,E <= 1000 giving 4 misassemblies:
    We consider it 2 local misassemblies with total length 416+607.

    Real Alignment 4: 2479002 2539931 | 144790 83870 | 60930 60921 | 99.77 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-169_pilon
      Extensive misassembly ( relocation, inconsistency = 46802 ) between these two alignments
    Real Alignment 5: 2431874 2432288 | 145117 144702 | 415 416 | 99.28 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-169_pilon
      Extensive misassembly ( relocation, inconsistency = -46838 ) between these two alignments
    Real Alignment 6: 2451046 2478582 | 172727 145247 | 27537 27481 | 99.79 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-169_pilon
      Gap between these two alignments (local misassembly). Inconsistency = 33
    Real Alignment 7: 2432201 2450891 | 191539 172849 | 18691 18691 | 99.99 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-169_pilon
      Extensive misassembly ( relocation, inconsistency = 1005169 ) between these two alignments
    Real Alignment 8: 1426513 1427119 | 192057 191452 | 607 606 | 98.52 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-169_pilon
      Extensive misassembly ( relocation, inconsistency = -1005139 ) between these two alignments
    Real Alignment 9: 2405671 2431741 | 218048 191968 | 26071 26081 | 99.9 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-169_pilon
      Gap between these two alignments (local misassembly). Inconsistency = 2
    Real Alignment 10: 2400323 2405226 | 223394 218491 | 4904 4904 | 100.0 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-169_pilon

    Example2:
    ABCD with sum(B,C) < N bp:

    Real Alignment 7: 2665021 2698092 | 197133 164062 | 33072 33072 | 99.99 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-125_pilon
      Extensive misassembly ( relocation, inconsistency = 1816608 ) between these two alignments
    Real Alignment 8: 848396 848466 | 197150 197080 | 71 71 | 98.59 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-125_pilon
      Extensive misassembly ( relocation, inconsistency = -1706970 ) between these two alignments
    Real Alignment 9: 2555227 2555349 | 197289 197167 | 123 123 | 97.56 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-125_pilon
      Extensive misassembly ( relocation, inconsistency = -109632 ) between these two alignments
    Real Alignment 10: 2654420 2664908 | 207728 197240 | 10489 10489 | 100.0 | NC_010079_Staphylococcus_aureus_subsp._aureus_USA300_TCH1516__complete_genome scaffold-125_pilon

	Author: Kristoffer Sahlin
			ksahlin@kth.se 
"""

import re
import argparse
import sys
import os


	
def is_local_middle(lines,i,N):
	genome_coordinate_diff = 0
	tot_len_almt_j = 0 
	tot_corrected_extensive_errors = 0
	j = 1
	#print lines[i-1].split('|')[0].strip().split(' ')
	left_reference_coords = map(lambda x: int(x), lines[i-1].split('|')[0].strip().split(' ')[-2:])
	left_scaf_coords = map(lambda x: int(x), lines[i-1].split('|')[1].strip().split(' '))

	while not re.match('CONTIG:(.)+' , lines[i+j].strip()):
		is_ext_error = re.match('Extensive misassembly (.)+' , lines[i+(j-1)].strip())
		if is_ext_error:
			tot_corrected_extensive_errors += 1

		right_reference_coords = map(lambda x: int(x), lines[i+j].split('|')[0].strip().split(' ')[-2:])
		right_scaffold_coords = map(lambda x: int(x), lines[i+j].split('|')[1].strip().split(' '))

		genome_coordinate_diff = min( abs(left_reference_coords[0]- right_reference_coords[1]) , 
										abs(left_reference_coords[1]- right_reference_coords[0]) )

		scaf_coordinate_diff =  min( abs(left_scaf_coords[0]- right_scaffold_coords[1]) , 
										abs(left_scaf_coords[1]- right_scaffold_coords[0]) )
		#print 'HHHHH:',genome_coordinate_diff, scaf_coordinate_diff, tot_len_almt_j
		local_error_seq_length = abs(genome_coordinate_diff - scaf_coordinate_diff)
		if  local_error_seq_length <= N and genome_coordinate_diff <= N and tot_len_almt_j <= N:
			#print 'local middle', tot_len_almt_j
			incorrect_seq = max(local_error_seq_length, tot_len_almt_j)
			return tot_corrected_extensive_errors, incorrect_seq, i+j

		elif tot_len_almt_j > N:
			return False, 0, i+j

		tot_len_almt_j += max( map(lambda x: int(x), lines[i+j].split('|')[2].strip().split(' ')))

		j+=2

	# Error started in the middle but continued until end of scaffold. Then
	# the only error sequence length to consiter is the piece from where the error started up 
	# until the end of the scaffold. Of course, it needs to be <= Nbp in order to be reported.
	if tot_len_almt_j <= N :
		incorrect_seq =  tot_len_almt_j
		#print 'several local end errors', tot_len_almt_j
		return tot_corrected_extensive_errors, incorrect_seq, i+j		
	else:
		return  False, 0, i+j

def is_local_end(lines,i,N):
	j = 1
	tot_len_almt_j = max( map(lambda x: int(x), lines[i+j].split('|')[2].strip().split(' ')))
	if tot_len_almt_j <= N:
		#print 'local end:', tot_len_almt_j
		return 1, tot_len_almt_j, i+j
	else:
		return False, 0, i+j


def is_local_beginning(lines,i,N):
	j = 1
	tot_len_almt_j = max( map(lambda x: int(x), lines[i-1].split('|')[2].strip().split(' ')))
	if tot_len_almt_j <= N:
		#print 'local beginning:', tot_len_almt_j
		return 1, tot_len_almt_j, i+j
	else:
		return False, 0, i+j

def is_local_2_aligmnets(lines,i,N):
	j = 1
	tot_len_almt_1 = max( map(lambda x: int(x), lines[i-1].split('|')[2].strip().split(' ')))
	tot_len_almt_2 = max( map(lambda x: int(x), lines[i+j].split('|')[2].strip().split(' ')))
	misassm_length = min(tot_len_almt_1,tot_len_almt_2)
	if misassm_length <= N:
		#print 'local one out of two almnts:', misassm_length
		return 1, misassm_length, i+j
	else:
		return False, 0, i+j


def get_sum_large_misassemblies(infile, N):
	lines = infile.readlines()
	i = 0
	tot_sum_extensive_to_local = 0
	tot_errors_corrected = 0
	while i < len(lines):
		misassm = re.match('Extensive misassembly (.)+' , lines[i].strip())
		if misassm:
			beginning = re.match('Real Alignment 1:(.)+' , lines[i-1].strip())
			left_length = max( map(lambda x: int(x), lines[i-1].split('|')[2].strip().split(' ')))
			#right_length = max( map(lambda x: int(x), lines[i+1].split('|')[2].strip().split(' ')))
			end = re.match('CONTIG:(.)+' , lines[i+3].strip())
			end_all = re.match('Analyzing(.)+' , lines[i+3].strip())

			# If contig scaffold only consists of two alignemnts
			# check it either one of them is local, i.e. less than N bp
			if beginning and (end or end_all):
				extensive_missassembly_corrected, mis_assm_length, i = is_local_2_aligmnets(lines,i,N)

			# if long beginning, check if middle alignment is local
			elif beginning and left_length > N:
				extensive_missassembly_corrected, mis_assm_length, i = is_local_middle(lines,i,N)

			# if short beginning check if first alignemtn is local
			elif beginning:
				extensive_missassembly_corrected, mis_assm_length, i = is_local_beginning(lines,i,N)

			# if end alignment, check if end alignment is local
			elif end or end_all:
				extensive_missassembly_corrected, mis_assm_length, i = is_local_end(lines,i,N)

			# misassembly somewhere in middle of contig/scaffold
			else:
				#print 'lolz'
		 		extensive_missassembly_corrected, mis_assm_length, i = is_local_middle(lines,i,N)

		 	if extensive_missassembly_corrected:
		 		#print 'NOOOO', mis_assm_length
		 		tot_errors_corrected += extensive_missassembly_corrected
		 		tot_sum_extensive_to_local += mis_assm_length


		i+=1

	return tot_sum_extensive_to_local, tot_errors_corrected


def get_sum_local_misassemblies(infile):
	# Gap between these two alignments (local misassembly). Inconsistency = 12
	l = ''.join(infile.readlines())
	tot_length = 0
	nr_local = 0
	local_misassemblies =  re.findall('\(local misassembly\). Inconsistency = [-\d]+', l)
	for hit in local_misassemblies:
		tot_length += abs(int(hit.split()[-1]))
		nr_local += 1
	return tot_length , nr_local

def get_sum_indels(infile):
	indel_length = int(infile.readlines()[-1].strip().split()[-1])
	return indel_length

def get_mismatches(infile):
	mismatch_length = int(infile.readlines()[-5].strip().split()[-1])
	return mismatch_length

def count_extensive_misassemblies(infile):
	l = ''.join(infile.readlines())
	return len(re.findall('Extensive misassembly (.)+', l))



def parse_quast_report(infile):
	out_string = ''
	lines = map( lambda line: line.strip().split(), infile)
	filtered_lines=[]
	#print lines[9]	
	
	
	if lines[10][0] == 'Reference' :
	    index2 = 10
	    index1 = 2
	else: 
	    index2 = 18
	    index1 = 5
	
	#print lines[9][index1]
	#print lines[index2][2]
	#print int(lines[9][index1])/float(lines[index2][2])
	
	#print index1
	#print index2
	
	if int(lines[9][index1])/float(lines[index2][2]) >= 0.75:
		if index2 == 10:
		    unaligned_row = 26
		    NGA_row = 34
		else: 
		    unaligned_row = 37
		    NGA_row = 46
	else:
		if index2 == 10:
		    unaligned_row = 24
		    NGA_row = 32
		else: 
		    unaligned_row = 35
		    NGA_row = 44
		
    
	for i,line in enumerate(lines):
		content = filter(lambda x: x != '',line )
		if unaligned_row == i:
			#print line
			unaligned_len = int(content[2])
		if NGA_row == i:
			NGA50 = int(content[1])
	return unaligned_len, NGA50



import re
import argparse
import sys

def fasta_iter(fasta_file):
    """
        Reads a fasta file into memory.

        Arguments:
        fasta_file - A python file object. The file should be in 
        fasta format.

        Returns:
            an iterator over accession, sequence.

    """  

    k = 0
    temp = []
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            k += 1
        elif line[0] == '>':
            temp = ''.join(temp)
            yield accession, temp
            temp = []
            accession = line[1:].strip().split()[0]
        else:
            temp.append(line.strip())
    
    temp = ''.join(temp)
    yield accession, temp

def get_gapstats(fasta_file):
	nr_gaps = 0
	tot_gap_length = 0
	gap_distr = []
	for acc, seq in fasta_iter(fasta_file):
		gap_list = re.findall('[Nn]+', seq)
		gap_lengths = map(lambda x: len(x), gap_list)
		gap_distr += gap_lengths
		tot_gap_length += sum(gap_lengths)
		nr_gaps += len(gap_lengths)

	#print 'nr_gaps\ttot_gap_length'
	#print '{0},{1}'.format(nr_gaps,tot_gap_length)
	return nr_gaps, tot_gap_length


def main(args):
	#if not os.path.exists(args.outpath):
	#	os.mkdir(args.outpath)

	err = False
	try:
		open(args.quast_stdout,'r')
	except IOError:
		#print '{0} does not exist. skipping...'.format(args.quast_stdout)
		err=True
	try:
		open(args.quast_misassemblies,'r')
	except IOError:
		#print '{0} does not exist. skipping...'.format(args.quast_misassemblies)
		err=True
	try:
		open(args.quast_report,'r')
	except IOError:
		#print '{0} does not exist. skipping...'.format(args.quast_report)
		err=True

	try:
		open(args.scaffolds,'r')
	except IOError:
		#sys.stderr.write('{0} does not exist. skipping...\n'.format(args.scaffolds))
		#print '{0},{1}'.format('-','-')
		err=True

	if not err:
		tot_extensive_reduced, nr_extensive_corrected = get_sum_large_misassemblies(open(args.quast_stdout,'r'),args.N)
		tot_length_local, nr_local = get_sum_local_misassemblies(open(args.quast_stdout,'r'))
		tot_length_indel = get_sum_indels(open(args.quast_misassemblies,'r'))
		tot_length_mismatch = get_mismatches(open(args.quast_misassemblies,'r'))
		tot_extensive_before = count_extensive_misassemblies(open(args.quast_stdout,'r'))
		unaligned_length, NGA50 = parse_quast_report(open(args.quast_report,'r'))
		number_of_gaps, total_gap_length = get_gapstats(open(args.scaffolds, 'r'))

		misassemblies = tot_extensive_before - nr_extensive_corrected
		erroneous_length =  tot_length_mismatch + tot_length_indel + tot_length_local + tot_extensive_reduced

		quality_string = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(misassemblies, erroneous_length, unaligned_length, NGA50, number_of_gaps, total_gap_length )
		print quality_string

		# print 'Number of "extensive errors" reduced from misassemblies to local misassemblies: {0}\nTotal extensive to local error length: {1}'.format(nr_extensive_corrected,tot_extensive_reduced)
		# print 'Total misassemblies before {0}. Total misassemblies after {1}'.format(tot_extensive_before, tot_extensive_before - nr_extensive_corrected)
		# print 'Local errors length: {0}\nTotal local errors found: {1}'.format(tot_length_local, nr_local)
		# print 'Total indel length: {0}'.format(tot_length_indel)
		# print 'Total mismatch length: {0}'.format(tot_length_mismatch)
		# error_seq_length = tot_length_mismatch+tot_length_indel+tot_length_local+tot_extensive_reduced

		# if args.supplementary:
		# 	print 'misassemblies,error-seq-length,unaligned-length,NGA50,len-mismatch,len-indel,len-local'
		# 	quality_string = '{0},{1},{2},{3},{4},{5},{6}'.format(tot_extensive_before - nr_extensive_corrected, error_seq_length, unaligned_len,NGA50,tot_length_mismatch,tot_length_indel,tot_length_local )
		# else:
		# 	print 'misassemblies,error-seq-length,unaligned-length,NGA50'
		# 	quality_string = '{0},{1},{2},{3}'.format(tot_extensive_before - nr_extensive_corrected, error_seq_length, unaligned_len,NGA50)

		# print quality_string
		# print >> open(os.path.join(args.outpath,'gap_filling_quality_eval.csv'),'w' ), quality_string

	else:
		print '-\t-\t-\t-\t-\t-'




if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Sum up total erronous sequence lengths reported by quast.")
	parser.add_argument('quast_stdout', type=str, help='A quast contig_reports/contig_reports.stdout file.')
	parser.add_argument('quast_misassemblies', type=str, help='A quast contig_reports/missassemblies.txt file')
	parser.add_argument('quast_report', type=str, help='A quast report.txt file.')
	parser.add_argument('scaffolds', type=str, help='A fasta file with scaffolds.')	
	#parser.add_argument('outpath', type=str, help='Folder for output.')

	parser.add_argument('--N', dest='N',type=int, default=1000, help='How large a large misassembly should be regarded as.')
	parser.add_argument('--skip_large_misassemblies', dest='skip_large_misassemblies' ,action='store_true', help='Skip calc length of large misassemblies.')
	parser.add_argument('--skip_indels', dest='skip_indels', action='store_true', help='Skip calc length of large misassemblies')
	parser.add_argument('--suppl', dest='supplementary', action='store_true', help='Skip calc length of large misassemblies')


	args = parser.parse_args()
	main(args)
