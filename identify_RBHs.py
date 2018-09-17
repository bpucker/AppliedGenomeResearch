### Boas Pucker ###
### v0.2 ###
### bpucker@cebitec.uni-bielefeld.de ###

### publicly available: https://www.researchgate.net/post/how_to_do_curating_miss_annotation_of_genome_analysis ###

__usage__ = """ python identify_RBHs.py\n
				
				--prefix <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS>\n
				--input1 <FULL_PATH_TO_INPUT_FILE1>\n
				--input2 <FULL_PATH_TO_INPUT_FILE2>\n
				--seq_type <'nucl'|'prot'>
				
				feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""

import re, os, sys
from operator import itemgetter


# --- end of imports --- #

def load_results_from_BLAST_result_file( BLAST_result_file, cutoff=0.99 ):
	"""! @brief load data from BLAST result file """
	
	data = {}
	
	with open( BLAST_result_file, "r" ) as f:
		line = f.readline()
		prev_query = line.split('\t')[0]
		hits = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_query:
				sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
				if len( sorted_hits ) > 1:
					if ( sorted_hits[-2]['score'] / sorted_hits[-1]['score'] ) < cutoff:
						data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
				else:
					data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
				hits = []
				prev_query = parts[0]
			hits.append( { 'query': parts[0], 'subject': parts[1], 'score': float( parts[-1] ) } )
			line = f.readline()
		sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
		if len( sorted_hits ) > 1:
			if ( sorted_hits[-2]['score'] / sorted_hits[-1]['score'] ) > cutoff:
				data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
		else:
			data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
	print "entries in data: " + str( len( data.keys() ) )
	return data


def compare_datasets( data1, data2, outputfile ):
	"""! @brief compares datasets and identifies bidirectional best hits """
	
	seq_IDs_of_interest = []
	
	counter = 0
	keys = data1.keys()
	with open( outputfile, "w" ) as out:
		out.write( "seq_IDs_of_file1\tseq_IDs_of_file2\n" )
		for key in keys:	#key=candidate gene
			try:
				value = data1[ key ]	#value=contig_ID
				try:
					other_value = data2[ value ]	#other_value=candidate_gene_ID
					if key == other_value:
						counter += 1
						out.write( key + '\t' + value + '\n' )
						seq_IDs_of_interest.append( value )
				except:
					pass
			except:
				pass
	print "number of matches: " + str( counter )
	return seq_IDs_of_interest


def load_multiple_fasta_file( fasta_file ):
	"""!@brief load content of multiple fasta file """
	
	content = {}
	with open( fasta_file, "r" ) as f:
		header = f.readline().strip()[1:]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				content.update( { header: seq } )
				header = line.strip()[1:]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		content.update( { header: seq } )
	return content


def write_seqs_of_interest_into_new_file( seq_of_interest_file, seq_IDs_of_interest, sequences ):
	"""! @brief write all sequences of interest into new output file """
	
	with open( seq_of_interest_file, "w" ) as out:
		for ID in seq_IDs_of_interest:
			seq = sequences[ ID ]
			out.write( '>' + ID + '\n' + seq + '\n' )


def main( parameters ):
	"""! @brief identifies RBHs between given data sets """
	
	prefix = parameters[ parameters.index( '--prefix' )+1 ]
	if prefix[-1] != '/':
		prefix += '/'
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	seq_file1 = parameters[ parameters.index( '--input1' )+1 ]
	seq_file2 = parameters[ parameters.index( '--input2' )+1 ]
	
	if not os.path.isfile( seq_file1 ):
		sys.exit( "ERROR: input file1 not detected!" )
	if not os.path.isfile( seq_file2 ):
		sys.exit( "ERROR: input file2 not detected!" )
	
	blast_type = parameters[ parameters.index( '--seq_type' )+1 ]
	
	RBH_file = prefix + "RBH_file.txt"
	
	seq_file1_db = prefix + "seq_file1_db"
	seq_file2_db= prefix + "seq_file2_db"
	
	seq_file1_blast_result_file = prefix + "seq_file1_blast_result_file.txt"
	seq_file2_blast_result_file = prefix + "seq_file2_blast_result_file.txt"
	
	# --- identify RBHs --- #
	print "db construction ... please wait!"
	if blast_type == "nucl":
		db_command1 = "makeblastdb -in " + seq_file1 + " -out " + seq_file1_db + " -dbtype 'nucl' -parse_seqids"
		db_command2 = "makeblastdb -in " + seq_file2 + " -out " + seq_file2_db + " -dbtype 'nucl' -parse_seqids"
	else:
		db_command1 = "makeblastdb -in " + seq_file1 + " -out " + seq_file1_db + " -dbtype 'prot' -parse_seqids"
		db_command2 = "makeblastdb -in " + seq_file2 + " -out " + seq_file2_db + " -dbtype 'prot' -parse_seqids"
	
	os.popen( db_command1 )
	os.popen( db_command2 )
	print " ... done."
	
	
	print "BLAST1 ... "
	if blast_type == "nucl":
		BLAST_command1 = "blastn -query " + seq_file1 + " -db " + seq_file2_db + " -out " + seq_file1_blast_result_file + " -outfmt 6 -evalue 0.0001 -max_target_seqs 5 -num_threads 8"
	else:
		BLAST_command1 = "blastp -query " + seq_file1 + " -db " + seq_file2_db + " -out " + seq_file1_blast_result_file + " -outfmt 6 -evalue 0.0001 -max_target_seqs 5 -num_threads 8"
	os.popen( BLAST_command1 )
	
	print "BLAST2 ... "
	if blast_type == "nucl":
		BLAST_command2 = "blastn -query " + seq_file2 + " -db " + seq_file1_db + " -out " + seq_file2_blast_result_file + " -outfmt 6 -evalue 0.0001 -max_target_seqs 5 -num_threads 8"
	else:
		BLAST_command2 = "blastp -query " + seq_file2 + " -db " + seq_file1_db + " -out " + seq_file2_blast_result_file + " -outfmt 6 -evalue 0.0001 -max_target_seqs 5 -num_threads 8"
	os.popen( BLAST_command2 )
	
	
	print "analyzing BLAST results ... please wait!"
	data1 = load_results_from_BLAST_result_file( seq_file1_blast_result_file )
	data2 = load_results_from_BLAST_result_file( seq_file2_blast_result_file )
	seq_IDs_of_interest = compare_datasets( data1, data2, RBH_file )
	
	# --- collect sequences of RBH set in new file --- #
	seq_of_interest_file = prefix + "seq_of_interest_file.fasta"
	sequences = load_multiple_fasta_file( seq_file2 )
	write_seqs_of_interest_into_new_file( seq_of_interest_file, seq_IDs_of_interest, sequences )


if __name__ == '__main__':
	
	if '--prefix' in sys.argv and '--input1' in sys.argv and '--input2' in sys.argv and '--seq_type' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"
