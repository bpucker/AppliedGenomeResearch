### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """
			python construct_DeSeq2_input.py\n
			--sample_sheet <FULL_PATH_TO_TABLE_WITH_INFORMATION_ABOUT_SAMPLES>\n
			--sample_dir <FULL_PATH_TO_DIRECTORY_CONTAINING_COUNT_TABLES>\n
			--output_dir <FULL_PATH_TO_OUTPUT_DIRECTORY(will_be_created_if_necessary)>\n
			
			bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
			"""


import glob, sys, re, os, random
from datetime import date
from operator import itemgetter

# --- end of imports --- #

def load_expression_values( filename ):
	"""! @brief load all expression values from featureCounts result file and return raw counts """
	
	expression_data = {}
	
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			try:
				ID = re.findall( "AT[12345CM]{1}G\d{5}", line )[0]
				current = random.choice( [ 1 ]*100 + [ 1000 ] ) * int( line.strip().split('\t')[-1] )
				expression_data.update( { ID: current } )
			except:
				pass	#print line
			line = f.readline()
		
	return expression_data


def load_sample_information( sample_information_file, data_dir, pseudo_replicates):
	"""! @brief load information about all samples into list of dictionaries """
	
	sample_information = []
	
	with open( sample_information_file, "r" ) as f:
		f.readline()	#header
		line = f.readline()
		while line:
			parts = line.strip().split(',')
			for i in range( pseudo_replicates ):
				sample_information.append( { 	'name': parts[1]+'_'+str(i+1),
												'expression':load_expression_values( data_dir + parts[2] + ".count_table.txt" ),
												'genotype': parts[1]
											} )
			line = f.readline()
	print "number of loaded sample information data sets: " + str( len( sample_information ) )
	return sample_information


def construct_sample_table( sample_information, sample_table ):
	"""! @brief construct sample table for DeSeq2 """
	
	# --- construct sample table --- #
	ordered_valid_samples = sorted( sample_information, key=itemgetter('genotype') )
	with open( sample_table, "w" ) as out:
		out.write( '\t'.join( [ "",  "genotype" ] ) + '\n' )
		for idx, sample in enumerate( ordered_valid_samples ):
			
			new_line = "\t".join( [ sample['name'],			#sample name
									sample['genotype']		#genotype
									] ) + '\n'
			out.write( new_line )
	
	return ordered_valid_samples


def construct_data_matrix( ordered_valid_samples, data_matrix_file ):
	"""! @brief construct the data matrix for DeSeq2 based on featureCount result files """
	
	with open( data_matrix_file, "w" ) as out:
		
		# --- construct header line (all sample IDs in order) --- #
		ordered_valid_ids = []
		for sample in ordered_valid_samples:
			ordered_valid_ids.append( sample['name'] )
		out.write( "\t".join( [ "" ] + ordered_valid_ids ) + '\n' )
		
		# --- construct data body (expression data per sample in one column) --- #
		gene_IDs = sorted( ordered_valid_samples[0]['expression'].keys() )
		for gene in gene_IDs:
			new_line = [ gene ]	#add geneID
			for sample in ordered_valid_samples:
				new_line.append( sample[ 'expression' ][ gene ] )
			out.write( '\t'.join( map( str, new_line ) ) + '\n' )


def main( arguments ):
	"""! @brief takes input and runs all functions """
	
	sample_dir = arguments[ arguments.index( '--sample_dir' ) +1 ]
	sample_information_file = arguments[ arguments.index( '--sample_sheet' ) +1 ]
	
	prefix = arguments[ arguments.index( '--output_dir' ) +1 ]
	
	pseudo_replicates = 3	#1 would prevent the production of pseudo replicates
	
	# --- nothing needs to be changed beyond this point --- #
	if prefix[-1] != '/':
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	sample_table = prefix + "clean_sample_table.txt"
	data_matrix_file = prefix + "clean_data_matrix.txt"
	
	sample_information = load_sample_information( sample_information_file, sample_dir, pseudo_replicates )
	ordered_valid_samples = construct_sample_table( sample_information, sample_table )
	construct_data_matrix( ordered_valid_samples, data_matrix_file )
	
	print "all done!"
	

if __name__ == '__main__':
	
	if '--sample_sheet' in sys.argv and '--sample_dir' in sys.argv and '--output_dir' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
