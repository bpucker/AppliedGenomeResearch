import re, os
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# --- end of imports --- #

def load_expression_values( filename, candidate_genes ):
	"""! @brief load all expression values from featureCounts result file and return normalized gene expression """
	
	# --- construct emty data structure for gene expression analysis --- #
	expression_data = {}
	for ID in candidate_genes:
		expression_data.update( { ID: 0 } )
	
	# --- load data from file --- #
	total_mapped_reads = 0
	with open( filename, "r" ) as f:
		f.readline()	#header1
		f.readline()	#header2
		line = f.readline()
		while line:
			ID = line.strip().split('\t')[0]
			current = int( line.strip().split('\t')[-1] )
			expression_data.update( { ID: current } )
			total_mapped_reads += current
			line = f.readline()
	
	# --- normalize values --- #
	total_mapped_reads = float( total_mapped_reads )
	for key in expression_data.keys():
		new_value = expression_data[ key ] / (total_mapped_reads / 1000000)	#TPM = Tags Per Million reads calculation
		del expression_data[ key ]
		expression_data.update( { key: new_value } )
		
	return expression_data


def construct_data_output_file( data, candidate_genes, outputfile, cutoff ):
	"""! @brief write expression values of all candidate genes into output file """
	
	datamatrix = []
	genes = []
	with open( outputfile, "w" ) as out:
		new_line = [ "gene" ]
		for each in data:
			new_line.append( each['id'] )
		tissues = new_line[1:]
		out.write( "\t".join( new_line ) + '\n' )
		for gene in sorted( candidate_genes ):
			new_line = [ gene ]
			for tissue in data:
				new_line.append( tissue['values'][ gene ] )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )
			if sum( new_line[1:] ) > cutoff:
				datamatrix.append( new_line[1:] )
				genes.append( gene )
	return genes, tissues, datamatrix


def construct_heatmap( datamatrix, genes, tissues, heatmap_file ):
	"""! @brief construct heatmap from given data matrix """
	
	print "number of genes for heatmap construction: " + str( len( genes ) )
	labels = []
	for gene in genes:
		labels.append( re.findall( "DN\d+", gene )[0] )
	
	df = DataFrame( datamatrix, index=labels, columns=tissues)
	
	sns.heatmap( df, linewidths=0.5, annot=True, annot_kws={'fontsize':5} )	#fmt="d"
	plt.subplots_adjust( left=0.2 )
	plt.savefig( heatmap_file, dpi=300  )


if __name__ == '__main__':
	
	candidate_gene_file = "<ADD_YOUR_PATH>/candidate_seqs.txt"	#file with up to 5 selected contigs/unigenes
	
	data_files = [ 	{ 'id': "WT", "file": "<YOUR_FILE_NAME_WT>.count_table" },
							{ 'id': "mut", "file": "<YOUR_FILE_NAME_MUT>.count_table" }
					]
	cutoff = 0 #minimal cumulative expression of gene in all tissues combined
	prefix = "<YOUR_OUTPUT_DIRECTORY>"
	
	
	if prefix[-1] != '/':
		prefix += "/"
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	# --- load data --- #
	candidate_seq_names = [ ]
	with open( candidate_gene_file, "r" ) as f:
		line = f.readline()
		while line:
			candidate_seq_names.append( line.strip() )
			line = f.readline()
	
	expression_data = []
	for entry in data_files:
		expression_data.append( { 'id': entry['id'], 'values': load_expression_values( entry['file'], candidate_seq_names ) } )
	
	# --- write normalized MADS box gene expression values into data output file --- #
	outputfile = prefix + "normalized_values.txt"
	genes, tissues, datamatrix = construct_data_output_file( expression_data, candidate_seq_names, outputfile, cutoff )
	
	heatmap_file = prefix + "heatmap.png"
	construct_heatmap( datamatrix, genes, tissues, heatmap_file )
	
	print "all done!"
