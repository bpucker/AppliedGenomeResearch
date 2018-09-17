### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """"
			python map_annotation.py\n
			--input_file <FULL_PATH_TO_INPUT_FILE>\n
			
			bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
			"""

import re, sys

# --- end of imports --- #

def load_annotation( annotation_file ):
	"""! @brief load all content from annotation file """
	
	mapping_table = {}
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 1:
				mapping_table.update( { parts[0].upper(): parts[1] } )
			line = f.readline()
	print "number of entries in mapping table: " + str( len( mapping_table.keys() ) )
	return mapping_table


def annotate_genes( mapping_table, data_input_file, output_file ):
	"""! @brief reads data from input file and maps annotation to geneID; results are written into output file """
	
	all_hit_genes = []
	with open( output_file, "w" ) as out:
		with open( data_input_file, "r" ) as f:
			line = f.readline()
			while line:
				try:
					gene = re.findall( "AT[12345CM]{1}G\d{5}", line.upper() )[0]
					all_hit_genes.append( gene )
					new_line = [ line.strip(), mapping_table[ gene ] ]
					out.write( '\t'.join( new_line ) + '\n' )
				except:
					new_line = [ line.strip(), "nothing_to_annotate" ]
					out.write( '\t'.join( new_line ) + '\n' )
				line = f.readline()
	print "number of annotated genes: " + str( len( all_hit_genes ) )
	return all_hit_genes


def main( arguments ):
	""""! @brief runs all functions for mapping of annotation """
	
	annotation_file = "/vol/agrcourse/scripts/all_gene_ids_with_function.csv"
	mapping_table = load_annotation( annotation_file )
	
	data_input_file = arguments[ arguments.index( '--input_file' ) + 1 ]
	output_file = data_input_file + "_annotated.txt"
	
	all_hit_genes = annotate_genes( mapping_table, data_input_file, output_file )
	
	print "all done!"
	

if __name__ == "__main__":
	
	if '--input_file' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
