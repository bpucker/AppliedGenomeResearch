### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###


__usage__ = """
					python fasta2gff.py\n
					--fasta <INPUT_FILE>\n
					--gff3 <OUTPUT_FILE>\n
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import os, re, sys

# --- end of imports --- #


def load_seqs_from_mult_fasta( filename ):
	"""! @brief load all contigs of assembly """
	
	sequences = {}
	
	with open( filename, "r" ) as f:
		header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				sequences.update( { header: seq } )
				header = line.strip()[1:].split(' ')[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	
	print "number of sequences in file: " + str( len( sequences.keys() ) )
	return sequences


def construct_gff( sequences, output_file ):
	"""! @brief construct gff3 file based on available sequences """
	
	with open( output_file, "w" ) as out:
		for key in sorted( sequences.keys() ):
			new_line = [ key, ".", "mRNA", "1", len( sequences[key] ), ".", "+", ".", "ID=" + key + ";Parent=" + key ]
			out.write( "\t".join( map( str, new_line ) ) + '\n' )
			

def main( arguments ):
	"""! @brief calls all functions of this script """
	
	input_fasta = arguments[ arguments.index( '--fasta' )+1 ]
	output_gff =arguments[ arguments.index( '--gff3' )+1 ]
	
	sequences = load_seqs_from_mult_fasta( input_fasta )
	construct_gff( sequences, output_gff )



if __name__ == '__main__':
	
	if '--fasta' in sys.argv and '--gff3' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"
