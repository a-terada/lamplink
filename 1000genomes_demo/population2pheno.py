#!/usr/bin/env python

"""
Generate a PED formatted file by regarding
the user-specified population as cases and
others as controls. 
"""

__author__ = "Aika TERADA"

import sys

PANEL_FILE = "integrated_call_samples.20130502.ALL.ped"

def readPanelFile( pos_pheno ):
	sample2anno_dict = {}
	f = open( PANEL_FILE, 'r' )
	line = ""; line_num = 0
	for line in f:
		line_num += 1
		if line_num < 2:
			continue
		s = line[:-1].split( '\t' )
		sample = s[1]; population = 1; gender = 0
		
		if s[6] == pos_pheno:
			population = 2
			
		gender = int( s[4] )
		if population > 0:
			sample2anno_dict[ sample ] = tuple( [gender, population] )
	f.close()
	return sample2anno_dict


def convertPheno( in_filename, sample2anno_dict ):
	f = open( in_filename, 'r' )
	line = ""
	for line in f:
		s = line[:-1].split( ' ' )
		for i in xrange( 0, 4 ):
			sys.stdout.write( "%s\t" % s[i] )
		# output gender and phenotype
		gender, population = sample2anno_dict[ s[0] ]
		sys.stdout.write( "%d\t%d" % ( gender, population ) )
		# output genotype
		for i in s[6:]:
			sys.stdout.write( "\t%s" % i )
		sys.stdout.write( "\n" )
	f.close()

def run( in_filename, pos_pheno ):
	try:
		sample2anno_dict = readPanelFile( pos_pheno )
		convertPheno( in_filename, sample2anno_dict )
		
	except IOError, e:
		sys.stderr.write( "Error: %e\n" )
	
	

if __name__ == "__main__":
	if len( sys.argv ) < 3:
		sys.stderr.write( "input ped file, positive_population.\n" )
		sys.exit()

	run( sys.argv[1], sys.argv[2] )
