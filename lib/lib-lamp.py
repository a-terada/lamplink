#!/usr/bin/env python

"""
Copyright (c) 2013-2015, LAMP development team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP development team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import sys, os.path, time, datetime
import transaction, lamp
import readFile
import frepattern.frequentPatterns as frequentPatterns
from optparse import OptionParser

import functions.functionsSuper as fs
import functions.functions4fisher as functions4fisher
import functions.functions4u_test as functions4u_test
import functions.functions4chi as functions4chi

set_opts = ("fisher", "u_test", "chi") # methods which used each test

__version__ = lamp.__version__


def makeResult( transaction_file, flag_file, threshold, set_method, max_comb, columnid2name, lam_star, k, \
	enrich_lst, transaction_list, func_f, alternative ):
	flag_size = -1
	if not set_method == "u_test":
		flag_size = func_f.getN1()
	result_strings=""
	
	# output setting
	result_strings="# LAMP ver. %s\t" % lamp.__version__
	result_strings=result_strings+"# item-file: %s\t" % (transaction_file)
	result_strings=result_strings+"# value-file: %s\t" % (flag_file)
	result_strings=result_strings+"# significance-level: %s\t" % threshold
	result_strings=result_strings+"# P-value computing procedure: %s" % (set_method)
	if alternative > 0:
		result_strings=result_strings+" (greater)\t"
	elif alternative < 0:
		result_strings=result_strings+" (less)\t"
	else:
		result_strings=result_strings+" (two.sided)\t"
	result_strings=result_strings+"# # of tested elements: %d, # of samples: %d" % ( len(columnid2name), len(transaction_list) )
	if flag_size > 0:
		result_strings=result_strings+", # of positive samples: %d" % flag_size
	result_strings=result_strings+"\t"
	result_strings=result_strings+"# Adjusted significance level: %.5g, " % (threshold/k) 
	result_strings=result_strings+"Correction factor: " + str(k) + " (# of target rows >= " + str(lam_star) + ")\t" 
	result_strings=result_strings+"# # of significant combinations: " + str(len(enrich_lst)) + "\t"
	# output header
	result_strings=result_strings+"!"
	if len(enrich_lst) > 0:
		enrich_lst.sort(lambda x,y:cmp(x[1], y[1]))
		rank = 0
		for l in enrich_lst:
			rank = rank + 1
			if flag_size >0:
				x=flag_size-l[3]
				y=l[2]-l[3]
				z=len(transaction_list)-(l[3]+x+y)
				result_strings=result_strings+"%d\t%d\t%d\t%d\t" % (l[3],x,y,z)
			result_strings=result_strings+"%d\t%.5g\t%.5g\t" % (rank, l[1], k*l[1])
			out_column = ""
			for i in l[0]:
				out_column = out_column + columnid2name[i-1] + ","
			result_strings=result_strings+"%s"%out_column[:-1]+"/"
	return result_strings

def version():
	return __version__

##
# Run multiple test.
# itemset_file: The file includes associations between TFs and genes.
#     Each line indicates a gene.
#     If gene is targeted by the TF, then value is 1, otherwise 0.
# flag_file: Each line indicates a gene. The column1 is gene name.
#     If gene has the feature, the column2 is 1. The other case 0.
# threshold: The statistical significance threshold.
# set_method: The procedure name for calibration p-value (fisher/u_test).
# max_comb: the maximal size which the largest combination size in tests set.
# alternative: hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
##
def run(transaction_file, flag_file, threshold, set_method, lcm_path, max_comb, log_file, alternative):
	# read 2 files and get transaction list
	sys.stderr.write( "Read input files ...\n" )
	transaction_list = set()
	try:
		transaction_list, columnid2name = readFile.readFiles(transaction_file, flag_file, ',')
		max_comb = lamp.convertMaxComb( max_comb, len(columnid2name) )
	except ValueError, e:
		return
	except KeyError, e:
		return
	
	# run multiple test
	transaction4lcm53 = transaction_file + ".4lcm53"
	# run
	try:
		outlog = open( log_file, 'w' )

		starttime = time.time()
		sys.stderr.write( "Compute the optimal correction factor ..." )
		fre_pattern, lam_star, max_lambda, correction_term_time, func_f \
					 = lamp.runMultTest(transaction_list, transaction4lcm53, threshold, \
                                        set_method, lcm_path, max_comb, outlog, alternative)
		k = fre_pattern.getTotal( lam_star )
		sys.stderr.write( " %s\n" % k )
		sys.stderr.write( "Compute P-values of testable combinations ...\n" )
		enrich_lst, finish_test_time \
					= lamp.fwerControll(transaction_list, fre_pattern, lam_star, max_lambda, \
								   threshold, func_f, columnid2name, outlog)
		
		outlog.close()
	except IOError, e:
		outlog.close()

	sys.stderr.write( "Output results ...\n" )
	# output result
	result = makeResult(transaction_file, flag_file, threshold, set_method, max_comb, \
                                  columnid2name, lam_star, k, enrich_lst, transaction_list, func_f, alternative )
	#print "result::::",result
	# output time cost
	sys.stdout.write("Time (sec.): Computing correction factor %.3f, Enumerating significant combinations %.3f, Total %.3f\n" \
					 % (correction_term_time-starttime, finish_test_time - correction_term_time, finish_test_time - starttime))

	return enrich_lst, k, lam_star, columnid2name,result


def lib_lamp(command_str):
	print "\n--------------LAMP-START------------------\n"
	#print command_str
	usage = "usage: %prog [options] transaction_file value_file significance_probability"
	p = OptionParser(usage = usage, version = "%s" % __version__)
	p.add_option('-p', '--pvalue', dest = "pvalue_procedure", help = "Choose the p-value calculation procedure from 'fiehser' (Fisher's exact test), 'chi' (Chi-square test) or 'u_test' (Mann-Whitney's U-test)")

	p.add_option('--lcm', dest = "lcm_path", default = "./lcm53/lcm", \
				 help = "Set LCM program path if you do not use ./lcm53/lcm")

	p.add_option('--max_comb', dest = "max_comb", default = "all", \
				 help = "Set the maximum size of combination to be tested.")
	
	p.add_option('-e', dest = "log_filename", default = "", help = "The file name to output log.\n")
    
	p.add_option('--alternative', dest = "alternative", default = "greater", help = "Indicate which alternative hypothesis is used. Select \"greater\", \"less\" or \"two.sided\"\n, and the default is \"greater\".")
	
	
	opts, args = p.parse_args(command_str.split())
	# check arguments
	if len(args) != 3:
		sys.stderr.write("Error: input [target-file], [expression-file] and [significance-level].\n")
		sys.exit()

	opts.max_comb = opts.max_comb.lower()
	if not opts.max_comb == "all":
		if (opts.max_comb.isdigit()):
			opts.max_comb = int(opts.max_comb)
		else:
			sys.stderr.write("Error: max_comb must be an integer value or all.\n")
			sys.exit()
		
	# check the file exist.
	if not os.path.isfile(args[0]):
		sys.stderr.write("IOError: No such file: \'" + args[0] + "\'\n")
		sys.exit()
	if not os.path.isfile(args[1]):
		sys.stderr.write("IOError: No such file: \'" + args[1] + "\'\n")
		sys.exit()
	try:
		sig_pro = float(args[2])
		if (sig_pro < 0) or (sig_pro > 1):
			raise ValueError
	except ValueError:
		sys.stderr.write("Error: significance probability must be a float value from 0.0 to 1.0.\n")
		sys.exit()

	# check the value of alternative hypothesis
	if opts.alternative == "greater":
		opts.alternative = 1
	elif opts.alternative == "less":
		opts.alternative = -1
	elif opts.alternative == "two.sided":
		opts.alternative = 0
	else:
		print opts.alternative
		sys.stderr.write( "Error: \"alternative\" should be one of {\"greater\", \"less\", \"two.sided\"}\n" )
		sys.exit()
	
	# change log file
	d = datetime.datetime.today()
	log_file = "lamp_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
	if len(opts.log_filename) > 0:
		log_file = opts.log_filename

	opts.delimiter = ','
	
	transaction_file = args[0]; flag_file = args[1]; threshold = float(args[2])
	enrich_lst, k, lam_star, columnid2name,result_str \
				= run(transaction_file, flag_file, threshold, opts.pvalue_procedure, \
					  opts.lcm_path, opts.max_comb, log_file, opts.alternative)
	return result_str

if __name__ == "__main__":
	lib_lamp( sys.argv[1] )
