/*
 * Copyright (c) 2013, LAMP development team
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the LAMP development team nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * File:   FastWYCore.cpp
 * Author: eiji
 *
 * Created on 2016/10/05, 11:26
 */

#include "FastWYCore.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <sstream>

#include <boost/version.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

#define RANDOM_EXT	0
#ifdef RANDOM_EXT
#include <fstream>
std::istream*	fin = nullptr;
#endif

double
duration_time(
		std::chrono::high_resolution_clock::time_point	begin_time,
		std::chrono::high_resolution_clock::time_point	end_time = std::chrono::high_resolution_clock::now()) {
	return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count() / 1000000000.0;
}

/**
 * Constructor
 */
FastWYCore::FastWYCore() : LampCore() {
	// initialize a pseudo-random number generator engine.
	std::random_device  seed_gen;
	mt_gen.seed(seed_gen());
}

/**
 * Destructor
 */
FastWYCore::~FastWYCore() {
}

/**
 * Run FastWY.
 *
 * @param transaction_file The file includes associations between TFs and genes.
 *                          Each line indicates a gene.
 *                          If gene is targeted by the TF, then value is 1, otherwise 0.
 * @param flag_file Each line indicates a gene. The column1 is gene name.
 *                   If gene has the feature, the column2 is 1. The other case 0.
 * @param threshold The statistical significance threshold.
 * @param k The number of permuted datasets.
 * @param set_method The procedure name for calibration p-value (fisher/u_test/chi).
 * @param max_comb the maximal size which the largest combination size in tests set.
 * @param log_file File name for logging.
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 * @param output if true then output result to stdout
 */
void
FastWYCore::run(
		const std::string&	transaction_file,
		const std::string&	flag_file,
		double				threshold,
		int					k,
		const std::string&	set_method,
		int					max_comb,
		const std::string&	log_file,
		int					alternative,
		bool				output) {
	clean();

	this->transaction_file = transaction_file;
	this->flag_file = flag_file;
	this->trans4lcm = transaction_file + ".4lcm53"; // the filename for outputting logs
	this->max_comb = max_comb;
	this->permute_num = k;
	this->threshold = threshold;
	this->set_method = set_method;
	this->alternative = alternative;
	
	// read 2 files and get transaction list
	std::cerr << "Read input files ..." << std::endl;
	try {
		readFile.readFiles(transaction_file, flag_file, ',');
		// If the alternative hypothesis is 'less',
		// the positive and negative of observe values are reversed,
		// and conduct the identical procedure to 'greater'.
		if (alternative < 0) {
			reverseValue(readFile.getTransaction_list(), set_method);
		}
		if ((int)readFile.getColumnid2name().size() <= max_comb) {
			max_comb = -1;
		}
	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to read input files.");
	}

	// run multiple test
	try {
		FILE*	fp_log = fopen(log_file.c_str(), "w");
		if (fp_log == nullptr) {
			throw std::string("Can't open file : " + log_file);
		}
		boost::iostreams::stream<boost::iostreams::file_descriptor_sink> outlog;
#if (BOOST_VERSION >= 104900)
		outlog.open(fileno(fp_log), boost::iostreams::close_handle);
#else
		outlog.open(fileno(fp_log), true);
#endif

		auto	starttime = std::chrono::high_resolution_clock::now();

		// generate null distribution
		std::cerr << "Calculate the minimum p-value distribution using the permutation test ..." << std::endl;
		outlog << "Calculate the minimum p-value distribution using the permutation test ..." << std::endl;
		max_lambda = maxLambda(readFile.getTransaction_list());
		fre_pattern = new LCM(max_lambda, fileno(fp_log));
		std::vector<min_p_t>	min_p_list;
		generateMinPDist(outlog, min_p_list);

		// adjusted significance level
		outlog << "Adjust significance level ..." << std::endl;
		std::vector<min_p_t>	sorted_min_p_list;
		double	adjusted_threshold = adjustedThreshold(min_p_list, sorted_min_p_list);
		outlog << "Adjusted significance level: " << adjusted_threshold << std::endl;

		auto	correction_term_time = std::chrono::high_resolution_clock::now();

		// enumerate combination whose P-value up to adjusted threshold
		outlog << "Calculate the p-values in the given data set ..." << std::endl;
		std::vector<enrich_t>	enrich_list;
		double					time_enumerate_freq;
		double					time_enumerate_total;
		enumerateSigComb(adjusted_threshold, outlog, enrich_list, time_enumerate_freq, time_enumerate_total);

		auto	finish_test_time = std::chrono::high_resolution_clock::now();

		// set the significant combinations to result
		setResult(enrich_list, adjusted_threshold, sorted_min_p_list);

		if (output) {
			// output the significant combinations
			outputResult();

			// output time cost
			std::cout << boost::format("Time (sec.): Computing correction factor %.3f") %
										duration_time(starttime, correction_term_time) <<
						boost::format(", Enumerating significant combinations %.3f") %
										time_enumerate_total <<
						boost::format(", Total %.3f") %
										duration_time(starttime, finish_test_time) <<
						std::endl;

			// output the minimum P-values
			outputMinP(min_p_list);
		}

		outlog.close();

	} catch (std::string &msg) {
		throw msg;
	} catch (...) {
		throw std::string("Error: An unexpected error occurred while trying to test.");
	}
}

/**
 * Generate the permuted values dataset.
 *
 * @param org2shuffled_list [out] the mapping from transaction ID from the row dataset to shuffled dataset.
 **/
void
FastWYCore::permuteData(
		std::vector<int>&	org2shuffled_list) {
	const auto&			transaction_list = readFile.getTransaction_list();
	std::vector<int>    random_index_list(transaction_list.size());
	std::iota(random_index_list.begin(), random_index_list.end(), 0);
	std::shuffle(random_index_list.begin(), random_index_list.end(), mt_gen);

#ifdef RANDOM_EXT
	fin->read((char*)random_index_list.data(), sizeof(int) * random_index_list.size());
#endif

	org2shuffled_list.resize(transaction_list.size());
	
	for (int i = 0; i < transaction_list.size(); ++i) {
		org2shuffled_list[transaction_list[random_index_list[i]]->getID()] = i;
	}
}

/**
 * Calculate the minimum p-value in the transaction_list shuffled by org2shuffled_list.
 *
 * @param org2shuffled_list [in] the mapping from transaction ID from the row dataset to shuffled dataset.
 * @param low_sup [out] minimum support.
 * @param freq_time [out] time to construct apriori.
 * @return minimum p-value.
 **/
double
FastWYCore::calculateMinimumPValue(
		const std::vector<int>&	org2shuffled_list,
		int&					low_sup,
		double&					freq_time) {
	double	min_p = 1.0;
	const Node::itemset_t*	min_p_pattern = nullptr; // the minimum p-value of the permuted set
	bool	flag = true;

	low_sup = fre_pattern->getMax_support();
	freq_time = 0;
	while (flag) {
		auto	starttime = std::chrono::high_resolution_clock::now(); // time to construct apriori
		fre_pattern->frequentPatterns(trans4lcm, low_sup, max_comb); // construct frequent patterns
		double	bound = calBound(low_sup);
		const auto&	cal_list = fre_pattern->getItemsetList(low_sup); // Itemset calculated its P-value
		freq_time += duration_time(starttime); // time to construct apriori
		std::vector<int>	flag_transaction_list; // transaction list which has all items in itemset.
		for (auto cal: cal_list) {
			flag_transaction_list.clear();
			flag_transaction_list.reserve(cal->tran_list->size());
			for (auto t: *(cal->tran_list)) {
				int	shuffled_id = org2shuffled_list[t];
				flag_transaction_list.push_back(shuffled_id);
			}
			double	score, stat_value;
			double	p = func_f->calPValue(flag_transaction_list, score, stat_value);
			if (p < min_p) {
				min_p = p;
				min_p_pattern = cal;
			}
		}
		// If the minimum p-value is less than the lower bound of P-value, finish the calculation.
		if ((min_p < bound) || (low_sup <= 1)) {
			flag = false;
		// If the minimum p-value is over than MASL, the minimum support is small and repeat the calculation.
		} else {
			--low_sup;
		}
	}

	return min_p;
}

/**
 * Generate a probability distribution of the minimum P-value using permuted datasets
 *
 * @param outlog [in] file object to output logs.
 * @param min_p_list [out] set of minimum p-values using permuted data.
 **/
void
FastWYCore::generateMinPDist(
		std::ostream&			outlog,
		std::vector<min_p_t>&	min_p_list) {
	auto	starttime = std::chrono::high_resolution_clock::now();

#ifdef RANDOM_EXT
	std::ostringstream	str;
	str << "random_" << permute_num << "x" << readFile.getTransaction_list().size() << ".dat";
	std::ifstream	in(str.str(), std::ifstream::binary);
	printf("%s\n", str.str().c_str());
	fin = &in;
#endif

	// Initialize the apriori and functinos using LAMP.
	runMultTest(readFile.getTransaction_list(),
				trans4lcm, threshold, set_method,
				max_comb, outlog, alternative, max_lambda);

	// calculate the set of minimum p-values using permuted data
	min_p_list.resize(permute_num);

	std::vector<int>	org2shuffled_list;
	// estimate the probability distribution of the minimum p-value using permuted datasets.
	for (int i = 0; i < permute_num; ++i) {
		permuteData(org2shuffled_list); // generate the permuted dataset.
		int		low_sup;
		double	freq_time, per_time;
		func_f->resetCalTime();
		double	min_p = calculateMinimumPValue(org2shuffled_list, low_sup, freq_time);
		int	total = fre_pattern->getTotal(low_sup);
		int	calTime = func_f->getCalTime();

		per_time = duration_time(starttime);

		min_p_list[i].min_p = min_p;
		min_p_list[i].low_sup = low_sup;
		min_p_list[i].total = total;
		min_p_list[i].freq_time = freq_time;
		min_p_list[i].per_time = per_time;
		min_p_list[i].cal_time = calTime;

		outlog << "[permute " << i <<
					"] minP " << min_p <<
					", minSupport " << low_sup <<
					", totalTest " << total <<
					", freqTime " << freq_time <<
					", totalTime " << per_time <<
					", #ofPvalue " << calTime << '\n';

		starttime = std::chrono::high_resolution_clock::now();
	}

#ifdef RANDOM_EXT
	fin = nullptr;
	in.close();
#endif
}

/**
 * Calculate the adjusted significance level
 *
 * @param min_p_list [in] the list of minimum P-values used by FastWY
 * @param sorted_min_p_list [out] the sorted list of minimum P-values used by FastWY
 * @return the adjusted significance level
 **/
double
FastWYCore::adjustedThreshold(
		const std::vector<min_p_t>&	min_p_list,
		std::vector<min_p_t>&		sorted_min_p_list) {
	// calculate the adjusted significance level
	int	min_p_index = std::max((int)std::floor(permute_num * threshold) - 1, 0);
	sorted_min_p_list = min_p_list;
	std::sort(sorted_min_p_list.begin(), sorted_min_p_list.end(),
			[](const min_p_t& x, const min_p_t& y) { return x.min_p < y.min_p; });
	double	adjusted_threshold = sorted_min_p_list[min_p_index].min_p;

	if (min_p_index + 1 >= min_p_list.size()) {
		return adjusted_threshold;
	}

	while (min_p_index > 0 &&
			adjusted_threshold == sorted_min_p_list[min_p_index + 1].min_p) {
		adjusted_threshold = sorted_min_p_list[--min_p_index].min_p;
	}
	return adjusted_threshold;
}

/**
 * Enumerate significant combinations (P-value <= adjusted threshold).
 *
 * @param adjusted_threshold [in] adjusted threshold for P-value
 * @param outlog [in] file object to output logs
 * @param enrich_lst [out] the list contains significant combinations
 * @param freq_time [out] the time to construct apriori
 * @param time_enumerate_total [out] the time to enumerate significant combinations
 **/
void
FastWYCore::enumerateSigComb(
		double					adjusted_threshold,
		std::ostream&			outlog,
		std::vector<enrich_t>&	enrich_lst,
		double&					freq_time,
		double&					time_enumerate_total) {
	// test the raw (non-permuted) dataset
	int		i = 0;
	bool	flag = true;
	int		low_sup = fre_pattern->getMax_support();

	auto	starttime = std::chrono::high_resolution_clock::now();

	std::vector<int>	flag_transaction_list; // transaction list which has all items in itemset.

	freq_time = 0;
	while (flag) {
		auto	freq_starttime = std::chrono::high_resolution_clock::now();
		fre_pattern->frequentPatterns(trans4lcm, low_sup, max_comb); // construct frequent patterns
		double	bound = calBound(low_sup);
		const auto&	cal_list = fre_pattern->getItemsetList(low_sup); // Itemset calculated its P-value
		freq_time += duration_time(freq_starttime);
		func_f->resetCalTime();
		for (auto cal: cal_list) {
			++i;
			outlog << "--- testing " << i << ": ";
			flag_transaction_list.clear();
			flag_transaction_list.reserve(cal->tran_list->size());
			for (auto t: *(cal->tran_list)) {
				flag_transaction_list.push_back(t);
			}
			double	score, stat_value;
			double	p = func_f->calPValue(flag_transaction_list, score, stat_value);
			outlog << "p " << p << ", stat_score " << score << '\n';
			if (p <= adjusted_threshold) {
				enrich_lst.push_back({cal->item_list, p, low_sup, score});
			}
		}
		// If the minimum p-value is less than MASL, finish the calculation.
		if ((adjusted_threshold < bound) || (low_sup <= 1)) {
			flag = false;
		// If the minimum p-value is over than MASL, the minimum support is small and repeat the calculation.
		} else {
			--low_sup;
		}
	}

	time_enumerate_total = duration_time(starttime);
}

/**
 * Adjusted threshold.
 *
 * @param pvalue [in] a raw P-value
 * @param sorted_min_p_list [in]
 * @param start_index [in,out] pvalue is compared with values from the start_index-th and the subsequent P-values of the sorted_min_p_list
 * @return an adjusted P-value
 **/
double
FastWYCore::adjustPval(
		double						pvalue,
		const std::vector<min_p_t>&	sorted_min_p_list,
		int&						start_index) {
	double	adj_pvalue;
	for (; start_index < sorted_min_p_list.size(); ++start_index) {
		double	min_p = sorted_min_p_list[start_index].min_p;
		if (min_p > pvalue) {
			break;
		}
	}
	if (pvalue < sorted_min_p_list[0].min_p) {
		adj_pvalue = -1;
	} else if (pvalue == sorted_min_p_list[start_index-1].min_p) {
		adj_pvalue = (double)start_index / (double)sorted_min_p_list.size();
	} else {
		adj_pvalue = (double)(start_index + 1) / (double)sorted_min_p_list.size();
	}

	return adj_pvalue;
}

/**
 * set results.
 *
 * @param enrich_list [in] list contains significant combinations
 * @param adjusted_threshold [in] adjusted threshold for P-value
 * @param sorted_min_p_list [in] the list of minimum P-values sorted used by FastWY. This list is sorted by P-values
 **/
void
FastWYCore::setResult(
		const std::vector<enrich_t>&	enrich_list,
		double							adjusted_threshold,
		const std::vector<min_p_t>&		sorted_min_p_list) {
	auto&	transaction_list = readFile.getTransaction_list();
	auto&	columnid2name = readFile.getColumnid2name();
	bool	is_utest = set_method.compare("u_test") == 0;

	std::ostringstream	str;
	result_info.clear();
	
	// output setting
	str << "# FastWY ver. " << getVersion() << std::endl;
	str << "# item-file: " << transaction_file << std::endl;
	str << "# value-file: " << flag_file << std::endl;
	str << "# significance-level: " << threshold <<
				", # of permuted datasets: " << permute_num << std::endl;
	str << "# P-value computing procedure: " << set_method;
	if (alternative > 0) {
		str << " (greater)" << std::endl;
	} else if (alternative < 0) {
		str << " (less)" << std::endl;
	} else {
		str << " (two sided)" << std::endl;
	}
	str << "# # of tested elements: " << columnid2name.size() <<
				", # of samples: " << transaction_list.size();
	if (!is_utest) {
		str << ", # of positive samples: " << func_f->getN1();
	}
	str << std::endl;
	str << boost::format("# Adjusted significance level: %.5g") % adjusted_threshold << std::endl;
	str << "# # of significant combinations: " << enrich_list.size() << std::endl;
	result_info = str.str();

	if (!enrich_list.empty()) {
		auto	copy_enrich_list = enrich_list;
		std::sort(copy_enrich_list.begin(),
					copy_enrich_list.end(),
					[](const enrich_t& x, const enrich_t& y)
							{ return x.p < y.p; });
		int		minp_index = 0;
		result_list.resize(copy_enrich_list.size());
		auto	iter = result_list.begin();
		for (const auto& l: copy_enrich_list) {
			iter->p = l.p;
			iter->adjusted_p = adjustPval(l.p, sorted_min_p_list, minp_index);
			iter->comb.clear();
			for (auto i = l.item_set->begin(); i != l.item_set->end(); ++i) {
				iter->comb.push_back(*columnid2name[*i - 1]);
			}
			iter->arity = (int)l.item_set->size();
			iter->target_rows = l.len;
			iter->stat_score = l.stat_score;
			++iter;
		}
	}
}

/**
 * Output results.
 **/
void
FastWYCore::outputResult() {
	bool	is_utest = set_method.compare("u_test") == 0;

	// output setting
	std::cout << result_info;

	// output header
	if (!result_list.empty()) {
		std::cout << "Rank\tRaw p-value\tAdjusted p-value\tCombination\tArity\t# of target rows\t";
		if (is_utest) {
			std::cout << "z-score" << std::endl;
		} else {
			std::cout << "# of positives in the targets" << std::endl;
		}
		int		rank = 0;
		double	smalllest_p = 1.0 / (double)permute_num;
		for (const auto& l : result_list) {
			++rank;
			std::cout << boost::format("%d\t%.5g\t") % rank % l.p;
			if (l.adjusted_p < 0) {
				std::cout << boost::format("< %.5g\t") % smalllest_p;
			} else {
				std::cout << boost::format("%.5g\t") % l.adjusted_p;
			}
			std::cout << l.comb.front();
			for (auto i = l.comb.begin() + 1; i != l.comb.end(); ++i) {
				std::cout << "," << *i;
			}
			std::cout << "\t" << l.arity << "\t" << l.target_rows << "\t";
			if (is_utest) {
				std::cout << boost::format("%.5g") % l.stat_score << '\n';
			} else {
				std::cout << (int)l.stat_score << '\n';
			}
		}
	}
}

/**
 * Output minimum P-value distribution.
 *
 * @param min_p_list [in] the list of minimum P-values used by FastWY
 */
void
FastWYCore::outputMinP(
		const std::vector<min_p_t>&	min_p_list) {
	// output minP distribution
	std::cout << "--- minimum P-values ---" << std::endl;
	std::cout << "[id]\tminP\tminSupport\t#ofCandComb\tfreqTime\ttotalTime\t#ofPvalueCalculation" << std::endl;
	int	i = 0;
	for (const auto& min_p: min_p_list) {
		std::cout << "[" << i++ << "]\t" <<
					boost::format("%.5g") % min_p.min_p << "\t" <<
					min_p.low_sup << "\t" <<
					min_p.total << "\t" <<
					boost::format("%.5g") % min_p.freq_time << "\t" <<
					boost::format("%.5g") % min_p.per_time << "\t" <<
					min_p.cal_time << '\n';
	}
}
