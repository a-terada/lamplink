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
 * File:   FastWYCore
 * Author: eiji
 *
 * Created on 2016/10/05, 11:19
 */
#ifndef FASTWYCORE_H
#define FASTWYCORE_H

#include "LampCore.h"

#include <random>
#include <vector>
#include <string>

#define __FASTWY_VER__  "2.0.0" " (LAMP ver." __LAMP_VER__ ")"

/** @file */
/** Execute FastWY function
 */
class FastWYCore : public LampCore {
public:
	/** This struct is fastwy result record.
	 */
	typedef struct {
		double		p;				/**< row p-value */
		double		adjusted_p;		/**< adjusted p-value */
		std::vector<std::string>
					comb;			/**< combination */
		int			arity;			/**< arity */
		int			target_rows;	/**< # of target rows */
		double		stat_score;		/**< # of positives in the targets / z-score */
	} stat_t;

private:
	/** This struct is minimum p-value in the permute_transaction_list
	 */
	typedef struct {
		double	min_p;		/**< minimum p-value */
		int		low_sup;	/**< minimum support */
		int		total;		/**< itemset size that support >= low_sup */
		double	freq_time;	/**< freq_time */
		double	per_time;	/**< per_time */
		int		cal_time;	/**< cal_time */
	} min_p_t;
	
public:
	FastWYCore();
	virtual ~FastWYCore();
	
	void		run(	const std::string&	transaction_file,
						const std::string&	flag_file,
						double				threshold,
						int					k,
						const std::string&	set_method,
						int					max_comb,
						const std::string&	log_file,
						int					alternative,
						bool				output = true);

	/**
	 * get result information text
	 * @return result information text
	 */
	const std::string&
				getResultInfo() const {
					return result_info;
				}
	/**
	 * get result list
	 * @return result list
	 */
	const std::vector<stat_t>&
				getResultList() const {
					return result_list;
				}
	/**
	 * get version number
	 * @return version string
	 */
	static std::string getVersion() { return __FASTWY_VER__; }

protected:
	void		permuteData(
						std::vector<int>&	org2shuffled_list);
	double		calculateMinimumPValue(
						const std::vector<int>&	org2shuffled_list,
						int&					low_sup,
						double&					freq_time);
	void		generateMinPDist(
						std::ostream&			outlog,
						std::vector<min_p_t>&	min_p_list);
	double		adjustedThreshold(
						const std::vector<min_p_t>&	min_p_list,
						std::vector<min_p_t>&		sorted_min_p_list);
	void		enumerateSigComb(
						double					adjusted_threshold,
						std::ostream&			outlog,
						std::vector<enrich_t>&	enrich_lst,
						double&					freq_time,
						double&					time_enumerate_total);
	double		adjustPval(
						double						pvalue,
						const std::vector<min_p_t>&	sorted_min_p_list,
						int&						start_index);
	void		setResult(
						const std::vector<enrich_t>&	enrich_list,
						double							adjusted_threshold,
						const std::vector<min_p_t>&		sorted_min_p_list);
	void		outputResult();
	void		outputMinP(
						const std::vector<min_p_t>&	min_p_list);

	std::string					transaction_file;	/**< The file that includes associations between TFs and genes. */
	std::string					flag_file;			/**< The gene expression file. */
	std::string					trans4lcm;			/**< The file for argument of LCM program. This file is made in LampCore.runMultTest method. */

	int							max_comb;			/**< The maximum arity limit. */
	int							permute_num;		/**< The number of permuted datasets. */
	double						threshold;			/**< The statistical significance threshold. */
	std::string					set_method;			/**< The procedure name for calibration p-value (fisher/chi/u_test). */
	int							alternative;		/**< hypothesis, 1 -> greater, 0 -> two sided, -1 -> less */

	int							max_lambda;			/**< The maximum lambda. */

	std::mt19937_64				mt_gen;				/**< A pseudo-random number generator engine.(a mersenne twister engine) */

	std::string					result_info;		/**< FastWY result information text */
	std::vector<stat_t>			result_list;		/**< FastWY result list */
};

#endif /* FASTWYCORE_H */
