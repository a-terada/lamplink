
//////////////////////////////////////////////////////////////////
//                                                              //
//          LAMPLINK (c) 2015 LAMP development team             //
//                                                              //
// This code is written to add options                          //
// (--lamp and --lamp-ld-removed) to PLINK v1.07                //
// (http://pngu.mgh.harvard.edu/~purcell/plink/)                //
// for the combinatorial detection with LAMP                    //
// (http://a-terada.github.io/lamp/).                           //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <vector>
#include <cassert>
#include <boost/format.hpp>

#include "plink.h"
#include "fisher.h"
#include "stats.h"
#include "helper.h"
#include "options.h"
#include "crandom.h"
#include "sets.h"
#include "perm.h"
#include "phase.h"

using namespace std;

extern ofstream LOG;

#include "genogroup.h"
#include "haplowindow.h"
#include "lamp.h"
#include "LampCore.h"


vector<string> split(string &str, char delim){
	vector<string> res;
	size_t current = 0, found;
	while((found = str.find_first_of(delim, current)) != string::npos){
		res.push_back(string(str, current, found - current));
		current = found + 1;
	}
	res.push_back(string(str, current, str.size() - current));
	return res;
}
string Replace( string String1, string String2, string String3 )
{
    std::string::size_type  Pos( String1.find( String2 ) );
    
    while( Pos != std::string::npos )
    {
        String1.replace( Pos, String2.length(), String3 );
        Pos = String1.find( String2, Pos + String3.length() );
    }
    return String1;
}

void Plink::LampAssocFull(Perm & perm)
{
	printf("Plink::LampAssocFull\n");
	printLOG("Full-model association for Lamp tests, minimum genotype count: --cell " +
			int2str(par::min_geno_cell) + "\n");

	vector<double> results(nl_all);
	//LampAssoc *outputla = new LampAssoc;
	vector<CSNP*>::iterator s = SNP.begin();
	int l=0;
	while ( s != SNP.end() )
	{ 
		//printf("test:%d\n",l);
		// In adaptive mode, possibly skip this test
		if (par::adaptive_perm && (!perm.snp_test[l]))
		{
			printf("In adaptive mode, possibly skip this test\n");
			s++;
			l++;
			continue;
		}

		int A11=0, A12=0, A22=0;
		int U11=0, U12=0, U22=0;

		/////////////////////////
		//printf("Autosomal or haploid?\n");
		bool X=false, haploid=false;
		if (par::chr_sex[locus[l]->chr]) X=true;
		else if (par::chr_haploid[locus[l]->chr]) haploid=true;

		//printf("Skip haploid markers\n");
		if (haploid)
		{
			s++;
			l++;
			continue;
		}
		/////////////////////////////
		// Iterate over individuals
		vector<bool>::iterator i1 = (*s)->one.begin();
		vector<bool>::iterator i2 = (*s)->two.begin();
		vector<Individual*>::iterator gperson = sample.begin();

		while ( gperson != sample.end() )
		{
			//printf("person!\n");
			// Phenotype for this person (i.e. might be permuted)
			Individual * pperson = (*gperson)->pperson;
			// SNP alleles
			bool s1 = *i1;
			bool s2 = *i2;
			//////////////////////////////////////////////////////////////////////////////////
			if ( ! pperson->missing  ) // if the phenotype is not missing value
			{
				//printf("%d,%d ",s1,s2);
				// Only consider diploid chromosomes
				if ( ! ( X && (*gperson)->sex ) )
				{
					//printf("%d,%d ",s1,s2);
					if ( pperson->aff )     // cases
					{
					  //printf("a,%d,%d \n",s1,s2);
						if ( ! s1 )
						{
							if ( ! s2 )     // Homozyg 00
								A11++;
							else            // Hetero  01
								A12++;
						}
						else                // Homozyg 11 or missing genotype
						  A22++;
						/*
						else if ( s2 )      // Homozyg 11
							A22++;
						*/
					}
					else
					{
					  //printf("u,%d,%d ",s1,s2);
						if ( !s1 )
						{
						  if ( !s2 )       // Homozyg 00
							U11++;
						  else             // Hetero  01
								U12++;
						}
						else                // Homozyg 11 or missing genotype
						  U22++;
						/*
						else if ( s2 )      // Homozyg 11
							U22++;
						*/
					}


				}
			}
			// Next individual
			gperson++;
			i1++;
			i2++;

		}
		//	printf("\n");
		//////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////
		// Calculate association statistics

		double obs_A = A11 + A12 + A22;
		double obs_U = U11 + U12 + U22;
		double obs_T = obs_A + obs_U;

		double obs_1 = 2*(A11+U11) + A12 + U12;
		double obs_2 = 2*(A22+U22) + A12 + U12;

		double obs_11 = A11+U11;
		double obs_12 = A12+U12;
		double obs_22 = A22+U22;

		bool invalid = false;
		if (A11 < par::min_geno_cell ||
				A12 < par::min_geno_cell ||
				A22 < par::min_geno_cell) invalid = true;
		else if (U11 < par::min_geno_cell ||
				U12 < par::min_geno_cell ||
				U22 < par::min_geno_cell) invalid = true;


		if ( par::trend_only )
			invalid = true;
		///////////////////////
		// Cochram-Armitage Trend test
		double CA = ( ( obs_U / obs_T * A12 ) - ( obs_A / obs_T * U12 ) )
			+ 2*( ( obs_U / obs_T * A22 ) - ( obs_A / obs_T * U22 ) ) ;

		double varCA = obs_A * obs_U
			* ( ( obs_T * ( obs_12 + 4*obs_22 )
						- ( obs_12+2*obs_22 ) * ( obs_12+2*obs_22 )  )
					/ (obs_T * obs_T * obs_T  )) ;
		double CA_chisq = (CA*CA) / varCA;
		double CA_p = chiprobP(CA_chisq,1);
		double mult_p, mult_chisq;

		///////////////////////
		// Multiplicative model
		double obs_A1 = 2*A11 + A12;
		double obs_A2 = 2*A22 + A12;
		double obs_U1 = 2*U11 + U12;
		double obs_U2 = 2*U22 + U12;

		if ( par::fisher_test )
		{
			table_t t;
			sizeTable(t,2,2);
			t[0][0] = (int)obs_A1;
			t[1][0] = (int)obs_A2;
			t[0][1] = (int)obs_U1;
			t[1][1] = (int)obs_U2;
			mult_p = fisher(t);

		}
		else
		{

			double exp_A1 = (obs_A * obs_1 ) / obs_T;    // note 2's cancelled for obs_A and obs_T
			double exp_A2 = (obs_A * obs_2 ) / obs_T;    // which are counts of individuals, not
			double exp_U1 = (obs_U * obs_1 ) / obs_T;    // alleles
			double exp_U2 = (obs_U * obs_2 ) / obs_T;

			mult_chisq =  ( ( obs_A1 - exp_A1 ) * ( obs_A1 - exp_A1 ) ) / exp_A1
				+ ( ( obs_A2 - exp_A2 ) * ( obs_A2 - exp_A2 ) ) / exp_A2
				+ ( ( obs_U1 - exp_U1 ) * ( obs_U1 - exp_U1 ) ) / exp_U1
				+ ( ( obs_U2 - exp_U2 ) * ( obs_U2 - exp_U2 ) ) / exp_U2;

			///////////////////////
			// Multiplicative model
			mult_p = chiprobP(mult_chisq,1);

		}


		double gen_p, dom_p, rec_p;
		gen_p = dom_p = rec_p = -9;
		double dom_chisq, rec_chisq, gen_chisq;


		if (!invalid)
		{
			//////////////////////////////////////////////////////////////
			// Standard chi-square test, or Fisher's exact
			if ( par::fisher_test )
			{
				////////////
				// General
				//
				table_t t;                                          
				sizeTable(t,3,2);
				t[0][0] = A11;
				t[1][0] = A12;
				t[2][0] = A22;
				t[0][1] = U11;
				t[1][1] = U12;
				t[2][1] = U22;
				gen_p = fisher(t);

				////////////
				//Dominant
				sizeTable(t,2,2);
				t[0][0] = A11+A12;
				t[1][0] = A22;
				t[0][1] = U11+U12;
				t[1][1] = U22;
				dom_p = fisher(t);

				/////////////
				// Recessive
				sizeTable(t,2,2);
				t[0][0] = A11;
				t[1][0] = A12+A22;
				t[0][1] = U11;
				t[1][1] = U12+U22;
				rec_p = fisher(t);
			}
			else
			{
				///////////////////////
				// General model
				double exp_A11 = (obs_A * obs_11 ) / obs_T;
				double exp_A12 = (obs_A * obs_12 ) / obs_T;
				double exp_A22 = (obs_A * obs_22 ) / obs_T;
				double exp_U11 = (obs_U * obs_11 ) / obs_T;
				double exp_U12 = (obs_U * obs_12 ) / obs_T;
				double exp_U22 = (obs_U * obs_22 ) / obs_T;

				gen_chisq =  ( ( A11 - exp_A11 ) * ( A11 - exp_A11 ) ) / exp_A11
					+ ( ( A12 - exp_A12 ) * ( A12 - exp_A12 ) ) / exp_A12
					+ ( ( A22 - exp_A22 ) * ( A22 - exp_A22 ) ) / exp_A22
					+ ( ( U11 - exp_U11 ) * ( U11 - exp_U11 ) ) / exp_U11
					+ ( ( U12 - exp_U12 ) * ( U12 - exp_U12 ) ) / exp_U12
					+ ( ( U22 - exp_U22 ) * ( U22 - exp_U22 ) ) / exp_U22;
				///////////////////////
				// Dominant (minor allele) (1) model
				dom_chisq =  ( ( (A11+A12) - (exp_A11+exp_A12) ) * ( (A11+A12) - (exp_A11+exp_A12) ) ) / (exp_A11+exp_A12)
					+ ( ( A22 - exp_A22 ) * ( A22 - exp_A22 ) ) / exp_A22
					+ ( ( (U11+U12) - (exp_U11+exp_U12) ) * ( (U11+U12) - (exp_U11+exp_U12) ) ) / (exp_U11+exp_U12)
					+ ( ( U22 - exp_U22 ) * ( U22 - exp_U22 ) ) / exp_U22;
				//////////////////////////////////////
				//Recessive (minor allele) (1) model

				rec_chisq =  ( ( (A22+A12) - (exp_A22+exp_A12) ) * ( (A22+A12) - (exp_A22+exp_A12) ) ) / (exp_A22+exp_A12)
					+ ( ( A11 - exp_A11 ) * ( A11 - exp_A11 ) ) / exp_A11
					+ ( ( (U22+U12) - (exp_U22+exp_U12) ) * ( (U22+U12) - (exp_U22+exp_U12) ) ) / (exp_U22+exp_U12)
					+ ( ( U11 - exp_U11 ) * ( U11 - exp_U11 ) ) / exp_U11;

				//////////////////////////////////
				// p-values and model comparisons
				gen_p = chiprobP(gen_chisq,2);
				dom_p = chiprobP(dom_chisq,1);
				rec_p = chiprobP(rec_chisq,1);
			}
		}
		////////////////////////////////////////////
		// Save best p-value for permutation test
		//////////////////////////
		// Save the desired result
		int best = 0 ;
		if (par::model_perm_best)
		{

			double best_p = mult_p;

			if (!invalid)
			{
				// Skip general model (i.e. just compare ALLELIC, DOM, REC

				//if (gen_p < best_p && gen_p >= 0 ) { best = 2; best_p = gen_p; }   // general
				if (dom_p < best_p && dom_p >= 0 ) { best_p = dom_p; }   // dom
				if (rec_p < best_p && rec_p >= 0 ) { best_p = rec_p; }   // rec
			}

			results[l] = 1-best_p;
		}
		else if ( par::model_perm_gen )
			results[l] = gen_p >= 0 ? 1-gen_p : -9 ;
		else if ( par::model_perm_dom )
			results[l] = dom_p >= 0 ? 1-dom_p : -9 ;
		else if ( par::model_perm_rec )
			results[l] = rec_p >= 0 ? 1-rec_p : -9;
		else if ( par::model_perm_trend )
			results[l] = CA_p >= 0 ? CA_chisq : -9;


		if ( ! par::trend_only )
		{
			double lOR;
			double SE;
			double ci_zt ;
			stringstream ss;
			LampAssoc *outputla = new LampAssoc;
			/////////////
			// Dominant
			//printf("lampassoc_snp:%s\n",locus[l]->name.c_str());
			//fprintf(stderr,"ci_zt:%f\n",ci_zt);
			if (par::model_perm_dom) 
			{
				//fprintf(stderr,"%d/%d\n",l,nl_all);
				outputla->chr=locus[l]->chr;
				outputla->snpname=locus[l]->name;
				outputla->allele1= locus[l]->allele1;
				outputla->allele2= locus[l]->allele2;
				outputla->test="DOM";
				outputla->AFF=int2str(A11+A12)+"/"+int2str(A22);
				outputla->UNAFF=int2str(U11+U12)+"/"+int2str(U22);
				if (A22==0 || U11+U12==0){
					outputla->OR="INF";
					outputla->CIL="NA";
					outputla->CIU="NA";
				}else{
					outputla->OR= dbl2str((((double) A11+(double)A12 ) 
								* (double)U22 )/((double)A22 
								* ( (double)U11+(double)U12 )));
					lOR=log((((double) A11+(double)A12 )
								* (double)U22 )/((double)A22
								* ( (double)U11+(double)U12 )));

					//fprintf(stderr,"%d,%d,%d,%d\n",A11+A12,U22,A22, U11+U12);
					if(A11+A12==0 ||  U22==0 || A22==0 ||  U11+U12==0){
						outputla->CIL="NA";
						outputla->CIU="NA";
					}else{
						SE=sqrt( 1/((double) A11+(double)A12 ) 
								+ 1/(double)U22 + 1/(double)A22 
								+ 1/( (double)U11+(double)U12 ) );
						ci_zt = ltqnorm( 1.0 - (1.0 - par::ci_level) / 2.0  );
						par::ci_zt = ci_zt;
						//fprintf(stderr,"lOR::%f\n",lOR);
						//fprintf(stderr,"SE::%f\n",SE);
						//fprintf(stderr,"par::ci_zt:%f\n",ci_zt);
						outputla->CIL=dbl2str(exp( lOR - par::ci_zt * SE ));
						outputla->CIU=dbl2str(exp( lOR + par::ci_zt * SE ));
					}
				}
				if (dom_p < -1)
				{
					if ( ! par::fisher_test )
					{
						outputla->chsq="NA";
						outputla->DF="NA";
					}
					outputla->P="NA";
				}
				else
				{
					if ( ! par::fisher_test )
					{
						//printf("ss!\n");
						outputla->chsq=dbl2str(dom_chisq);
						outputla->DF="1";
					}
					outputla->P=dbl2str(dom_p);
				}
			}
			/////////////
			// Recessive
			if (par::model_perm_rec)
			{
				outputla->chr=locus[l]->chr;
				outputla->snpname=locus[l]->name;
				outputla->allele1= locus[l]->allele1;
				outputla->allele2= locus[l]->allele2;
				outputla->test="REC";
				outputla->AFF=int2str(A11)+"/"+int2str(A12+A22);
				outputla->UNAFF=int2str(U11)+"/"+int2str(U12+U22);
				if (A12+A22==0 || U11==0){
					outputla->OR="INF";
					outputla->CIL="INF";
					outputla->CIU="INF";
				}else{
					outputla->OR=dbl2str( ((double)A11 * ((double)U12+(double)U22) )/( ((double)A12+(double)A22) * (double)U11 ));
					lOR=log(((double)A11 * ((double)U12+(double)U22) )/( ((double)A12+(double)A22) * (double)U11 ));
					if(A11==0 ||  U12+U22==0 || A12+A22==0 ||  U11==0){
						outputla->CIL="NA";
						outputla->CIU="NA";
					}else{
						SE=sqrt( 1/(double)A11 
								+ 1/((double)U12+(double)U22) 
								+ 1/((double)A12+(double)A22) + 1/(double)U11 );
						ci_zt = ltqnorm( 1 - (1 - par::ci_level) / 2  );
						outputla->CIL=dbl2str(exp( lOR - par::ci_zt * SE ));
						outputla->CIU=dbl2str(exp( lOR + par::ci_zt * SE ));
					}
				}
				if (rec_p < -1)
				{
					if ( ! par::fisher_test )
					{
						outputla->chsq="NA";
						outputla->DF="NA";
					}
					outputla->P="NA";
				}
				else
				{
					if ( ! par::fisher_test )
					{
						outputla->chsq=dbl2str(rec_chisq);
						outputla->DF="1";
					}
					outputla->P=dbl2str(rec_p);
				}
			}
			//printf("lampassoc_snp:%s\n",(outputla->snpname).c_str());
			lampassoc.push_back(outputla);
		}
		s++;
		l++;
	}
}

void Plink::makeLampInput()
{
	printf("Plink::makeLampInput\n");
	vector<Individual*>::iterator gperson = sample.begin();
	ofstream iASC,vASC;
	int l = 0 ;
	string ifname=par::output_file_name + ".item";
	string vfname=par::output_file_name + ".value";
	iASC.open(ifname.c_str(),ios::out);
	vASC.open(vfname.c_str(),ios::out);

	vASC<<"#IID,affection status\n";
	iASC << "#IID,";
	int sep_flag=0;
	vector<CSNP*>::iterator s = SNP.begin();
	int n=0;
	while(s !=SNP.end())
	{
		if (locus[n]->freq < par::MAF_UPPER)
		{
			if (s!=SNP.begin() && sep_flag==1)
			{
				iASC << "," ;
				sep_flag =0;
			}
			iASC << locus[n]->name ;
			sep_flag=1;
		}
		s++;
		n++;
	}
	iASC << "\n" ;
	int snp_count=0;
	while ( gperson != sample.end() )
	{
		// Phenotype for this person (i.e. might be permuted)
		Individual * pperson = (*gperson)->pperson;
		if (!pperson->missing)
		{
			//printf("pperson->fid:%s,pperson->iid:%s\n",pperson->fid.c_str(),pperson->iid.c_str());
			vASC<< pperson->fid << "_" << pperson->iid <<",";
			if (!par::assoc_utest) {
			  if (pperson->aff)
			    vASC<< "1\n";
			  else
			    vASC<< "0\n";
			}
			else {
			  vASC<< pperson->phenotype;
			  vASC<< "\n";
			}
			//iASC << "\n" ;
			iASC << pperson->fid << "_" << pperson->iid << ",";
			vector<CSNP*>::iterator s = SNP.begin();
			n=0;
			sep_flag=0;
			//printf("ID::%s_%s",pperson->fid.c_str(),pperson->iid.c_str());
			snp_count=0;
			while ( s != SNP.end() )
			{	
				vector<bool>::iterator i1 = (*s)->one.begin();
				vector<bool>::iterator i2 = (*s)->two.begin();
				// SNP alleles
				bool s1 = *(i1 + l)  ;
				bool s2 = *(i2 + l)  ;
				//	printf("%d,%d/",s1,s2);
				if ( locus[n]->freq < par::MAF_UPPER )
				{
					//	printf("%d,%d/",s1,s2);
					if (s != SNP.begin() && sep_flag==1 )
					{
						iASC << "," ;
						sep_flag=0;
					}
					if ( ! s1 )
					{
						if ( ! s2 )
						{// Homozyg 00 
							iASC << "1"	;
						}
						else
						{// Hetero  01
							//printf("Hetero\n");
							if(par::model_perm_dom) iASC << "1";
							else if(par::model_perm_rec) iASC << "0";
						}
					}
					else if ( s2 )      // Homozyg 11
						iASC << "0";
					else
						iASC << "0";// missing
					sep_flag=1;
					snp_count++;
				}		
				s++;
				n++;
			}
			//printf("\n");
			//printLOG("lamp_snps:"+int2str(snp_count)+"\n");
			iASC << "\n";
		}
		gperson++;
		l++;
	}
	iASC.close();
	vASC.close();
	printLOG("lamp_snps:"+int2str(snp_count)+"\n");
	if(snp_count==0)
		error("lamp_snps= 0 SNPs! \n");
}
void Plink::DoLamp()
{
	LampCore lampCore;
	std::string set_method;
	// run LAMP
	try {
		std::string log_file;
		time_t now;
		std::time(&now);
		struct tm* d = localtime(&now);
		int year = d->tm_year + 1900;
		int month = d->tm_mon + 1;
		int day = d->tm_mday;
		int hour = d->tm_hour;
		int minute = d->tm_min;
		int second = d->tm_sec;
		log_file = "lamp_log_";// = "lamp_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
		log_file += std::to_string(year) + std::to_string(month) + std::to_string(day) + "_";
		log_file += std::to_string(hour) + std::to_string(minute) + std::to_string(second) + ".txt";
		if (par::fisher_test) {
			set_method = "fisher";
		} else if (par::assoc_utest) {
			set_method = "u_test";
		} else {
			set_method = "chi";
		}
		std::string ifile = par::output_file_name + ".item";
		std::string vfile = par::output_file_name + ".value";
		lampCore.run(ifile, vfile, par::SGLEV, set_method, -1, log_file, par::lamp_alternative);
	} catch (...) {
		error("Lamp error. Please check inputfile.\n");
	}
	
	// convert results
	try {
		int flag_size = lampCore.getN1();
		int tran_size = lampCore.getTransactionSize();
		int k = lampCore.getK();
		int lam_star = lampCore.getLam_star();
		std::vector<LampCore::enrich_t*> enrich_lst = lampCore.getEnrich_lst();
		const std::vector<std::string*> columnid2name = lampCore.getColumnid2name();
		printLOG("-----lamp::info--------\n");
		printLOG("# LAMP ver. " + LampCore::getVersion() + "\n");
		printLOG("# item-file: " + par::output_file_name + ".item" + "\n");
		printLOG("# value-file: " + par::output_file_name + ".value" + "\n");
		printLOG("# significance-level: " + to_string(par::SGLEV) + "\n"); 
		printLOG("# P-value computing procedure: " + set_method);
		if (0 < par::lamp_alternative) {
			printLOG(" (greater)");
		}
		else if (par::lamp_alternative < 0) {
			printLOG(" (less)");
		}
		else {
			printLOG(" (two.sided)");
		}
		printLOG("\n");
		printLOG("# # of tested elements: " + to_string(columnid2name.size()) + ", # of samples: " + to_string(tran_size) );
		if (0 < flag_size)
			printLOG(", # of positive samples: " + flag_size);
		printLOG("\n");
		printLOG("# Adjusted significance level: " + (boost::format("%.5g") % (par::SGLEV/k)).str() + ", "); 
		printLOG("Correction factor: " + to_string(k) + " (# of target rows >= " + to_string(lam_star) + ")\n");
		printLOG("# # of significant combinations: " + to_string(enrich_lst.size()) + "\n");

		// if (0 < enrich_lst.size() && 0 < flag_size) {
		if (0 < enrich_lst.size()) {
			std::sort(enrich_lst.begin(), enrich_lst.end(), LampCore::cmpEnrich);
			int rank = 0;
			for (LampCore::enrich_t* l : enrich_lst) {
				rank = rank + 1;
				Lamp * res_lamp = new Lamp;
				int x = flag_size - l->stat_score;
				int y = l->p - l->len;
				int z = lampCore.getTransactionSize() - (l->stat_score + x + y);
				res_lamp->n11 = l->stat_score;
				res_lamp->n10 = x;
				res_lamp->n01 = y;
				res_lamp->n00 = z;
				res_lamp->Rank = to_string(rank);
				res_lamp->Raw_P_value = (boost::format("%.5g") % (l->p)).str();
				res_lamp->Ajusted_P_value = (boost::format("%.5g") % (k * l->p)).str();
				std::string out_column = "";
				for (int i : *l->item_set) {
					out_column += *columnid2name[i-1] + ",";
				}
				out_column.erase( --out_column.end() );
				res_lamp->COMB = out_column;
				lamp.push_back(res_lamp);
			}
		} else {
		  printLOG("No significant SNP combinations were detected.\n");
		}
	} catch (...) {
		error("Lamp error. Please check setting.\n");
	}
}

void Plink::setLampformat()
{
	lamplink_format["CHR"]=4;
	lamplink_format["SNP"]=par::pp_maxsnp;
	lamplink_format["A1"]=4;
	lamplink_format["A2"]=4;
	lamplink_format["TEST"]=8;
	lamplink_format["AFF"]=14;
	lamplink_format["UNAFF"]=14;
	lamplink_format["CHISQ"]=12;
	lamplink_format["DF"]=4;
	lamplink_format["P"]=12;
	lamplink_format["OR"]=12;
	lamplink_format["LXX"]=12;
	lamplink_format["UXX"]=12;
	lamplink_format["COMB"]=8;

	lamp_format["COMBID"]=10;
	lamp_format["Raw_P"]=12;
	lamp_format["Adjusted_P"]=12;
	lamp_format["OR"]=12;
	lamp_format["LXX"]=12;
	lamp_format["UXX"]=12;
}


void Plink::outputLamplink()
{
	ofstream ASC,lASC;
	string outputfile=par::output_file_name + ".lamplink";
	string loutputfile=par::output_file_name + ".lamp";
	lASC.open(loutputfile.c_str(),ios::out);
	ASC.open(outputfile.c_str(),ios::out);
	printf("Plink::outputLamplink\n");
	//make_header
	ASC << setw(4) <<"CHR" <<" "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(4) << "A1" << " "
		<< setw(4) << "A2" << " "
		<< setw(8) << "TEST" << " "
		<< setw(14) << "AFF" << " "
		<< setw(14) << "UNAFF" << " " ;
	if ( ! par::fisher_test )
		ASC << setw(12) << "CHISQ" << " "
			<< setw(4) << "DF" << " ";
	ASC << setw(12) << "P" << " " 
		<<setw(12)<<"OR"<<" ";
	if (par::display_ci)
		ASC <<setw(12) << string("L"+dbl2str(par::ci_level*100)) <<" "
			<<setw(12)<<string("U"+dbl2str(par::ci_level*100)) <<" ";

	//make_header lampfile
	lASC << setw(10) <<"COMBID" <<" " 
		<<setw(12)<< "Raw_P" <<" "
		<<setw(12)<<"Adjusted_P" <<" ";
	if (par::display_ci)
		lASC	<<setw(12)<<"OR"
			<<setw(12) << string("L"+dbl2str(par::ci_level*100)) <<" "
			<<setw(12)<<string("U"+dbl2str(par::ci_level*100)) <<" ";
	lASC << "COMB" <<"\n"; 
	map<string,string> SNP_COMB ;
	vector<Lamp *>::iterator glamp=lamp.begin();
	int n=1;
	//make 
	while(glamp != lamp.end())
	{
		SNP_COMB[string("COMB"+int2str(n))]=(*glamp)->COMB;
		ASC<<setw(8)<<string("COMB"+int2str(n))<<" ";
		//make lamp file
		lASC << setw(10) <<"COMB"+int2str(n) <<" "
			<<setw(12)<< (*glamp)->Raw_P_value <<" "
			<<setw(12)<< (*glamp)->Ajusted_P_value<<" ";
		if (par::display_ci)
		{
			double lOR;
			double SE;
			double ci_zt ;
			fprintf(stderr,"%d,%d,%d,%d",(*glamp)->n11,(*glamp)->n00,(*glamp)->n10,(*glamp)->n01);
			if((*glamp)->n10==0 || (*glamp)->n01==0)
			{
				lASC    <<setw(12)<<"INF" <<" "
					<<setw(12) <<"NA" <<" "
					<< setw(12)<<"NA"<<" ";
			}else{
				lASC <<setw(12)<< dbl2str( ((double)(*glamp)->n11 * (double)(*glamp)->n00 )
						/((double)(*glamp)->n10 * (double)(*glamp)->n01)) <<" ";
				if((*glamp)->n11==0 || (*glamp)->n00==0){
					lASC <<setw(12) <<"NA" <<" "<< setw(12)<<"NA"<<" ";
				}else{
					SE=sqrt(1/(double)(*glamp)->n11 + 1/(double)(*glamp)->n10 + 1/(double)(*glamp)->n01 + 1/(double)(*glamp)->n00 );
					ci_zt = ltqnorm( 1 - (1 - par::ci_level) / 2  );
					lOR=log(((double)(*glamp)->n11 * (double)(*glamp)->n00)/((double)(*glamp)->n10 * (double)(*glamp)->n01));
					lASC <<setw(12) <<dbl2str(exp( lOR - par::ci_zt * SE )) <<" "
						<< setw(12)<<dbl2str(exp( lOR + par::ci_zt * SE )) <<" ";
				}
			}
		}
		lASC << (*glamp)->COMB <<"\n";
		n++;
		glamp++;
	}
	lASC.close();
	ASC<<"\n";
	int snpcounter=0;
	vector<LampAssoc *>::iterator glampassoc=lampassoc.begin();
	//make_data
	while(glampassoc != lampassoc.end())
	{
		ASC << setw(4) << int2str((*glampassoc)->chr) << " "
			<<setw(par::pp_maxsnp)<< (*glampassoc)->snpname <<" "
			<<setw(4)<< (*glampassoc)->allele1 <<" "
			<<setw(4)<< (*glampassoc)->allele2 <<" "
			<<setw(8)<< (*glampassoc)->test <<" "
			<<setw(14)<<(*glampassoc)->AFF <<" "
			<<setw(14)<<(*glampassoc)->UNAFF<<" ";
		if(! par::fisher_test)
			ASC<< setw(12) << (*glampassoc)->chsq 
				<<setw(4)<<(*glampassoc)->DF <<" ";
		ASC << setw(12) << (*glampassoc)->P << " "
			<<setw(12)<< (*glampassoc)->OR <<" ";
		if (par::display_ci)
			ASC <<setw(12) << (*glampassoc)->CIL <<" "
				<<setw(12)<< (*glampassoc)->CIU <<" ";
		n=1;

		while( SNP_COMB.find(string("COMB"+int2str(n))) != SNP_COMB.end())
		{
			string target_comb = SNP_COMB[string("COMB"+int2str(n))];
			size_t current=0,found;
			int res=0;
			vector<string> comb_snps =split(target_comb,',');
			vector<string>::iterator icomb_snp=comb_snps.begin();
			while(icomb_snp!=comb_snps.end())
			{
				if((*icomb_snp).compare((*glampassoc)->snpname)==0)res=1;
				icomb_snp++;
			}
			if(res==1) ASC<<setw(8)<<"1"<<" ";
			else ASC<<setw(8)<<"0"<<" ";
			n++;
		}
		ASC<<"\n";
		glampassoc++;
		snpcounter++;
	}
	ASC.close();
}

void Plink::readLamplinkfile()
{
	//FILE *fLAMP;
	string fname= par::combfilename+".lamplink";
	checkFileExists(fname);
	printLOG("read_lamplinkfile.........."+fname+"\n");
	ifstream fLAMP;
	fLAMP.open(fname.c_str());
	fLAMP.clear();
	int l_counter=0;
	int sep_count;
	while(!fLAMP.eof())
	{
		readLamplink *inputlp = new readLamplink;
		string sline;
		int col_count=0;
		getline(fLAMP,sline);
		if (sline=="") continue;
		string buf="";
		stringstream ss(sline);
		while(ss>>buf)
		{
			if(l_counter ==0){
				if(buf.compare("COMB1")==0) sep_count=col_count;
				//printf("header:%s\n",buf.c_str());
				inputlp->header.push_back(buf);
			}else{
				//printf("col:%d,,buf:%s\n",col_count,buf.c_str());
				if(col_count>sep_count-1)
				{
					//printf("combination_flag:%s\n",buf.c_str());
					inputlp->combination.push_back(buf);
				}
				else inputlp->item.push_back(buf);
			}
			col_count++;
		}
		l_counter=1;
		readlamplink.push_back(inputlp);
	}
	printf("sep_count:::%d\n",sep_count);
}
void Plink::readLampfile(){
	//FILE *rLAMP;
	string fname=par::combfilename+".lamp";
	checkFileExists(fname);
	printLOG("read_lampfile.........."+fname+"\n");
	ifstream fLAMP;
	fLAMP.open(fname.c_str());
	fLAMP.clear();
	int comb_counter=0;
	int sep_count;
	while(!fLAMP.eof())
	{
		readLamp *inputl =new readLamp;
		string sline;
		int col_count=0;
		getline(fLAMP,sline);
		if (sline=="") continue;
		string buf="";
		stringstream ss(sline);
		while(ss>>buf)
		{
			if(comb_counter ==0)
			{
				if(buf.compare("COMB")==0 && buf.compare("COMBID")!=0) sep_count=col_count;
				inputl->header.push_back(buf);
			}else{
				if(col_count==sep_count){
					inputl->combination = buf;
				}
				else inputl->item.push_back(buf);
			}
			col_count++;
		}
		comb_counter++;
		readlamp.push_back(inputl);
	}
	combcount=comb_counter-1;
	printf("combination:%d\n",combcount);
	//return comb_counter;
}

void Plink::LampLDcalc()
{
	ofstream rmlamp,rmlamplink;
	string f = par::output_file_name + ".lamp";
	rmlamp.open(f.c_str(),ios::out);
	set<string>remove_combid;
	set<int> ldAnchorSet;
	map<string,int> mlocus;
	//fprintf(stderr,"start_LampLDCalc\n");
	for (int l=0; l<nl_all; l++)
		mlocus.insert( make_pair( locus[l]->name , l ) );
	vector<readLamp*>::iterator ir = readlamp.begin();
	vector<string>lheader=(*ir)->header;
	//make_header4lampfile
	vector<string>::iterator ih =lheader.begin();

	while(ih!=lheader.end())
	{
		if((*ih).compare("COMB")!=0)
		{
			if (string((*ih),0,1).compare("U")==0) rmlamp<<setw(lamp_format["UXX"])<<(*ih)<<" ";
			else if (string((*ih),0,1).compare("L")==0) rmlamp<<setw(lamp_format["LXX"])<<(*ih)<<" ";
			else rmlamp<<setw(lamp_format[*ih])<<(*ih)<<" ";
		}
		else
			rmlamp<<(*ih);
		ih++;
	}
	rmlamp<<"\n";
	set<double> r_set;
	while( ir != readlamp.end() )
	{
		if(ir != readlamp.begin())
		{
			vector<string> snps =split((*ir)->combination,',') ;
			vector<string>::iterator snp = snps.begin();
			if(snps.size()>1)
			{
				ldAnchorSet.clear();
				while(snp!=snps.end())
				{
					map<string,int>::iterator i = mlocus.find( *snp );

					if ( i != mlocus.end() )
						ldAnchorSet.insert( i->second );
					snp++;
				}
				int end = nl_all;
				r_set.clear();
				for (int l1=0; l1<end; l1++)
				{
					int start=0;
					if ( ldAnchorSet.find ( l1 ) == ldAnchorSet.end() )
						continue;
					for (int l2=start; l2<nl_all; l2++)
					{
						if ( locus[l2]->chr != locus[l1]->chr )
							continue;
						if ( ldAnchorSet.find ( l2 ) == ldAnchorSet.end() )
							continue;
						if (l1 < l2)
						{
							double r = correlation2SNP(l1,l2, false, false);
							r *= r;
							fprintf(stdout,"l1:%s\tl2:%s::%f\n",(locus[l1]->name).c_str(),(locus[l2]->name).c_str(),r);
							r_set.insert(r);
						}
					}
				}
				if(r_set.lower_bound(par::LAMP_R2)!=r_set.end())
				{
					vector<string>::iterator rmid=(*ir)->item.begin();
					printLOG("remove_combid:"+(*rmid)+"\n");
					remove_combid.insert((*rmid));
					ir++;
					continue;
				}
			}
			//for(vector<string>::iterator it = snps.begin(); it != snps.end(); it++ )
			//fprintf(stderr,"%s",((*ir)->combination).c_str());
			//fprintf(stderr,"::OK\n");
			int l=0;
			vector<string>::iterator ih =lheader.begin();
			while(ih!=lheader.end())		 
			{
				vector<string>::iterator item=(*ir)->item.begin();
				if(item+l != (*ir)->item.end()){
					if (string((*ih),0,1).compare("U")==0) rmlamp<<setw(lamp_format["UXX"])<<*(item+l)<<" ";
					else if (string((*ih),0,1).compare("L")==0) rmlamp<<setw(lamp_format["LXX"])<<*(item+l)<<" ";
					else rmlamp<<setw(lamp_format[*ih])<<*(item+l)<<" ";
				}else{
					rmlamp<<(*ir)->combination<<" ";
				}
				ih++;
				l++;
			}
			rmlamp<<"\n";
		}
		ir++;
	}
	rmlamp.close();
	f=par::output_file_name + ".lamplink";
	rmlamplink.open(f.c_str(),ios::out);
	vector<readLamplink*>::iterator irp = readlamplink.begin();
	vector<string>lpheader=(*irp)->header;
	vector<string>::iterator ihp =lpheader.begin();
	while(ihp!=lpheader.end())
	{
		if(string((*ihp),0,4).compare("COMB")!=0)
		{
			if(string((*ihp)).compare("UNAFF")!=0 &&
					(string((*ihp),0,1).compare("U")==0 ||string((*ihp),0,1).compare("L")==0))
			{
				if (string((*ihp),0,1).compare("U")==0) rmlamplink<<setw(lamplink_format["UXX"])<<(*ihp)<<" ";
				else if (string((*ihp),0,1).compare("L")==0) rmlamplink<<setw(lamplink_format["LXX"])<<(*ihp)<<" ";

			}
			else
				rmlamplink<<setw(lamplink_format[*ihp])<<(*ihp)<<" ";
		}
		else
		{
			if(remove_combid.find(*ihp)==remove_combid.end())
				rmlamplink<<setw(lamplink_format["COMB"])<<(*ihp)<<" ";
		}
		ihp++;
	}
	rmlamplink<<"\n";
	//printf("header_size:%d\n",(*irp)->header.size());
	while( irp != readlamplink.end() )
	{
		vector<string>::iterator ihp =lpheader.begin();
		int l=0,m=0;
		if(irp!=readlamplink.begin())
		{
			vector<string>::iterator item=(*irp)->item.begin();
			vector<string>::iterator comb=(*irp)->combination.begin();
			while(ihp !=lpheader.end())
			{
				//printf("l::,%d////m;;,%d\n",l,m);
				//vector<string>::iterator item=(*irp)->item.begin();
				//vector<string>::iterator comb=(*irp)->combination.begin();
				//if(remove_combid.find(*ihp)!=remove_combid.end())
				//{
				//printf("header%d::%s\n",l,(*ihp).c_str());
				if(item+l != (*irp)->item.end())
				{
					if(string((*ihp)).compare("UNAFF")!=0 &&
							(string((*ihp),0,1).compare("U")==0 ||string((*ihp),0,1).compare("L")==0))
					{
						if (string((*ihp),0,1).compare("U")==0) rmlamplink<<setw(lamplink_format["UXX"])<<*(item+l)<<" ";
						else if (string((*ihp),0,1).compare("L")==0) rmlamplink<<setw(lamplink_format["LXX"])<<*(item+l)<<" ";

					}
					else
						rmlamplink<<setw(lamplink_format[*ihp])<<*(item+l)<<" ";
					l++;
				}
				else
				{
					vector<string>::iterator comb=(*irp)->combination.begin();
					//	printf("comb_size:%d - %d\n",(*irp)->combination.size(),m);
					//	printf("combID%d::%s\n",m,(*ihp).c_str());	
					//i++;
					//	printf("lamplink_combflag::%s\n",(*(comb+m)).c_str());
					if( (comb+m) !=(*irp)->combination.end() )
					{
						if(remove_combid.find(*ihp)==remove_combid.end()){
							rmlamplink << setw(lamplink_format["COMB"])<<(*(comb + m)) <<" ";
						}
						m++;
					}
				}
				ihp++;
			}
			rmlamplink<<"\n";
			}
			irp++;
		}
		rmlamplink.close();	



	}

	void Plink::calcAssociationWithPermutation4Lamp(Perm & perm)
	{
		// SNP-major mode analyses?

		if (par::assoc_glm)
		{
			if (par::SNP_major) 
				SNP2Ind();
		}
		else if (!par::SNP_major) 
			Ind2SNP();


		//////////////////////////////
		// Profile-based set test?

		//   if ( par::set_score )
		//     pS->profileTestInitialise();


		//////////////////////////////
		// LD-clump within each set?

		if ( par::set_test && par::set_r2 )
		{
			printLOG("Performing LD-based set test, with parameters:\n");
			printLOG("     r-squared  (--set-r2)   = " + dbl2str( par::set_r2_val ) + "\n" );
			printLOG("     p-value    (--set-p)    = " + dbl2str( chiprobP(par::set_chisq_threshold,1) ) + "\n" );
			printLOG("     max # SNPs (--set-max)  = " + int2str( par::set_max ) + "\n" );
			pS->makeLDSets();
		}

		//////////////////////////////
		// Step-wise set tests

		if ( par::set_step )
		{
			vector_t r = pS->fitStepwiseModel();
			shutdown();
		}

		// Basic association testing results
		vector<double> original;

		//////////////////////
		// Empirical p-values

		////////////////////
		// Number of tests

		int ntests = nl_all;

		if ( par::set_test && ( par::set_r2 || par::set_score ) ) 
			ntests = pS->snpset.size();

		perm.setTests(ntests);


		// Case/control missingness test statistics
		vector<double> missing;

		// Observed marginals
		int aff;
		int unf;

		// Odds ratio
		vector<double> odds(nl_all);

		//////////////////////////////////////////
		// Working vectors for assoc_test_alt_perm

		vector<int> a1;
		vector<int> a2;
		vector<int> a0;

		// Expected values for the 2x2 test
		vector<double> exp_afffreq1;
		vector<double> exp_afffreq2;
		vector<double> exp_unffreq1;
		vector<double> exp_unffreq2;

		if (par::assoc_test_alt_perm)
		{
			a1.resize(nl_all);
			a2.resize(nl_all);
			a0.resize(nl_all);
			exp_afffreq1.resize(nl_all);
			exp_afffreq2.resize(nl_all);
			exp_unffreq1.resize(nl_all);
			exp_unffreq2.resize(nl_all);
		}


		////////////////////////////////
		// Set up permutation structure 
		// (we need to perform this step
		//  whether or not we also 
		//  subsequently permute)

		perm.setPermClusters(*this);
		perm.originalOrder();

		////////////////////////////////////////////////////////
		// Perform a test of missingness by case/control status

		if (par::test_missing && par::bt)
			original = testMiss(perm,true);


		// If we didn't know how many values to expect back, 
		// resize now (i.e. from haplotype tests)

		if ( par::test_hap_GLM )
		{
			ntests = original.size();      
			perm.setTests(ntests);
		}


		////////////////////////////
		// If no permutation, then 
		// leave now

		if (!par::permute) 
			return;


		////////////////////////
		// Score original sets

		vector<int> setsigsize;

		if (par::set_test) 
		{
			if ( par::set_r2 || par::set_score )
			{
				// Score...

				// original = pS->profileTestScore();

				original = pS->fitLDSetTest(original,true);

				// ...and save # of significant SNPs	  
				setsigsize.clear();
				for (int i=0; i<pS->profileSNPs.size(); i++)
					setsigsize.push_back( pS->s_min[i] );
			}
			else
				pS->cumulativeSetSum_WITHLABELS(*this,original);
		}


		/////////////////////////////
		// Ordered/rank permutation?

		if (par::mperm_rank)
			perm.setOriginalRanking(original);


		//////////////////////
		// Verbose dumping?

		if (par::mperm_save_all)
			printLOG("Dumping all permutation statistics to [ "
					+ par::output_file_name+".mperm.dump.all ]\n");
		else if (par::mperm_save_best)
			printLOG("Dumping best permutation statistics to [ "
					+ par::output_file_name+".mperm.dump.best ]\n");



		//////////////////////
		// Begin permutations

		bool finished = par::replicates == 0 ? true : false;

		while(!finished)
		{

			// Store permuted results

			vector<double> pr(ntests);


			if (par::perm_genedrop)
			{
				if (par::perm_genedrop_and_swap)
					perm.permuteInCluster();
				perm.geneDrop();
			}
			else
				perm.permuteInCluster();

			if (par::test_missing)
				pr = testMiss(perm,false);
			else if ((!par::assoc_test_alt_perm) 
					|| par::qt 
					|| par::full_model_assoc 
					|| par::CMH_test_1
					|| par::assoc_glm)
			{

				if (par::qt)
				{
					if (par::assoc_glm)
						pr = glmAssoc(false,perm);
					else if ( par::test_hap_GLM )
						pr = haplo->phaseAllHaplotypes(false,perm);
					else
						pr = testQAssoc(false,perm);

				}
				else if (par::full_model_assoc)
					pr = fullModelAssoc(false,perm);
				else if (par::assoc_glm)
					pr = glmAssoc(false,perm);
				else if (par::CMH_test_1)
					pr = calcMantelHaenszel_2x2xK(perm, false);
				else if ( par::test_hap_GLM )
					pr = haplo->phaseAllHaplotypes(false,perm);
				else 
					pr = testAssoc(aff,unf,
							a1,a2,a0,
							odds,	 
							exp_afffreq1, exp_afffreq2,
							exp_unffreq1, exp_unffreq2,
							perm,
							false);
			}
			else
			{

				/////////////////////////
				// For binary traits only

				//  -------------
				//  | A | B | E |  aff
				//  ------------- 
				//  | C | D | F |  unf
				//  -------------
				//   a1  a2  a0

				// a1 most likely to be common, followed by a2, then a0
				// save aff+unf (do not alter by locus)
				// and a1,a2,a0 marginals (which do alter by locus) 
				// then we only need count A and B in each subsequent replicate:
				// 

				int A, B, C, D, M;

				///////////////////////////////
				// Iterate over SNPs

				vector<CSNP*>::iterator s = SNP.begin();
				int l=0;

				while ( s != SNP.end() )
				{	

					// In adaptive mode, possibly skip this test
					if (par::adaptive_perm && (!perm.snp_test[l])) 
					{
						s++;
						l++;
						continue;
					}


					///////////////// 	      
					// clear counts

					D=M=0;

					///////////////// 	      
					// Autosomal or haploid?

					bool X=false, haploid=false;
					if (par::chr_sex[locus[l]->chr]) X=true;
					else if (par::chr_haploid[locus[l]->chr]) haploid=true;

					/////////////////////////////
					// Iterate over individuals

					vector<bool>::iterator i1 = (*s)->one.begin();
					vector<bool>::iterator i2 = (*s)->two.begin();
					vector<Individual*>::iterator gperson = sample.begin();

					while ( gperson != sample.end() )
					{

						// Phenotype for this person (i.e. might be permuted)
						Individual * pperson = (*gperson)->pperson;

						// SNP alleles

						bool s1 = *i1;
						bool s2 = *i2;

						if ( ! pperson->missing )
						{
							if (! pperson->aff ) // unaffected 
							{	      
								if ( haploid || ( X && (*gperson)->sex ) ) 
								{ 
									if ( s2 )
										D++;               // (hemi, one count)
									else if ( s1 ) M++;  // (missing, one count)
								}
								else
								{
									if ( s2 )
									{
										if (!s1) D++;    // (het, one A count)
										else D+=2;       // (hom, two B count)
									}
									else if ( s1 ) M+=2; // (missing, two B count)
								}
							}
						}

						// Next individual
						gperson++;
						i1++;
						i2++;

					}


					// reconstruct rest of 2x2 table
					C = unf - D - M;
					A = a1[l] - C;
					B = a2[l] - D;

					pr[l] = ( (A - exp_afffreq1[l]) * ( A - exp_afffreq1[l] ) ) / exp_afffreq1[l] 
						+ ( (C - exp_unffreq1[l]) * ( C - exp_unffreq1[l] ) ) / exp_unffreq1[l] 
						+ ( (B - exp_afffreq2[l]) * ( B - exp_afffreq2[l] ) ) / exp_afffreq2[l] 
						+ ( (D - exp_unffreq2[l]) * ( D - exp_unffreq2[l] ) ) / exp_unffreq2[l];


					// Next SNP
					s++;
					l++;

				}
			}
		}
	}
