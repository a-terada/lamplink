
//////////////////////////////////////////////////////////////////
//                                                              //
//          LAMPLINK (c) 2015-2017 LAMP development team        //
//                                                              //
// This code is written to add options                          //
// (--lamp and --lamp-ld-removed) to PLINK v1.07                //
// (http://zzz.bwh.harvard.edu/plink/)                          //
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
#include <memory>
#include <cassert>
#include <boost/format.hpp>

#include "plink.h"
#include "options.h"
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
#include "FastWYCore.h"

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
	printLOG("Full-model association for Lamp tests\n");

	lampassoc.clear();

	int	transaction_size = 0;
	int	n1_count = 0;
	for (const auto gperson : sample) {
		const Individual*	pperson = gperson->pperson;
		if (!pperson->missing) {
			if (pperson->aff) {
				++n1_count;
			}
			++transaction_size;
		}
	}
	if (par::lamp_alternative < 0) {	//less
		n1_count = transaction_size - n1_count;
	}

	std::unique_ptr<FunctionsSuper>	tFunc;
	if (par::assoc_utest) {
		tFunc.reset(new Functions4u_test(transaction_size, n1_count, par::lamp_alternative));
	} else if (par::fisher_test) {
		tFunc.reset(new Functions4fisher(transaction_size, n1_count, std::abs(par::lamp_alternative)));
	} else {
		tFunc.reset(new Functions4chi(transaction_size, n1_count, std::abs(par::lamp_alternative)));
	}

	auto newLampAccoc = [](const vector<Locus*>& locus, int l) {
		LampAssoc *outputla = new LampAssoc;
		outputla->chr = locus[l]->chr;
		outputla->snpname = locus[l]->name;
		outputla->allele1 = locus[l]->allele1;
		outputla->allele2 = locus[l]->allele2;
		if (par::model_perm_dom) {
			outputla->test="DOM";
		} else {
			outputla->test="REC";
		}
		return outputla;
	};

	bool	reverseValueB = par::lamp_alternative < 0;
	double	pvalue, score, stat_value;
	int		l = 0;
	for (const auto s : SNP) {
		auto	i1 = s->one.begin();
		auto	i2 = s->two.begin();
		
		LampAssoc *outputla = newLampAccoc(locus, l);
		if (par::assoc_utest) {
			std::vector<double>	group_x;
			std::vector<double>	group_y;
			for (const auto gperson : sample) {
				bool	s1 = *i1++;
				bool	s2 = *i2++;
				const Individual*	pperson = gperson->pperson;
				if (!pperson->missing) {
					bool	s;
					if (!s1) {
						if (!s2) {		// Homozyg 00
							s = true;
						} else {		// Hetero  01
							if (par::model_perm_dom) {
								s = true;
							} else { // par::model_perm_rec
								s = false;
							}
						}
					} else if (s2) {	// Homozyg 11
						s = false;
					} else {			// missing
						s = false;
					}
					if (s) {
						group_x.push_back(pperson->phenotype);
					} else {
						group_y.push_back(pperson->phenotype);
					}
				}
			}
			std::sort(group_x.begin(), group_x.end());
			std::sort(group_y.begin(), group_y.end());
			pvalue = static_cast<Functions4u_test*>(tFunc.get())->calPValue(group_x, group_y, score, stat_value);
			outputla->U = dbl2str(stat_value);
			outputla->P = dbl2str(pvalue);
			lampassoc.push_back(outputla);
		} else {
			auto	calcOddsFunc = [](double ovalues[2][2], bool reverse, double& ciL, double& ciU) {
						double	odds;
						ciL = std::numeric_limits<double>::quiet_NaN();
						ciU = std::numeric_limits<double>::quiet_NaN();
						if (reverse) {
							odds = ovalues[0][1] / ovalues[0][0] *
									ovalues[1][0] / ovalues[1][1];
						} else {
							odds = ovalues[0][0] / ovalues[0][1] *
									ovalues[1][1] / ovalues[1][0];
						}
						if (ovalues[0][0] != 0 && ovalues[1][1] != 0 &&
							ovalues[0][1] != 0 && ovalues[1][0] != 0) {
							double	lOR = std::log(odds);
							double	SE = std::sqrt(1 / ovalues[0][0] +
													1 / ovalues[0][1] +
													1 / ovalues[1][0] +
													1 / ovalues[1][1]);
							ciL = std::exp(lOR - par::ci_zt * SE);
							ciU = std::exp(lOR + par::ci_zt * SE);
						}
						return odds;
					};

			double	ovalues[2][2] = {{0, 0}, {0, 0}};
			double	ovalues_r[2][2] = {{0, 0}, {0, 0}};
			for (const auto gperson : sample) {
				bool	s1 = *i1++;
				bool	s2 = *i2++;
				const Individual*	pperson = gperson->pperson;
				if (!pperson->missing) {
					int	s, s_r;
					if (!s1) {
						if (!s2) {		// Homozyg 00
							s = 0;
							s_r = 1;
						} else {		// Hetero  01
							if (par::model_perm_dom) {
								s = s_r = 0;
							} else { // par::model_perm_rec
								s = s_r = 1;
							}
						}
					} else if (s2) {	// Homozyg 11
						s = 1;
						s_r = 0;
					} else {			// missing
						s = s_r = 1;
					}
					if (pperson->aff == !reverseValueB) {
						ovalues  [s  ][0] += 1;
						ovalues_r[s_r][0] += 1;
					} else {
						ovalues  [s  ][1] += 1;
						ovalues_r[s_r][1] += 1;
					}
				}
			}

			double	ciL, ciU;
			double	odds, odds_r;

			if (par::lamp_allele_test) {
				odds = calcOddsFunc(ovalues, reverseValueB, ciL, ciU);
				odds_r = calcOddsFunc(ovalues_r, reverseValueB, ciL, ciU);
				if (!std::isnan(odds_r) &&
					(std::isnan(odds) || odds_r > odds)) {
					ovalues[0][0] = ovalues_r[0][0];
					ovalues[0][1] = ovalues_r[0][1];
					ovalues[1][0] = ovalues_r[1][0];
					ovalues[1][1] = ovalues_r[1][1];
					auto	i1 = s->one.begin();
					auto	i2 = s->two.begin();
					for (; i1 != s->one.end(); ++i1, ++i2) {
						if (*i1 == *i2) {
							*i1 = !(*i1);
							*i2 = !(*i2);
						}
					}
					locus[l]->allele1 = outputla->allele2;
					locus[l]->allele2 = outputla->allele1;
					outputla->allele1 = locus[l]->allele1;
					outputla->allele2 = locus[l]->allele2;
					outputla->test = outputla->test + "!";
				}
			}

			pvalue = tFunc->calPValue(ovalues, score, stat_value);
			double	chi = 0;
			if (!par::fisher_test) {
				chi = stat_value;
			}

			odds = calcOddsFunc(ovalues, reverseValueB, ciL, ciU);

			outputla->OR = std::isnan(odds) ? "NA" : std::isinf(odds) ? "INF" : dbl2str(odds);
			outputla->CIL = std::isnan(ciL) ? "NA" : std::isinf(ciL) ? "INF" : dbl2str(ciL);
			outputla->CIU = std::isnan(ciU) ? "NA" : std::isinf(ciU) ? "INF" : dbl2str(ciU);

			if (reverseValueB) {
				outputla->AFF=int2str(ovalues[0][1])+"/"+int2str(ovalues[1][1]);
				outputla->UNAFF=int2str(ovalues[0][0])+"/"+int2str(ovalues[1][0]);
			} else {
				outputla->AFF=int2str(ovalues[0][0])+"/"+int2str(ovalues[1][0]);
				outputla->UNAFF=int2str(ovalues[0][1])+"/"+int2str(ovalues[1][1]);
			}
			outputla->chsq = dbl2str(chi);
			outputla->DF = "1";
			outputla->P = dbl2str(pvalue);
			lampassoc.push_back(outputla);
		}
		++l;
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
	LampCore	lampCore;
	FastWYCore	fastWYCore;
	std::string set_method;

	lamp.clear();

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
		if (!par::fastwy) {
			log_file = "lamp_log_";// = "lamp_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
		} else {
			log_file = "fastwy_log_";
		}
		log_file += std::to_string(year) + std::to_string(month) + std::to_string(day) + "_";
		log_file += std::to_string(hour) + std::to_string(minute) + std::to_string(second) + ".txt";
		log_file = par::output_file_name + "." + log_file;
		if (par::fisher_test) {
			set_method = "fisher";
		} else if (par::assoc_utest) {
			set_method = "u_test";
		} else {
			set_method = "chi";
		}
		std::string ifile = par::output_file_name + ".item";
		std::string vfile = par::output_file_name + ".value";
		if (!par::fastwy) {
			lampCore.run(ifile, vfile, par::SGLEV, set_method, -1, log_file, par::lamp_alternative);
		} else {
			fastWYCore.run(ifile, vfile, par::SGLEV, par::fastwy_nperm, set_method, -1, log_file, par::lamp_alternative, false);
		}
	} catch (std::string &msg) {
		error("Lamp error. Please check inputfile. " + msg + "\n");
	} catch (...) {
		error("Lamp error. Please check inputfile.\n");
	}
	
	// convert results
	try {
		if (!par::fastwy) {
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
				printLOG(", # of positive samples: " + to_string(flag_size));
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
					int y = l->len - l->stat_score;
					int z = lampCore.getTransactionSize() - (l->stat_score + x + y);
					if (par::lamp_alternative < 0) {
						z = flag_size - y;
						x = lampCore.getTransactionSize() - (l->stat_score + y + z);
					}
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
		} else {
			const auto&	result_list = fastWYCore.getResultList();
			int flag_size = fastWYCore.getN1();
			printLOG("-----fastWY::info--------\n");
			printLOG(fastWYCore.getResultInfo());
			if (result_list.size() > 0) {
				int rank = 0;
				for (const auto& l : result_list) {
					rank = rank + 1;
					Lamp * res_lamp = new Lamp;
					int x = flag_size - l.stat_score;
					int y = l.target_rows - l.stat_score;
					int z = fastWYCore.getTransactionSize() - (l.stat_score + x + y);
					if (par::lamp_alternative < 0) {
						z = flag_size - y;
						x = fastWYCore.getTransactionSize() - (l.stat_score + y + z);
					}
					res_lamp->n11 = l.stat_score;
					res_lamp->n10 = x;
					res_lamp->n01 = y;
					res_lamp->n00 = z;
					res_lamp->Rank = to_string(rank);
					res_lamp->Raw_P_value = (boost::format("%.5g") % (l.p)).str();
					res_lamp->Ajusted_P_value = (boost::format("%.5g") % (l.adjusted_p)).str();
					std::string out_column = "";
					for (const auto& i : l.comb) {
						out_column += i + ",";
					}
					out_column.erase( --out_column.end() );
					res_lamp->COMB = out_column;
					lamp.push_back(res_lamp);
				}
			} else {
			  printLOG("No significant SNP combinations were detected.\n");
			}
		}
	} catch (std::string &msg) {
		error("Lamp error. Please check setting. " + msg + "\n");
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
	lamplink_format["U"]=12;
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
	string loutputfile;
	if (par::fastwy) {
		loutputfile=par::output_file_name + ".fastwy";
	} else {
		loutputfile=par::output_file_name + ".lamp";
	}
	lASC.open(loutputfile.c_str(),ios::out);
	ASC.open(outputfile.c_str(),ios::out);
	printf("Plink::outputLamplink\n");
	//make_header
	ASC << setw(lamplink_format["CHR"]) <<"CHR" <<" "
		<< setw(lamplink_format["SNP"]) << "SNP" << " "
		<< setw(lamplink_format["A1"]) << "A1" << " "
		<< setw(lamplink_format["A2"]) << "A2" << " "
		<< setw(lamplink_format["TEST"]) << "TEST" << " ";
	if (!par::assoc_utest) {
		ASC << setw(lamplink_format["AFF"]) << "AFF" << " "
			<< setw(lamplink_format["UNAFF"]) << "UNAFF" << " " ;
		if ( ! par::fisher_test )
			ASC << setw(lamplink_format["CHISQ"]) << "CHISQ" << " "
				<< setw(lamplink_format["DF"]) << "DF" << " ";
		ASC << setw(lamplink_format["P"]) << "P" << " "
			<<setw(lamplink_format["OR"])<< "OR"<< " ";
		if (par::display_ci)
			ASC <<setw(lamplink_format["LXX"]) << string("L"+dbl2str(par::ci_level*100)) <<" "
				<<setw(lamplink_format["UXX"])<<string("U"+dbl2str(par::ci_level*100)) <<" ";
	} else {
		ASC << setw(lamplink_format["U"]) << "U" << " ";
		ASC << setw(lamplink_format["P"]) << "P" << " ";
	}

	//make_header lampfile
	lASC << setw(lamp_format["COMBID"]) <<"COMBID" <<" "
		<<setw(lamp_format["Raw_P"])<< "Raw_P" <<" "
		<<setw(lamp_format["Adjusted_P"])<<"Adjusted_P" <<" ";
	if (par::display_ci)
		lASC	<<setw(lamp_format["OR"])<<"OR" << " "
			<<setw(lamp_format["LXX"]) << string("L"+dbl2str(par::ci_level*100)) <<" "
			<<setw(lamp_format["UXX"])<<string("U"+dbl2str(par::ci_level*100)) <<" ";
	lASC << "COMB" <<"\n"; 
	map<string,string> SNP_COMB ;
	vector<Lamp *>::iterator glamp=lamp.begin();
	int n=1;
	//make 
	while(glamp != lamp.end())
	{
		SNP_COMB[string("COMB"+int2str(n))]=(*glamp)->COMB;
		ASC<<setw(lamplink_format["COMB"])<<string("COMB"+int2str(n))<<" ";
		//make lamp file
		lASC << setw(lamp_format["COMBID"]) <<"COMB"+int2str(n) <<" "
			<<setw(lamp_format["Raw_P"])<< (*glamp)->Raw_P_value <<" "
			<<setw(lamp_format["Adjusted_P"])<< (*glamp)->Ajusted_P_value<<" ";
		if (par::display_ci)
		{
			double lOR;
			double SE;
			double ci_zt ;
//			fprintf(stderr,"%d,%d,%d,%d",(*glamp)->n11,(*glamp)->n00,(*glamp)->n10,(*glamp)->n01);
			if((*glamp)->n10==0 || (*glamp)->n01==0)
			{
				lASC    <<setw(lamp_format["OR"])<<"INF" <<" "
					<<setw(lamp_format["LXX"]) <<"NA" <<" "
					<< setw(lamp_format["UXX"])<<"NA"<<" ";
			}else{
				lASC <<setw(lamp_format["OR"])<< dbl2str( ((double)(*glamp)->n11 * (double)(*glamp)->n00 )
						/((double)(*glamp)->n10 * (double)(*glamp)->n01)) <<" ";
				if((*glamp)->n11==0 || (*glamp)->n00==0){
					lASC <<setw(lamp_format["LXX"]) <<"NA" <<" "<< setw(lamp_format["UXX"])<<"NA"<<" ";
				}else{
					SE=sqrt(1/(double)(*glamp)->n11 + 1/(double)(*glamp)->n10 + 1/(double)(*glamp)->n01 + 1/(double)(*glamp)->n00 );
					ci_zt = ltqnorm( 1 - (1 - par::ci_level) / 2  );
					lOR=log(((double)(*glamp)->n11 * (double)(*glamp)->n00)/((double)(*glamp)->n10 * (double)(*glamp)->n01));
					lASC <<setw(lamp_format["LXX"]) <<dbl2str(exp( lOR - par::ci_zt * SE )) <<" "
						<< setw(lamp_format["UXX"])<<dbl2str(exp( lOR + par::ci_zt * SE )) <<" ";
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
		ASC << setw(lamplink_format["CHR"]) << int2str((*glampassoc)->chr) << " "
			<<setw(lamplink_format["SNP"])<< (*glampassoc)->snpname <<" "
			<<setw(lamplink_format["A1"])<< (*glampassoc)->allele1 <<" "
			<<setw(lamplink_format["A2"])<< (*glampassoc)->allele2 <<" "
			<<setw(lamplink_format["TEST"])<< (*glampassoc)->test <<" ";
		if (!par::assoc_utest) {
			ASC <<setw(lamplink_format["AFF"])<<(*glampassoc)->AFF <<" "
				<<setw(lamplink_format["UNAFF"])<<(*glampassoc)->UNAFF<<" ";
			if(! par::fisher_test)
				ASC<< setw(lamplink_format["CHISQ"]) << (*glampassoc)->chsq
					<<setw(lamplink_format["DF"])<<(*glampassoc)->DF <<" ";
			ASC << setw(lamplink_format["P"]) << (*glampassoc)->P << " "
				<<setw(lamplink_format["OR"])<< (*glampassoc)->OR <<" ";
			if (par::display_ci)
				ASC <<setw(lamplink_format["LXX"]) << (*glampassoc)->CIL <<" "
					<<setw(lamplink_format["UXX"])<< (*glampassoc)->CIU <<" ";
		} else {
			ASC << setw(lamplink_format["U"]) << (*glampassoc)->U << " ";
			ASC << setw(lamplink_format["P"]) << (*glampassoc)->P << " ";
		}
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
			if(res==1) ASC<<setw(lamplink_format["COMB"])<<"1"<<" ";
			else ASC<<setw(lamplink_format["COMB"])<<"0"<<" ";
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
