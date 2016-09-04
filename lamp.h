

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

#ifndef __LAMP_H__
#define __LAMP_H__

//#include <Python.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <new>

#include "zed.h"

class Lamp 
{
	public:
	int n11;
	int n10;
	int n01;
	int n00;
	string Rank;
	string Raw_P_value;
	string Ajusted_P_value;
	string COMB;

};

class LampAssoc
{
	public:
	int chr;
	string snpname;
	string allele1;
	string allele2;
	string test;
	string AFF;
	string UNAFF;
	string chsq;
	string DF;
	string P;
	string OR;
	string CIL;
	string CIU;
};


class readLamplink
{
	public:
	vector <string> header;
	vector <string> item; 
	// len(item)=11 (--model-dom or --model-rec) or (--fisher --ci )
	// len(item)=9  (--model-dom or --model-rec) --fisher
	// len(item)=13 (--model-dom or --model-rec) --ci X.XX
	vector <string> combination; //COMB1 COMB2 COMB3....etc
};

class readLamp
{
	public:
	vector <string> header;
	vector <string> item; //first item is COMBID
	//len(item)==5 --lamp
	//len(item)==7 --lamp --ci X.XX
	string combination;

};

#endif
