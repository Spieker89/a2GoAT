#ifndef     ___HistoManu_h___
#define	    ___HistoManu_h___

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 
#include <sstream>

class	HistoManu{
protected:
static double cutPromptMin,cutPromptMax,cutSideMax,cutSideMin,backgroundSubstractionFactor;
public:
HistoManu();
virtual ~HistoManu();
static void    InitCutss(double PromptMin, double PromptMax, double RandMin, double RandMax);
static bool   IsPromptt(double value);
static bool    IsRandomm(double value);

static void FillTH1_timeweighted(TH1F *h1, double fillvalue,double timing);
static void FillTH1_timeandvalueweighted(TH1F *h1, double fillvalue,double timing, double weight);
static void FillTH2_timeweighted(TH2F *h1, double fillvalue,double fillvalue2,double timing);
static void FillTH2_timeandvalueweighted(TH2F *h1, double fillvalue,double fillvalue2, double timing, double weight);
static void FillTH3_timeweighted(TH3F *h1, double fillvalue,double fillvalue2, double fillvalue3,double timing);
static void FillTH3_timeandvalueweighted(TH3F *h1, double fillvalue,double fillvalue2,double fillvalue3,double timing, double weight);

};
#endif

