#include "HistoManu.h"
#include "GTreeTagger.h"

using namespace std;



Double_t HistoManu::cutPromptMin=1000;
Double_t HistoManu::cutPromptMax=1000;
Double_t HistoManu::cutSideMax=10000;
Double_t HistoManu::cutSideMin=1000;
Double_t HistoManu::backgroundSubstractionFactor=1000;



void	HistoManu::InitCuts(const Double_t PromptMin,const Double_t PromptMax,const Double_t RandMin, const Double_t RandMax){
	    cutPromptMin    = PromptMin;
	    cutPromptMax    = PromptMax;
    	    cutSideMin    = RandMin;
	    cutSideMax    = RandMax;
	    backgroundSubstractionFactor = (PromptMax - PromptMin)/(2*(cutSideMax - cutSideMin));
}


Bool_t    HistoManu::IsPrompt(const Double_t value)
{
   if ((value >= cutPromptMin) && (value <= cutPromptMax))
       return kTRUE;
   return kFALSE;
}

Bool_t    HistoManu::IsRandom(const Double_t value)
{

	if((value > (-cutSideMax) && value < (-cutSideMin)) ||(value > cutSideMin && value < cutSideMax)){	
            return kTRUE;}else{
		    return kFALSE;}
}


HistoManu::HistoManu(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup):
    TH1F(name, title, nbinsx, xlow, xup)
{
TH1::Sumw2();
}

HistoManu::~HistoManu()
{
}



								
