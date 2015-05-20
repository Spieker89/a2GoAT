#include <HistoManu.h>

using namespace std;


HistoManu::HistoManu()//constructor
{
}


HistoManu::~HistoManu()//destructor
{
}

double HistoManu::cutPromptMin=1000;
double HistoManu::cutPromptMax=1000;
double HistoManu::cutSideMax=10000;
double HistoManu::cutSideMin=1000;
double HistoManu::backgroundSubstractionFactor=1000;



void HistoManu::InitCutss(double PromptMin, double PromptMax, double RandMin, double RandMax){
	    cutPromptMin    = PromptMin;
	    cutPromptMax    = PromptMax;
    	    cutSideMin    = RandMin;
	    cutSideMax    = RandMax;
	    backgroundSubstractionFactor = (PromptMax - PromptMin)/(2*(cutSideMax - cutSideMin));
		cout << backgroundSubstractionFactor << endl;
}


bool HistoManu::IsPromptt(double value)
{
   if ((value >= cutPromptMin) && (value <= cutPromptMax)){
       return 1;}else{
	   return 0;}
	
}
	
bool HistoManu::IsRandomm(double value)
{

	if((value > (-cutSideMax) && value < (-cutSideMin)) ||(value > cutSideMin && value < cutSideMax)){	
            return 1;}else{
		    return 0;}
}

void HistoManu::FillTH1_timeweighted(TH1F *h1, double fillvalue,double timing)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue);
	}


}

void HistoManu::FillTH1_timeandvalueweighted(TH1F *h1, double fillvalue,double timing, double weight)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,weight*backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,weight);
	}


}

void HistoManu::FillTH2_timeweighted(TH2F *h1, double fillvalue,double fillvalue2,double timing)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2);
	}


}

void HistoManu::FillTH2_timeandvalueweighted(TH2F *h1, double fillvalue,double fillvalue2, double timing, double weight)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,weight*backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2,weight);
	}


}

void HistoManu::FillTH3_timeweighted(TH3F *h1, double fillvalue,double fillvalue2, double fillvalue3,double timing)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3,backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3);
	}


}

void HistoManu::FillTH3_timeandvalueweighted(TH3F *h1, double fillvalue,double fillvalue2,double fillvalue3,double timing, double weight)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3,weight*backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3,weight);
	}


}

								
