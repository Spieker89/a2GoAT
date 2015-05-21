#ifndef     ___HistoManu_h___
#define	    ___HistoManu_h___

#include "TH1.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 
#include <sstream>

#define GHBS_folderName         "BackgroundSubstraction"
#define GHBS_randFolderName     "RandomWindows"
#define GHBS_randNameSuffix     "_Rand"
#define GHBS_randTitleSuffix    " Rand "
#define GHBS_randSumNameSuffix  "_RandSum"
#define GHBS_randSumTitleSuffix " RandSum"
#define GHBS_promptNameSuffix   "_Prompt"
#define GHBS_promptTitleSuffix  " Prompt"


#include <TDirectory.h>


class  HistoManu : public TH1D
{
private:

    static  Double_t    cutPromptMin;
    static  Double_t    cutPromptMax;
    static  Double_t    backgroundSubstractionFactor;


public:
    HistoManu(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
    virtual ~HistoManu();

    static  void    InitCuts(Double_t PromptMin, Double_t PromptMax, Double_t RandMin, Double_t RandMax);
    static  Bool_t  IsPrompt(const Double_t value);
    static  Bool_t  IsRandom(const Double_t value);

    virtual Int_t	Fill(double fillvalue, double timing, Double_t weight=1)   {if(IsPrompt(timing)) TH1D::Fill(fillvalue,weight); else if(IsRandom(timing)) TH1D::Fill(fillvalue, weight*backgroundSubstractionFactor);}

};


#endif

