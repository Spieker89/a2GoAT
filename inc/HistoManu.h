#ifndef     ___HistoManu_h___
#define	    ___HistoManu_h___

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
#include <TH1.h>


class  HistoManu : public TH1F
{
protected:


public:


    static  Double_t    cutPromptMin;
    static  Double_t    cutPromptMax;
    static  Double_t    cutSideMax;
    static  Double_t    cutSideMin;
    static  Double_t    backgroundSubstractionFactor;


    HistoManu();
    HistoManu(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
    virtual ~HistoManu();

    static void    InitCuts(double PromptMin, double PromptMax, double RandMin, double RandMax);
    static  Bool_t  IsPrompt(const Double_t value);
    static  Bool_t  IsRandom(const Double_t value);

    virtual Int_t	Fill(double fillvalue, double timing, double weight=1)   {if(IsPrompt(timing)) return TH1F::Fill(fillvalue,weight); else if(IsRandom(timing)) return TH1F::Fill(fillvalue, -1*backgroundSubstractionFactor*weight);}





};


#endif

