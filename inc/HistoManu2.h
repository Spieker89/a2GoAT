#ifndef     ___HistoManu2_h___
#define	    ___HistoManu2_h___

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
#include <TH2.h>
#include <TH1.h>
#include "HistoManu.h"

class  HistoManu2 : public TH2F  
{
private:


public:
    HistoManu2();
    HistoManu2(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup,Int_t nbinsy,const  Double_t ylow, const Double_t yup);
    virtual ~HistoManu2();

    virtual Int_t	Fill(Double_t x);

    virtual Int_t	Fill(double fillvalue, double fillvalue2, double timing, double weight=1)   {if(HistoManu::IsPrompt(timing)) return TH2F::Fill(fillvalue,fillvalue2, weight); else if(HistoManu::IsRandom(timing)) return TH2F::Fill(fillvalue,fillvalue2, -1*HistoManu::backgroundSubstractionFactor*weight);}

    virtual HistoManu*    ProjectionX(const char* name = "_px", Int_t firstybin = 0, Int_t lastybin = -1, Option_t* option = "");
    virtual HistoManu*    ProjectionY(const char* name = "_py", Int_t firstxbin = 0, Int_t lastxbin = -1, Option_t* option = "");

   virtual  Int_t Write()   {return TH2F::Write();}



};


#endif

