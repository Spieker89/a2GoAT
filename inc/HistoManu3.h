#ifndef     ___HistoManu3_h___
#define	    ___HistoManu3_h___

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
#include <TH3.h>
#include "HistoManu.h"
#include "HistoManu2.h"

class  HistoManu3 : public TH3F
{
private:


public:
    HistoManu3();
    HistoManu3(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup,Int_t nbinsy,const  Double_t ylow, const Double_t yup, Int_t nbinsz, const Double_t zlow, const Double_t zup);
    virtual ~HistoManu3();

    virtual Int_t	Fill(Double_t x);
    virtual Int_t	Fill(Double_t x, Double_t y);

    virtual Int_t	Fill(double fillvalue, double fillvalue2, double fillvalue3, double timing, double weight=1)   {if(HistoManu::IsPrompt(timing))return TH3F::Fill(fillvalue,fillvalue2,fillvalue3, weight); else if(HistoManu::IsRandom(timing)) return TH3F::Fill(fillvalue,fillvalue2,fillvalue3,-1*HistoManu::backgroundSubstractionFactor*weight);}

    virtual HistoManu2*    ProjectionXY(const char* name = "_pxy", Int_t firstzbin = 0, Int_t lastzbin = -1, Option_t* option = "");


   virtual  Int_t Write()   {return TH3F::Write();}



};


#endif

