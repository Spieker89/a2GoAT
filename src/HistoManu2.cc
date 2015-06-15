#include "HistoManu2.h"
#include "GTreeTagger.h"

using namespace std;




HistoManu2::HistoManu2(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup):
       TH2F(name, title, nbinsx, xlow, xup,nbinsy, ylow, yup) 
{
TH1::Sumw2();
}

HistoManu2::~HistoManu2()
{
}


Int_t   HistoManu2::Fill(const Double_t x)
{
    std::cout << "ERROR: You tried to fill a 2 dim. TH2F with only 1 value." << std::endl;
    return 0;
}

HistoManu*    HistoManu2::ProjectionX(const char* name, Int_t firstybin, Int_t lastybin, Option_t* option)
{
return HistoManu2::ProjectionX(name, firstybin, lastybin, option);

    
}	

HistoManu*    HistoManu2::ProjectionY(const char* name, Int_t firstybin, Int_t lastybin, Option_t* option)
{
return HistoManu2::ProjectionY(name, firstybin, lastybin, option);

    
}	

							
