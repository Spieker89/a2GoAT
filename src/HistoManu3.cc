#include "HistoManu3.h"
#include "GTreeTagger.h"

using namespace std;




HistoManu3::HistoManu3(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup):
       TH3F(name, title, nbinsx, xlow, xup,nbinsy, ylow, yup,nbinsz,zlow,zup) 
{
TH1::Sumw2();
}

HistoManu3::~HistoManu3()
{
}


Int_t   HistoManu3::Fill(const Double_t x)
{
    std::cout << "ERROR: You tried to fill a 3 dim. TH3F with only 1 value." << std::endl;
    return 0;
}

Int_t   HistoManu3::Fill(const Double_t x, Double_t y)
{
    std::cout << "ERROR: You tried to fill a 3 dim. TH3F with only 1 value." << std::endl;
    return 0;
}


HistoManu2*    HistoManu3::ProjectionXY(const char* name, Int_t firstzbin, Int_t lastzbin, Option_t* option)
{
return HistoManu3::ProjectionXY(name, firstzbin, lastzbin, option);


}

							
