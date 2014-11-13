#include "GCheckProton.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"



GCheckProtonHist::GCheckProtonHist(const char* name, const char* title, Bool_t linkHistogram) :
    protonAngeDiff(TString(name).Append("_ProtonAngleDiff"), TString(title).Append(" Proton Angle Diff."), 1000, 0, 100, 48, linkHistogram),
    protonCoplanarity(TString(name).Append("_Coplanarity"), TString(title).Append(" Coplanarity"), 3600, 0, 360, 48, linkHistogram)
{

}

GCheckProtonHist::~GCheckProtonHist()
{

}

void    GCheckProtonHist::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        protonAngeDiff.PrepareWriteList(arr, TString(name).Append("_prAng").Data());
        protonCoplanarity.PrepareWriteList(arr, TString(name).Append("_copl").Data());
    }
}






GCheckProton::GCheckProton(const char* name, const char* title, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    raw(TString(name).Append("_Raw"), TString(title).Append(" Raw"), kFALSE),
    cutProtonAngle(TString(name).Append("_cutProtonAngle"), TString(title).Append(" Cut Proton Angle"), kFALSE),
    cutCoplanarity(TString(name).Append("_cutCoplanarity"), TString(title).Append(" Cut Coplanarity"), kFALSE),
    cutBoth(TString(name).Append("_cutBoth"), TString(title).Append(" Cut Both"), kFALSE),
    CutProtonAngleDiff(10)
{
    CutProtonCoplanarity[0] = 160;
    CutProtonCoplanarity[1] = 200;
}

GCheckProton::~GCheckProton()
{

}

void   GCheckProton::CalcResult()
{
    raw.CalcResult();
    cutProtonAngle.CalcResult();
    cutCoplanarity.CalcResult();
    cutBoth.CalcResult();
}

Bool_t  GCheckProton::Check(const GTreeMeson& meson, const GTreeParticle& proton, const TLorentzVector& beamAndTarget, const Double_t taggerTime)
{
    if(proton.GetNParticles()!=1)
        return kFALSE;

    Bool_t      passed          = kFALSE;
    Double_t    helpAngleDiff   = TMath::RadToDeg()*(beamAndTarget-meson.Particle(0)).Angle(proton.Particle(0).Vect());
    Double_t    helpCoplanarity = TMath::RadToDeg()*TMath::Abs(meson.Particle(0).Phi()-proton.Particle(0).Phi());

    raw.Fill(helpAngleDiff, taggerTime);
    if(helpCoplanarity>CutProtonCoplanarity[0] && helpCoplanarity<CutProtonCoplanarity[1])
    {
        cutCoplanarity.Fill(helpAngleDiff, taggerTime);
        if(helpAngleDiff<CutProtonAngleDiff)
        {
            passed = kTRUE;
            cutProtonAngle.Fill(helpAngleDiff, taggerTime);
            cutBoth.Fill(helpAngleDiff, taggerTime);
        }
    }
    else if(helpAngleDiff<CutProtonAngleDiff)
        cutProtonAngle.Fill(helpAngleDiff, taggerTime);

    return passed;
}

void    GCheckProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        GHistWriteList* folder  = arr->GetDirectory("raw");
        raw.PrepareWriteList(folder, TString(name).Append("_raw").Data());

        folder  = arr->GetDirectory("cutCoplanarity");
        cutCoplanarity.PrepareWriteList(folder, TString(name).Append("_cutCopl").Data());


        folder  = arr->GetDirectory("cutProtonAngle");
        cutProtonAngle.PrepareWriteList(folder, TString(name).Append("cutPrAng").Data());


        folder  = arr->GetDirectory("cutBoth");
        cutBoth.PrepareWriteList(folder, TString(name).Append("_cutBoth").Data());
    }
}

void    GCheckProton::Reset(Option_t* option)
{
    raw.Reset(option);
    cutProtonAngle.Reset(option);
    cutCoplanarity.Reset(option);
    cutBoth.Reset(option);
}

void    GCheckProton::SetCuts(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity)
{
    CutProtonAngleDiff      = maxProtonAngleDiff;
    CutProtonCoplanarity[0] = minCoplanarity;
    CutProtonCoplanarity[1] = maxCoplanarity;
}
