#include "GCheckProton.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"



GCheckProtonHist::GCheckProtonHist(const char* name, const char* title, Bool_t linkHistogram) :
    protonAngeDiff(TString(name).Append("_ProtonAngleDiff"), TString(title).Append(" Proton Angle Diff."), 1000, 0, 100, 48, linkHistogram),
    protonAngeDiffSmalest(TString(name).Append("_ProtonAngleDiffSmalest"), TString(title).Append(" Smalest Proton Angle Diff."), 1000, 0, 100, 48, linkHistogram),
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
        protonAngeDiffSmalest.PrepareWriteList(arr, TString(name).Append("_prAngMin").Data());
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

Bool_t  GCheckProton::Check(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger)
{
    if(proton.GetNParticles()>1)
        return kFALSE;
    if(tagger.GetNTagged() == 0)
        return kFALSE;

    Bool_t      passedAngleDiff   = kFALSE;
    Double_t    helpAngleDiff;
    Double_t    smalestAngleDiff    = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(0)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
    Double_t    smalestAngleDiffTaggerTime     = tagger.GetTagged_t(0);
    Double_t    smalestAngleDiffTaggerBin      = tagger.GetTagged_ch(0);
    Double_t    helpCoplanarity     = TMath::RadToDeg()*TMath::Abs(meson.Particle(0).Phi()-proton.Particle(0).Phi());

    if(helpCoplanarity>CutProtonCoplanarity[0] && helpCoplanarity<CutProtonCoplanarity[1])
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
            if(helpAngleDiff < smalestAngleDiff)
                smalestAngleDiff            = helpAngleDiff;
            raw.Fill(helpAngleDiff);
            cutCoplanarity.Fill(helpAngleDiff);
            if(helpAngleDiff<CutProtonAngleDiff)
            {
                passedAngleDiff = kTRUE;
                cutProtonAngle.Fill(helpAngleDiff);
                cutBoth.Fill(helpAngleDiff);
            }
        }
        raw.Fill(smalestAngleDiff, helpCoplanarity);
        cutCoplanarity.Fill(smalestAngleDiff, helpCoplanarity);
        if(passedAngleDiff == kTRUE)
        {
            cutProtonAngle.Fill(smalestAngleDiff, helpCoplanarity);
            cutBoth.Fill(smalestAngleDiff, helpCoplanarity);
            return kTRUE;
        }
    }
    else
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
            if(helpAngleDiff < smalestAngleDiff)
            {
                smalestAngleDiff            = helpAngleDiff;
                smalestAngleDiffTaggerTime  = tagger.GetTagged_t(i);
                smalestAngleDiffTaggerBin   = tagger.GetTagged_ch(i);
            }
            raw.Fill(helpAngleDiff);
            if(helpAngleDiff<CutProtonAngleDiff)
            {
                passedAngleDiff = kTRUE;
                cutProtonAngle.Fill(helpAngleDiff);
            }
        }
        raw.Fill(smalestAngleDiff, helpCoplanarity);
        if(passedAngleDiff == kTRUE)
            cutProtonAngle.Fill(smalestAngleDiff, helpCoplanarity);
    }
    return kFALSE;
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
