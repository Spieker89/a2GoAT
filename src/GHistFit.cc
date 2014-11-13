#include "GHistFit.h"
#include "GTreeMeson.h"
#include "GTreeTagger.h"


GHistFitStruct::GHistFitStruct(const char* name, const char* title, const Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    im(TString(name).Append("_IM"), TString(title).Append(" inv Mass"), 2000, 0, 2000, 48, kFALSE),
    ChiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 10, 48, kFALSE),
    ConfidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE)
{
}

void    GHistFitStruct::Fill(const GFitStruct& fit)
{
    im.Fill(fit.im);
    ChiSq.Fill(fit.ChiSq);
    ConfidenceLevel.Fill(fit.ConfidenceLevel);
}

void    GHistFitStruct::Fill(const GFitStruct& fit, const Double_t taggerTime)
{
    im.Fill(fit.im, taggerTime);
    ChiSq.Fill(fit.ChiSq);
    ConfidenceLevel.Fill(fit.ConfidenceLevel, taggerTime);
}

void    GHistFitStruct::Fill(const GFitStruct& fit, const Double_t taggerTime, const Int_t taggerChannel)
{
    im.Fill(fit.im, taggerTime, taggerChannel);
    ChiSq.Fill(fit.ChiSq);
    ConfidenceLevel.Fill(fit.ConfidenceLevel, taggerTime);
}

void    GHistFitStruct::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
    ChiSq.PrepareWriteList(arr, TString(name).Append("ChiSq").Data());
    ConfidenceLevel.PrepareWriteList(arr, TString(name).Append("ConfidenceLevel").Data());
}








GHistFit3Constraints::GHistFit3Constraints(const char* name, const char* title, const Bool_t IsEtap, const Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(IsEtap),
    raw(TString(name).Append("_raw"), TString(title).Append(" Raw"), kFALSE),
    CutChiSq(TString(name).Append("_CutChiSq"), TString(title).Append(" Cut ChiSq"), kFALSE),
    CutConfidenceLevel(TString(name).Append("_CutConfLev"), TString(title).Append(" Cut ConfLev"), kFALSE),
    CutBoth(TString(name).Append("_CutBoth"), TString(title).Append(" Cut Both"), kFALSE),
    ConfidenceLevelCut(0.1),
    ChiSqCut(10)
{
}

GHistFit3Constraints::~GHistFit3Constraints()
{

}

void    GHistFit3Constraints::Fill(const GFit3Constraints& fit)
{
    raw.Fill(fit.GetResult());
    if(fit.GetResult().ChiSq<ChiSqCut)
    {
        CutChiSq.Fill(fit.GetResult());
        if(fit.GetResult().ConfidenceLevel>ConfidenceLevelCut)
        {
            CutConfidenceLevel.Fill(fit.GetResult());
            CutBoth.Fill(fit.GetResult());
        }
    }
    else if(fit.GetResult().ConfidenceLevel>ConfidenceLevelCut)
        CutConfidenceLevel.Fill(fit.GetResult());
}
void    GHistFit3Constraints::Fill(const GFit3Constraints& fit, const Double_t taggerTime)
{
    raw.Fill(fit.GetResult(), taggerTime);
    if(fit.GetResult().ChiSq<ChiSqCut)
    {
        CutChiSq.Fill(fit.GetResult(), taggerTime);
        if(fit.GetResult().ConfidenceLevel>ConfidenceLevelCut)
        {
            CutConfidenceLevel.Fill(fit.GetResult(), taggerTime);
            CutBoth.Fill(fit.GetResult(), taggerTime);
        }
    }
    else if(fit.GetResult().ConfidenceLevel>ConfidenceLevelCut)
        CutConfidenceLevel.Fill(fit.GetResult(), taggerTime);
}
void    GHistFit3Constraints::Fill(const GFit3Constraints& fit, const Double_t taggerTime, const Int_t taggerChannel)
{
    raw.Fill(fit.GetResult(), taggerTime, taggerChannel);
    if(fit.GetResult().ChiSq<ChiSqCut)
    {
        CutChiSq.Fill(fit.GetResult(), taggerTime, taggerChannel);
        if(fit.GetResult().ConfidenceLevel>ConfidenceLevelCut)
        {
            CutConfidenceLevel.Fill(fit.GetResult(), taggerTime, taggerChannel);
            CutBoth.Fill(fit.GetResult(), taggerTime, taggerChannel);
        }
    }
    else if(fit.GetResult().ConfidenceLevel>ConfidenceLevelCut)
        CutConfidenceLevel.Fill(fit.GetResult(), taggerTime, taggerChannel);
}

void    GHistFit3Constraints::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory("raw");
    raw.PrepareWriteList(folder, TString(name).Append("_raw").Data());
    folder  = arr->GetDirectory("CutChiSq");
    CutChiSq.PrepareWriteList(folder, TString(name).Append("CutChiSq").Data());
    folder  = arr->GetDirectory("CutConfidenceLevel");
    CutConfidenceLevel.PrepareWriteList(folder, TString(name).Append("CutConfidenceLevel").Data());
    folder  = arr->GetDirectory("CutBoth");
    CutBoth.PrepareWriteList(folder, TString(name).Append("CutBoth").Data());
}





/*

GHistFit4Constraints::GHistFit4Constraints(const Bool_t IsEtap) :
    GHistFit3Constraints(6, 4, 0)
{
}

GHistFit4Constraints::~GHistFit4Constraints()
{
}

void    GHistFit4Constraints::InitFit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
   GHistFit3Constraints::InitFit(meson);

   Int_t    index[6] = {0,1,2,3,4,5};
   fitter.AddSubMissMassConstraint(beamAndTarget, 6, index, MASS_PROTON);
}

Bool_t    GHistFit4Constraints::Fit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
    InitFit(meson, beamAndTarget);

    return GHistFit3Constraints::Fit(meson);
}*/
