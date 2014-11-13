#include "GHistFit.h"
#include "GTreeMeson.h"
#include "GTreeTagger.h"


GHistFitStruct::GHistFitStruct(const char* name, const char* title, const Bool_t linkHistogram = kTRUE)   :
    im(TString(name).Append("_IM"), TString(title).Append(" inv Mass"), 2000, 0, 2000, 48, linkHistogram),
    ChiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 10, 48, linkHistogram),
    ConfidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, linkHistogram)
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



GHistFit3Constraints::GHistFit3Constraints(const char* name, const char* title, const Bool_t IsEtap, const Bool_t linkHistogram = kTRUE) :
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
