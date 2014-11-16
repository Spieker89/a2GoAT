#include "GHistFit.h"
#include "GTreeMeson.h"
#include "GTreeTagger.h"


GHistFitPullParticle::GHistFitPullParticle(const char* name, const char* title, const Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    Pull_Px(TString(name).Append("Pull_Px"), TString(title).Append("Pull_Px"), 1000, -5, 5, 48, kFALSE),
    Pull_Py(TString(name).Append("Pull_Py"), TString(title).Append("Pull_Py"), 1000, -5, 5, 48, kFALSE),
    Pull_Pz(TString(name).Append("Pull_Pz"), TString(title).Append("Pull_Pz"), 1000, -5, 5, 48, kFALSE),
    Pull_E(TString(name).Append("Pull_E"), TString(title).Append("Pull_E"), 1000, -5, 5, 48, kFALSE)
{

}

void    GHistFitPullParticle::Fill(const Double_t* pulls)
{
    Pull_Px.Fill(pulls[0]);
    Pull_Py.Fill(pulls[0]);
    Pull_Pz.Fill(pulls[0]);
    Pull_E.Fill(pulls[0]);
}

void    GHistFitPullParticle::Fill(const Double_t* pulls, const Double_t taggerTime)
{
    Pull_Px.Fill(pulls[0], taggerTime);
    Pull_Py.Fill(pulls[0], taggerTime);
    Pull_Pz.Fill(pulls[0], taggerTime);
    Pull_E.Fill(pulls[0], taggerTime);
}

void    GHistFitPullParticle::Fill(const Double_t* pulls, const Double_t taggerTime, const Int_t taggerChannel)
{
    Pull_Px.Fill(pulls[0], taggerTime, taggerChannel);
    Pull_Py.Fill(pulls[0], taggerTime, taggerChannel);
    Pull_Pz.Fill(pulls[0], taggerTime, taggerChannel);
    Pull_E.Fill(pulls[0], taggerTime, taggerChannel);
}

void    GHistFitPullParticle::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    Pull_Px.PrepareWriteList(arr, TString(name).Append("Pull_Px").Data());
    Pull_Py.PrepareWriteList(arr, TString(name).Append("Pull_Py").Data());
    Pull_Pz.PrepareWriteList(arr, TString(name).Append("Pull_Pz").Data());
    Pull_E.PrepareWriteList(arr, TString(name).Append("Pull_E").Data());
}

void    GHistFitPullParticle::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    Pull_Px.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Pull_Py.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Pull_Pz.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Pull_E.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}







GHistFitPull6Photons::GHistFitPull6Photons(const char* name, const char* title, const Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    p0(TString(name).Append("Photon0"), TString(title).Append("Photon0"), kFALSE),
    p1(TString(name).Append("Photon1"), TString(title).Append("Photon1"), kFALSE),
    p2(TString(name).Append("Photon2"), TString(title).Append("Photon2"), kFALSE),
    p3(TString(name).Append("Photon3"), TString(title).Append("Photon3"), kFALSE),
    p4(TString(name).Append("Photon4"), TString(title).Append("Photon4"), kFALSE),
    p5(TString(name).Append("Photon5"), TString(title).Append("Photon5"), kFALSE)
{

}

void    GHistFitPull6Photons::Fill(const GFitStruct& fit)
{
    p0.Fill(&fit.PullPhotons[0]);
    p1.Fill(&fit.PullPhotons[4]);
    p2.Fill(&fit.PullPhotons[8]);
    p3.Fill(&fit.PullPhotons[12]);
    p2.Fill(&fit.PullPhotons[16]);
    p3.Fill(&fit.PullPhotons[24]);
}

void    GHistFitPull6Photons::Fill(const GFitStruct& fit, const Double_t taggerTime)
{
    p0.Fill(&fit.PullPhotons[0], taggerTime);
    p1.Fill(&fit.PullPhotons[4], taggerTime);
    p2.Fill(&fit.PullPhotons[8], taggerTime);
    p3.Fill(&fit.PullPhotons[12], taggerTime);
    p2.Fill(&fit.PullPhotons[16], taggerTime);
    p3.Fill(&fit.PullPhotons[24], taggerTime);
}

void    GHistFitPull6Photons::Fill(const GFitStruct& fit, const Double_t taggerTime, const Int_t taggerChannel)
{
    p0.Fill(&fit.PullPhotons[0], taggerTime, taggerChannel);
    p1.Fill(&fit.PullPhotons[4], taggerTime, taggerChannel);
    p2.Fill(&fit.PullPhotons[8], taggerTime, taggerChannel);
    p3.Fill(&fit.PullPhotons[12], taggerTime, taggerChannel);
    p2.Fill(&fit.PullPhotons[16], taggerTime, taggerChannel);
    p3.Fill(&fit.PullPhotons[24], taggerTime, taggerChannel);
}

void    GHistFitPull6Photons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    p0.PrepareWriteList(arr, TString(name).Append("Photon0").Data());
    p1.PrepareWriteList(arr, TString(name).Append("Photon1").Data());
    p2.PrepareWriteList(arr, TString(name).Append("Photon2").Data());
    p3.PrepareWriteList(arr, TString(name).Append("Photon3").Data());
    p4.PrepareWriteList(arr, TString(name).Append("Photon4").Data());
    p5.PrepareWriteList(arr, TString(name).Append("Photon5").Data());
}

void    GHistFitPull6Photons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    p0.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    p1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    p2.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    p3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    p4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    p5.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}








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
    ChiSq.Fill(fit.ChiSq, taggerTime);
    ConfidenceLevel.Fill(fit.ConfidenceLevel, taggerTime);
}

void    GHistFitStruct::Fill(const GFitStruct& fit, const Double_t taggerTime, const Int_t taggerChannel)
{
    im.Fill(fit.im, taggerTime, taggerChannel);
    ChiSq.Fill(fit.ChiSq, taggerTime);
    ConfidenceLevel.Fill(fit.ConfidenceLevel, taggerTime);
}

void    GHistFitStruct::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
    ChiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
    ConfidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfidenceLevel").Data());
}

void    GHistFitStruct::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    ChiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    ConfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}










GHistFit::GHistFit(const char* name, const char* title, const Bool_t IsEtap, const Bool_t linkHistogram) :
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

GHistFit::~GHistFit()
{

}

void    GHistFit::Fill(const GFit3Constraints& fit)
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
void    GHistFit::Fill(const GFit3Constraints& fit, const Double_t taggerTime)
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
void    GHistFit::Fill(const GFit3Constraints& fit, const Double_t taggerTime, const Int_t taggerChannel)
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

void    GHistFit::PrepareWriteList(GHistWriteList* arr, const char* name)
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

void    GHistFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CutChiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CutConfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CutBoth.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}




/*

GHistFit4Constraints::GHistFit4Constraints(const Bool_t IsEtap) :
    GHistFit(6, 4, 0)
{
}

GHistFit4Constraints::~GHistFit4Constraints()
{
}

void    GHistFit4Constraints::InitFit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
   GHistFit::InitFit(meson);

   Int_t    index[6] = {0,1,2,3,4,5};
   fitter.AddSubMissMassConstraint(beamAndTarget, 6, index, MASS_PROTON);
}

Bool_t    GHistFit4Constraints::Fit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
    InitFit(meson, beamAndTarget);

    return GHistFit::Fit(meson);
}*/
