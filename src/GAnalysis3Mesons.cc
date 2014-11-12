#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const Bool_t IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(IsEtap),
    hist_raw(TString(name).Append("_Raw"), TString(title).Append(" Raw Data")),
    hist_SubImCut(TString(name).Append("_SubImCut"), TString(title).Append(" Sub inv. Mass Cut")),
    hist_MmCut(TString(name).Append("_MmCut"), TString(title).Append(" Sub mis. Mass Cut"))
{
    if(IsEtap)
        SetCutSubIM(0, 497, 697);
    else
        SetCutSubIM(0, 110, 160);
    SetCutSubIM(1, 110, 160);
    SetCutSubIM(2, 110, 160);

    SetCutMM(838, 1038);
}

GAnalysis3Mesons::~GAnalysis3Mesons()
{

}

void   GAnalysis3Mesons::CalcResult()
{
    hist_raw.CalcResult();
    hist_SubImCut.CalcResult();
    hist_MmCut.CalcResult();
}

Bool_t GAnalysis3Mesons::Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Bool_t  found = kFALSE;

    hist_raw.Fill(meson, tagger, CreateHistogramsForTaggerBinning);

    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
       (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
       (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        hist_SubImCut.Fill(meson, tagger, CreateHistogramsForTaggerBinning);

        Double_t    mm;
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
            if(mm>cutMM[0] && mm<cutMM[1])
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MmCut.Fill(meson.Particle(0).M(), mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MmCut.Fill(meson.Particle(0).M(), mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));
                found = kTRUE;
            }
        }

    }
    return found;
}

void    GAnalysis3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        GHistWriteList* folder  = arr->GetDirectory("raw");
        hist_raw.PrepareWriteList(folder, TString(name).Append("_raw").Data());
        folder  = arr->GetDirectory("SubIM_Cut");
        hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
        folder  = arr->GetDirectory("MM_Cut");
        hist_MmCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
    }
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    hist_MmCut.Reset(option);
}

void    GAnalysis3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MmCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GAnalysis3Mesons::SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max)
{
    cutSubIM[2*subNumber] = min;
    cutSubIM[(2*subNumber)+1] = max;
}

void    GAnalysis3Mesons::SetCutMM(const Double_t min, const Double_t max)
{
    cutMM[0] = min;
    cutMM[1] = max;
}











GAnalysis3MesonsProton::GAnalysis3MesonsProton(const char* name, const char* title, const Bool_t IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    hist_meson(name, title, IsEtap, kFALSE),
    fit_meson(TString(name).Append("_fit").Data(), TString(title).Append(" kin. Fit").Data(), IsEtap, kFALSE),
    check_meson_proton(TString(name).Append("_checkProton").Data(), TString(title).Append(" Check Proton").Data(), kFALSE),
    hist_meson_proton(TString(name).Append("_proton").Data(), TString(title).Append(" Proton").Data(), IsEtap, kFALSE),
    fit_meson_proton(TString(name).Append("_proton_fit").Data(), TString(title).Append(" Proton kin. Fit").Data(), IsEtap, kFALSE)
{

}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void   GAnalysis3MesonsProton::CalcResult()
{
    hist_meson.CalcResult();
    fit_meson.CalcResult();
    check_meson_proton.CalcResult();
    hist_meson_proton.CalcResult();
    fit_meson_proton.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    if(proton.GetNParticles()>0)
    {
        if(check_meson_proton.Check(meson, proton, tagger) == kTRUE)
        {
            if(hist_meson_proton.Fill(meson, tagger) == kTRUE)
                fit_meson_proton.Fit(meson, tagger, CreateHistogramsForTaggerBinning);
            return;
        }
        if(hist_meson.Fill(meson, tagger) == kTRUE)
            fit_meson.Fit(meson, tagger, CreateHistogramsForTaggerBinning);
    }
}

void    GAnalysis3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        GHistWriteList* folder  = arr->GetDirectory(TString(name));
        GHistWriteList* subFolder  = folder->GetDirectory("WithoutProton");
        hist_meson.PrepareWriteList(subFolder, name);
        subFolder  = folder->GetDirectory("WithProton");
        hist_meson_proton.PrepareWriteList(subFolder, TString(name).Append("_proton").Data());
    }
    else
    {
        if(hist_meson.IsEtap()==kTRUE)
        {
            GHistWriteList* folder  = arr->GetDirectory("etap");
            GHistWriteList* subFolder  = folder->GetDirectory("WithoutProton");
            GHistWriteList* subsubFolder  = subFolder->GetDirectory("etap");
            hist_meson.PrepareWriteList(subsubFolder, "etap");
            subsubFolder  = subFolder->GetDirectory("fit");
            fit_meson.PrepareWriteList(subsubFolder, "etap_fit");
            subFolder  = folder->GetDirectory("WithProton");
            subsubFolder  = subFolder->GetDirectory("checkProton");
            check_meson_proton.PrepareWriteList(subsubFolder, "etap_proton_fit");
            subsubFolder  = subFolder->GetDirectory("etap");
            hist_meson_proton.PrepareWriteList(subsubFolder, "etap_proton");
            subsubFolder  = subFolder->GetDirectory("fit");
            fit_meson_proton.PrepareWriteList(subsubFolder, "etap_proton_fit");
        }
        else
        {
            GHistWriteList* folder  = arr->GetDirectory("eta");
            GHistWriteList* subFolder  = folder->GetDirectory("WithoutProton");
            GHistWriteList* subsubFolder  = subFolder->GetDirectory("eta");
            hist_meson.PrepareWriteList(subsubFolder, "eta");
            subsubFolder  = subFolder->GetDirectory("fit");
            fit_meson.PrepareWriteList(subsubFolder, "eta_fit");
            subFolder  = folder->GetDirectory("WithProton");
            subsubFolder  = subFolder->GetDirectory("checkProton");
            check_meson_proton.PrepareWriteList(subsubFolder, "eta_proton_fit");
            subsubFolder  = subFolder->GetDirectory("eta");
            hist_meson_proton.PrepareWriteList(subsubFolder, "eta_proton");
            subsubFolder  = subFolder->GetDirectory("fit");
            fit_meson_proton.PrepareWriteList(subsubFolder, "eta_proton_fit");
        }
    }
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    hist_meson.Reset(option);
    fit_meson.Reset(option);
    check_meson_proton.Reset(option);
    hist_meson_proton.Reset(option);
    fit_meson_proton.Reset(option);
}

void    GAnalysis3MesonsProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_meson.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit_meson.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_meson_proton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit_meson_proton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void	GAnalysis3MesonsProton::SetHistMeson(const Double_t sub0_min, const Double_t sub0_max,
                                             const Double_t sub1_min, const Double_t sub1_max,
                                             const Double_t sub2_min, const Double_t sub2_max,
                                             const Double_t mm_min, const Double_t mm_max)
{
    hist_meson.SetCutSubIM(0, sub0_min, sub0_max);
    hist_meson.SetCutSubIM(1, sub1_min, sub1_max);
    hist_meson.SetCutSubIM(2, sub2_min, sub2_max);
    hist_meson.SetCutMM(mm_min, mm_max);
}

void    GAnalysis3MesonsProton::SetFitMeson(const Double_t fit3_CutConfidenceLevel, const Double_t fit4_CutConfidenceLevel)
{
    fit_meson.SetConfidenceLevelCut(fit3_CutConfidenceLevel, fit4_CutConfidenceLevel);
}

void	GAnalysis3MesonsProton::SetCheckProton(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity)
{
    check_meson_proton.SetCuts(maxProtonAngleDiff, minCoplanarity, maxCoplanarity);
}

void	GAnalysis3MesonsProton::SetHistMesonProton(const Double_t sub0_min, const Double_t sub0_max,
                                                   const Double_t sub1_min, const Double_t sub1_max,
                                                   const Double_t sub2_min, const Double_t sub2_max,
                                                   const Double_t mm_min, const Double_t mm_max)
{
    hist_meson_proton.SetCutSubIM(0, sub0_min, sub0_max);
    hist_meson_proton.SetCutSubIM(1, sub1_min, sub1_max);
    hist_meson_proton.SetCutSubIM(2, sub2_min, sub2_max);
    hist_meson_proton.SetCutMM(mm_min, mm_max);
}

void    GAnalysis3MesonsProton::SetFitMesonProton(const Double_t fit3_CutConfidenceLevel, const Double_t fit4_CutConfidenceLevel)
{
    fit_meson_proton.SetConfidenceLevelCut(fit3_CutConfidenceLevel, fit4_CutConfidenceLevel);
}
