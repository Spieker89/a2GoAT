#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    fit3(_IsEtap),
    hist_SubImCut(TString(name).Append("_SubImCut"), TString(title).Append(" Sub inv. Mass Cut"), kFALSE),
    hist_SubImCut_fit3(TString(name).Append("_SubImCut_fit3"), TString(title).Append(" Sub inv. Mass Cut fit3"), _IsEtap, kFALSE),
    hist_MmCut(TString(name).Append("_MmCut"), TString(title).Append(" Sub mis. Mass Cut"), kFALSE),
    hist_MmCut_fit3(TString(name).Append("_MmCut_fit3"), TString(title).Append(" Sub mis. Mass Cut fit3"), _IsEtap, kFALSE)
{
    if(_IsEtap==kTRUE)
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
    hist_SubImCut.CalcResult();
    hist_MmCut.CalcResult();
}

Bool_t  GAnalysis3Mesons::Fill(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
        (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
        (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        mm  = (beamAndTarget-meson.Particle(0)).M();

        hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2);
        if(fit3.Fit(meson)==kTRUE)
            hist_SubImCut_fit3.Fill(fit3);

        if(mm>cutMM[0] && mm<cutMM[1])
        {
            hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2);
            return kTRUE;
        }
    }
    return kFALSE;
}
Bool_t  GAnalysis3Mesons::Fill(const GTreeMeson& meson, const TLorentzVector& beamAndTarget, const Double_t taggerTime)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
        (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
        (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        mm  = (beamAndTarget-meson.Particle(0)).M();

        hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, taggerTime);
        if(fit3.Fit(meson)==kTRUE)
            hist_SubImCut_fit3.Fill(fit3, taggerTime);

        if(mm>cutMM[0] && mm<cutMM[1])
        {
            hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, taggerTime);
            return kTRUE;
        }
    }
    return kFALSE;
}
Bool_t  GAnalysis3Mesons::Fill(const GTreeMeson& meson, const TLorentzVector& beamAndTarget, const Double_t taggerTime, const Int_t taggerChannel)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
        (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
        (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        mm  = (beamAndTarget-meson.Particle(0)).M();

        hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, taggerTime, taggerChannel);
        if(fit3.Fit(meson)==kTRUE)
            hist_SubImCut_fit3.Fill(fit3, taggerTime, taggerChannel);

        if(mm>cutMM[0] && mm<cutMM[1])
        {
            hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, taggerTime, taggerChannel);
            return kTRUE;
        }
    }
    return kFALSE;
}

Bool_t GAnalysis3Mesons::Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Bool_t  found = kFALSE;

    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
        (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
        (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();

            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            else
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));
            bool    fitDone = false;
            if(fit3.Fit(meson)==kTRUE)
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_SubImCut_fit3.Fill(fit3, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_SubImCut_fit3.Fill(fit3, tagger.GetTagged_t(i));
                fitDone = true;
            }

            if(mm>cutMM[0] && mm<cutMM[1])
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                if(fitDone)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_MmCut_fit3.Fill(fit3, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        hist_MmCut_fit3.Fill(fit3, tagger.GetTagged_t(i));
                }
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
        GHistWriteList* folder  = arr->GetDirectory("SubIM_Cut");
        hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
        GHistWriteList* fitfolder  = folder->GetDirectory("fit3");
        hist_SubImCut_fit3.PrepareWriteList(fitfolder, TString(name).Append("_fit3").Data());
        folder  = arr->GetDirectory("MM_Cut");
        hist_MmCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
        fitfolder  = folder->GetDirectory("fit3");
        hist_SubImCut_fit3.PrepareWriteList(fitfolder, TString(name).Append("_fit3").Data());
    }
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_SubImCut.Reset(option);
    hist_MmCut.Reset(option);
}

void    GAnalysis3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
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
    check_meson_proton(TString(name).Append("_checkProton").Data(), TString(title).Append(" Check Proton").Data(), kFALSE),
    hist_meson_proton(TString(name).Append("_proton").Data(), TString(title).Append(" Proton").Data(), IsEtap, kFALSE)
{

}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void   GAnalysis3MesonsProton::CalcResult()
{
    hist_meson.CalcResult();
    check_meson_proton.CalcResult();
    hist_meson_proton.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    if(proton.GetNParticles()>0)
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            if(check_meson_proton.Check(meson, proton, tagger.GetVectorProtonTarget(i), tagger.GetTagged_t(i)) == kTRUE)
            {
                hist_meson_proton.Fill(meson, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                return;
            }
            hist_meson.Fill(meson, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
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
            subFolder  = folder->GetDirectory("WithProton");
            subsubFolder  = subFolder->GetDirectory("checkProton");
            check_meson_proton.PrepareWriteList(subsubFolder, "etap_proton_fit");
            subsubFolder  = subFolder->GetDirectory("etap");
            hist_meson_proton.PrepareWriteList(subsubFolder, "etap_proton");
        }
        else
        {
            GHistWriteList* folder  = arr->GetDirectory("eta");
            GHistWriteList* subFolder  = folder->GetDirectory("WithoutProton");
            GHistWriteList* subsubFolder  = subFolder->GetDirectory("eta");
            hist_meson.PrepareWriteList(subsubFolder, "eta");
            subFolder  = folder->GetDirectory("WithProton");
            subsubFolder  = subFolder->GetDirectory("checkProton");
            check_meson_proton.PrepareWriteList(subsubFolder, "eta_proton_fit");
            subsubFolder  = subFolder->GetDirectory("eta");
            hist_meson_proton.PrepareWriteList(subsubFolder, "eta_proton");
        }
    }
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    hist_meson.Reset(option);
    check_meson_proton.Reset(option);
    hist_meson_proton.Reset(option);
}

void    GAnalysis3MesonsProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_meson.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_meson_proton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
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
    //fit_meson.SetConfidenceLevelCut(fit3_CutConfidenceLevel, fit4_CutConfidenceLevel);
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
    //fit_meson_proton.SetConfidenceLevelCut(fit3_CutConfidenceLevel, fit4_CutConfidenceLevel);
}
