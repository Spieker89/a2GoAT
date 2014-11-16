#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    fit3(_IsEtap),
    fit4(_IsEtap),
    fit3Beam(_IsEtap),
    fit4Beam(_IsEtap),
    hist_SubImCut(TString(name).Append("_SubImCut"), TString(title).Append(" Sub inv. Mass Cut"), kFALSE),
    hist_SubImCut_fit3(TString(name).Append("_SubImCut_fit3"), TString(title).Append(" Sub inv. Mass Cut fit3"), _IsEtap, kFALSE),
    hist_SubImCut_fit4(TString(name).Append("_SubImCut_fit4"), TString(title).Append(" Sub inv. Mass Cut fit4"), _IsEtap, kFALSE),
    hist_SubImCut_fit3Beam(TString(name).Append("_SubImCut_fit3Beam"), TString(title).Append(" Sub inv. Mass Cut fit3Beam"), _IsEtap, kFALSE),
    hist_SubImCut_fit4Beam(TString(name).Append("_SubImCut_fit4Beam"), TString(title).Append(" Sub inv. Mass Cut fit4Beam"), _IsEtap, kFALSE),
    hist_MmCut(TString(name).Append("_MmCut"), TString(title).Append(" Sub mis. Mass Cut"), kFALSE),
    hist_MmCut_fit3(TString(name).Append("_MmCut_fit3"), TString(title).Append(" Sub mis. Mass Cut fit3"), _IsEtap, kFALSE),
    hist_MmCut_fit4(TString(name).Append("_MmCut_fit4"), TString(title).Append(" Sub mis. Mass Cut fit4"), _IsEtap, kFALSE),
    hist_MmCut_fit3Beam(TString(name).Append("_MmCut_fit3Beam"), TString(title).Append(" Sub mis. Mass Cut fit3Beam"), _IsEtap, kFALSE),
    hist_MmCut_fit4Beam(TString(name).Append("_MmCut_fit4Beam"), TString(title).Append(" Sub mis. Mass Cut fit4Beam"), _IsEtap, kFALSE)
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
    hist_SubImCut_fit3.CalcResult();
    hist_SubImCut_fit4.CalcResult();
    hist_SubImCut_fit3Beam.CalcResult();
    hist_SubImCut_fit4Beam.CalcResult();
    hist_MmCut.CalcResult();
    hist_MmCut_fit3.CalcResult();
    hist_MmCut_fit4.CalcResult();
    hist_MmCut_fit3Beam.CalcResult();
    hist_MmCut_fit4Beam.CalcResult();
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

        bool    fit3Done = false;
        if(fit3.Fit(meson)==kTRUE)
        {
            hist_SubImCut_fit3.Fill(fit3, taggerTime);
            fit3Done = true;
        }
        bool    fit4Done = false;
        if(fit4.Fit(meson, beamAndTarget)==kTRUE)
        {
            hist_SubImCut_fit4.Fill(fit4, taggerTime);
            fit4Done = true;
        }
        bool    fit3BeamDone = false;
        if(fit3Beam.Fit(meson, beamAndTarget)==kTRUE)
        {
            hist_SubImCut_fit3Beam.Fill(fit3Beam, taggerTime);
            fit3BeamDone = true;
        }
        bool    fit4BeamDone = false;
        if(fit4Beam.Fit(meson, beamAndTarget)==kTRUE)
        {
            hist_SubImCut_fit4Beam.Fill(fit4Beam, taggerTime);
            fit4BeamDone = true;
        }

        if(mm>cutMM[0] && mm<cutMM[1])
        {
            hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, taggerTime);

            if(fit3Done)
                hist_MmCut_fit3.Fill(fit3, taggerTime);
            if(fit4Done)
                hist_MmCut_fit4.Fill(fit4, taggerTime);
            if(fit3BeamDone)
                hist_MmCut_fit3Beam.Fill(fit3Beam, taggerTime);
            if(fit4BeamDone)
                hist_MmCut_fit4Beam.Fill(fit4Beam, taggerTime);


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

        bool    fit3Done = false;
        if(fit3.Fit(meson)==kTRUE)
        {
            hist_SubImCut_fit3.Fill(fit3, taggerTime, taggerChannel);
            fit3Done = true;
        }
        bool    fit4Done = false;
        if(fit4.Fit(meson, beamAndTarget)==kTRUE)
        {
            hist_SubImCut_fit4.Fill(fit4, taggerTime, taggerChannel);
            fit4Done = true;
        }
        bool    fit3BeamDone = false;
        if(fit3Beam.Fit(meson, beamAndTarget)==kTRUE)
        {
            hist_SubImCut_fit3Beam.Fill(fit3Beam, taggerTime, taggerChannel);
            fit3BeamDone = true;
        }
        bool    fit4BeamDone = false;
        if(fit4Beam.Fit(meson, beamAndTarget)==kTRUE)
        {
            hist_SubImCut_fit4Beam.Fill(fit4Beam, taggerTime, taggerChannel);
            fit4BeamDone = true;
        }

        if(mm>cutMM[0] && mm<cutMM[1])
        {
            hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, taggerTime, taggerChannel);

            if(fit3Done)
                hist_MmCut_fit3.Fill(fit3, taggerTime, taggerChannel);
            if(fit4Done)
                hist_MmCut_fit4.Fill(fit4, taggerTime, taggerChannel);
            if(fit3BeamDone)
                hist_MmCut_fit3Beam.Fill(fit3Beam, taggerTime, taggerChannel);
            if(fit4BeamDone)
                hist_MmCut_fit4Beam.Fill(fit4Beam, taggerTime, taggerChannel);

            return  kTRUE;
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

            bool    fit3Done = false;
            if(fit3.Fit(meson)==kTRUE)
            {
                hist_SubImCut_fit3.Fill(fit3, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                fit3Done = true;
            }
            bool    fit4Done = false;
            if(fit4.Fit(meson, tagger.GetVectorProtonTarget(i))==kTRUE)
            {
                hist_SubImCut_fit4.Fill(fit4, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                fit4Done = true;
            }
            bool    fit3BeamDone = false;
            if(fit3Beam.Fit(meson, tagger.GetVectorProtonTarget(i))==kTRUE)
            {
                hist_SubImCut_fit3Beam.Fill(fit3Beam, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                fit3BeamDone = true;
            }
            bool    fit4BeamDone = false;
            if(fit4Beam.Fit(meson, tagger.GetVectorProtonTarget(i))==kTRUE)
            {
                hist_SubImCut_fit4Beam.Fill(fit4Beam, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                fit4BeamDone = true;
            }

            if(mm>cutMM[0] && mm<cutMM[1])
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MmCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                if(fit3Done)
                    hist_MmCut_fit3.Fill(fit3, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                if(fit4Done)
                    hist_MmCut_fit4.Fill(fit4, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                if(fit3BeamDone)
                    hist_MmCut_fit3Beam.Fill(fit3Beam, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                if(fit4BeamDone)
                    hist_MmCut_fit4Beam.Fill(fit4Beam, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));

                found   = kTRUE;
            }
        }

    }
    return found;
}

void    GAnalysis3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
    GHistWriteList* fitfolder  = folder->GetDirectory("fit3");
    hist_SubImCut_fit3.PrepareWriteList(fitfolder, TString(name).Append("_fit3").Data());
    fitfolder  = folder->GetDirectory("fit4");
    hist_SubImCut_fit4.PrepareWriteList(fitfolder, TString(name).Append("_fit4").Data());
    fitfolder  = folder->GetDirectory("fit3Beam");
    hist_SubImCut_fit3Beam.PrepareWriteList(fitfolder, TString(name).Append("_fit3Beam").Data());
    fitfolder  = folder->GetDirectory("fit4Beam");
    hist_SubImCut_fit4Beam.PrepareWriteList(fitfolder, TString(name).Append("_fit4Beam").Data());
    folder  = arr->GetDirectory("MM_Cut");
    hist_MmCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
    fitfolder  = folder->GetDirectory("fit3");
    hist_SubImCut_fit3.PrepareWriteList(fitfolder, TString(name).Append("_fit3").Data());
    fitfolder  = folder->GetDirectory("fit4");
    hist_SubImCut_fit4.PrepareWriteList(fitfolder, TString(name).Append("_fit4").Data());
    fitfolder  = folder->GetDirectory("fit3Beam");
    hist_SubImCut_fit3Beam.PrepareWriteList(fitfolder, TString(name).Append("_fit3Beam").Data());
    fitfolder  = folder->GetDirectory("fit4Beam");
    hist_SubImCut_fit4Beam.PrepareWriteList(fitfolder, TString(name).Append("_fit4Beam").Data());
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_SubImCut.Reset(option);
    hist_SubImCut_fit3.Reset(option);
    hist_SubImCut_fit4.Reset(option);
    hist_SubImCut_fit3Beam.Reset(option);
    hist_SubImCut_fit4Beam.Reset(option);
    hist_MmCut.Reset(option);
    hist_MmCut_fit3.Reset(option);
    hist_MmCut_fit4.Reset(option);
    hist_MmCut_fit3Beam.Reset(option);
    hist_MmCut_fit4Beam.Reset(option);
}

void    GAnalysis3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit3Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MmCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MmCut_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MmCut_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MmCut_fit3Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MmCut_fit4Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
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
        Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
        Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
        Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

        if((sub_im_0>hist_meson_proton.cutSubIM[0] && sub_im_0<hist_meson_proton.cutSubIM[1]) &&
            (sub_im_1>hist_meson_proton.cutSubIM[2] && sub_im_1<hist_meson_proton.cutSubIM[3]) &&
            (sub_im_2>hist_meson_proton.cutSubIM[4] && sub_im_2<hist_meson_proton.cutSubIM[5]))
        {
            for(int i=0; i<tagger.GetNTagged(); i++)
            {
                check_meson_proton.Check(meson, proton, tagger.GetVectorProtonTarget(i), tagger.GetTagged_t(i));
                hist_meson_proton.Fill(meson, tagger.GetVectorProtonTarget(i), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                return;
            }
        }
    }
    else
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
            hist_meson.Fill(meson, tagger.GetVectorProtonTarget(i), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
    }
}

void    GAnalysis3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

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
