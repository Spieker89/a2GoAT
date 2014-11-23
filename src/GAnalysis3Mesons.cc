#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    hist_raw(TString(name).Append("raw"), TString(title).Append("raw"), kFALSE),
    hist_SubImCut(TString(name).Append("SubImCut"), TString(title).Append("SubImCut"), kFALSE),
    hist_MMCut(TString(name).Append("MMCut"), TString(title).Append("MMCut"), kFALSE),
    fit3(kTRUE),
    fit4(kTRUE),
    fit3Beam(kTRUE),
    fit4Beam(kTRUE),
    fit3_im(TString(name).Append("fit3_im"), TString(title).Append("fit3_im"), 2000, 0, 2000, 48, kFALSE),
    fit3_cs(TString(name).Append("fit3_cs"), TString(title).Append("fit3_cs"), 1000, 0, 1000, 48, kFALSE),
    fit3_cl(TString(name).Append("fit3_cl"), TString(title).Append("fit3_cl"), 1000, 0, 1, 48, kFALSE),
    fit3_Pull(TString(name).Append("fit3_Pull"), TString(title).Append("fit3_Pull"), 100, -5, 5, 24, 0, 24, kFALSE),
    fit4_im(TString(name).Append("fit4_im"), TString(title).Append("fit4_im"), 2000, 0, 2000, 48, kFALSE),
    fit4_cs(TString(name).Append("fit4_cs"), TString(title).Append("fit4_cs"), 1000, 0, 1000, 48, kFALSE),
    fit4_cl(TString(name).Append("fit4_cl"), TString(title).Append("fit4_cl"), 1000, 0, 1, 48, kFALSE),
    fit4_Pull(TString(name).Append("fit4_Pull"), TString(title).Append("fit4_Pull"), 100, -5, 5, 24, 0, 24, kFALSE),
    fit3Beam_im(TString(name).Append("fit3Beam_im"), TString(title).Append("fit3Beam_im"), 2000, 0, 2000, 48, kFALSE),
    fit3Beam_cs(TString(name).Append("fit3Beam_cs"), TString(title).Append("fit3Beam_cs"), 1000, 0, 1000, 48, kFALSE),
    fit3Beam_cl(TString(name).Append("fit3Beam_cl"), TString(title).Append("fit3Beam_cl"), 1000, 0, 1, 48, kFALSE),
    fit3Beam_Pull(TString(name).Append("fit3Beam_Pull"), TString(title).Append("fit3Beam_Pull"), 100, -5, 5, 28, 0, 28, kFALSE),
    fit4Beam_im(TString(name).Append("fit4Beam_im"), TString(title).Append("fit4Beam_im"), 2000, 0, 2000, 48, kFALSE),
    fit4Beam_cs(TString(name).Append("fit4Beam_cs"), TString(title).Append("fit4Beam_cs"), 1000, 0, 1000, 48, kFALSE),
    fit4Beam_cl(TString(name).Append("fit4Beam_cl"), TString(title).Append("fit4Beam_cl"), 1000, 0, 1, 48, kFALSE),
    fit4Beam_Pull(TString(name).Append("fit4Beam_Pull"), TString(title).Append("fit4Beam_Pull"), 100, -5, 5, 28, 0, 28, kFALSE)
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
    hist_raw.CalcResult();
    hist_SubImCut.CalcResult();
    fit3_im.CalcResult();
    fit3_cs.CalcResult();
    fit3_cl.CalcResult();
    fit3_Pull.CalcResult();
    fit4_im.CalcResult();
    fit4_cs.CalcResult();
    fit4_cl.CalcResult();
    fit4_Pull.CalcResult();
    fit3Beam_im.CalcResult();
    fit3Beam_cs.CalcResult();
    fit3Beam_cl.CalcResult();
    fit3Beam_Pull.CalcResult();
    fit4Beam_im.CalcResult();
    fit4Beam_cs.CalcResult();
    fit4Beam_cl.CalcResult();
    fit4Beam_Pull.CalcResult();
}

void    GAnalysis3Mesons::Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        if(CreateHistogramsForTaggerBinning==kTRUE)
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        else
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));
    }

    if((sub_im_0>500 && sub_im_0<580) &&
        (sub_im_1>100 && sub_im_1<170) &&
        (sub_im_2>100 && sub_im_2<170))
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();if(CreateHistogramsForTaggerBinning==kTRUE)
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            else
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

            if(mm>850 && mm<1025)
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                fit3.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5));

                if(fit3.Solve()==kTRUE)
                {
                    if(fit3.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit3_im.Fill(fit3.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit3_im.Fill(fit3.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit3_cs.Fill(fit3.GetChi2(), tagger.GetTagged_t(i));
                    fit3_cl.Fill(fit3.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<24; i++)
                        fit3_Pull.Fill(fit3.GetPull(i), i);
                }

                fit4.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit3Beam.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit4Beam.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));

                if(fit4.Solve()>0)
                {
                    if(fit4.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit4_im.Fill(fit4.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit4_im.Fill(fit4.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit4_cs.Fill(fit4.GetChi2(), tagger.GetTagged_t(i));
                    fit4_cl.Fill(fit4.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<24; i++)
                        fit4_Pull.Fill(fit4.GetPull(i), i);
                }
                if(fit3Beam.Solve()>0)
                {
                    if(fit3Beam.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit3Beam_im.Fill(fit3Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit3Beam_im.Fill(fit3Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit3Beam_cs.Fill(fit3Beam.GetChi2(), tagger.GetTagged_t(i));
                    fit3Beam_cl.Fill(fit3Beam.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<28; i++)
                        fit3Beam_Pull.Fill(fit3Beam.GetPull(i), i);
                }
                if(fit4Beam.Solve()>0)
                {
                    if(fit4Beam.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit4Beam_im.Fill(fit4Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit4Beam_im.Fill(fit4Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit4Beam_cs.Fill(fit4Beam.GetChi2(), tagger.GetTagged_t(i));
                    fit4Beam_cl.Fill(fit4Beam.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<28; i++)
                        fit4Beam_Pull.Fill(fit4Beam.GetPull(i), i);
                }

            }
        }
    }
}

void    GAnalysis3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* h  = arr->GetDirectory("WithoutProton");
    GHistWriteList* folder  = h->GetDirectory("Raw");
    hist_raw.PrepareWriteList(folder, TString(name).Append("_Raw").Data());
    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
    folder  = h->GetDirectory("fit3");
    fit3_im.PrepareWriteList(folder, TString(name).Append("_fit3_im").Data());
    fit3_cs.PrepareWriteList(folder, TString(name).Append("_fit3_cs").Data());
    fit3_cl.PrepareWriteList(folder, TString(name).Append("_fit3_cl").Data());
    fit3_Pull.PrepareWriteList(folder, TString(name).Append("_fit3_Pull").Data());
    folder  = h->GetDirectory("fit4");
    fit4_im.PrepareWriteList(folder, TString(name).Append("_fit4_im").Data());
    fit4_cs.PrepareWriteList(folder, TString(name).Append("_fit4_cs").Data());
    fit4_cl.PrepareWriteList(folder, TString(name).Append("_fit4_cl").Data());
    fit4_Pull.PrepareWriteList(folder, TString(name).Append("_fit4_Pull").Data());
    folder  = h->GetDirectory("fit3Beam");
    fit3Beam_im.PrepareWriteList(folder, TString(name).Append("_fit3Beam_im").Data());
    fit3Beam_cs.PrepareWriteList(folder, TString(name).Append("_fit3Beam_cs").Data());
    fit3Beam_cl.PrepareWriteList(folder, TString(name).Append("_fit3Beam_cl").Data());
    fit3Beam_Pull.PrepareWriteList(folder, TString(name).Append("_fit3Beam_Pull").Data());
    folder  = h->GetDirectory("fit4Beam");
    fit4Beam_im.PrepareWriteList(folder, TString(name).Append("_fit4Beam_im").Data());
    fit4Beam_cs.PrepareWriteList(folder, TString(name).Append("_fit4Beam_cs").Data());
    fit4Beam_cl.PrepareWriteList(folder, TString(name).Append("_fit4Beam_cl").Data());
    fit4Beam_Pull.PrepareWriteList(folder, TString(name).Append("_fit4Beam_Pull").Data());
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    fit3_im.Reset(option);
    fit3_cs.Reset(option);
    fit3_cl.Reset(option);
    fit3_Pull.Reset(option);
    fit4_im.Reset(option);
    fit4_cs.Reset(option);
    fit4_cl.Reset(option);
    fit4_Pull.Reset(option);
    fit3Beam_im.Reset(option);
    fit3Beam_cs.Reset(option);
    fit3Beam_cl.Reset(option);
    fit3Beam_Pull.Reset(option);
    fit4Beam_im.Reset(option);
    fit4Beam_cs.Reset(option);
    fit4Beam_cl.Reset(option);
    fit4Beam_Pull.Reset(option);
}

void    GAnalysis3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
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
    isEtap(IsEtap),
    hist_raw(TString(name).Append("raw"), TString(title).Append("raw"), kFALSE),
    hist_SubImCut(TString(name).Append("SubImCut"), TString(title).Append("SubImCut"), kFALSE),
    hist_MMCut(TString(name).Append("MMCut"), TString(title).Append("MMCut"), kFALSE),
    fit3(kTRUE),
    fit4(kTRUE),
    fit3Beam(kTRUE),
    fit4Beam(kTRUE),
    fit4Proton(kTRUE),
    fit3BeamProton(kTRUE),
    fit4BeamProton(kTRUE),
    fit3_im(TString(name).Append("fit3_im"), TString(title).Append("fit3_im"), 2000, 0, 2000, 48, kFALSE),
    fit3_cs(TString(name).Append("fit3_cs"), TString(title).Append("fit3_cs"), 1000, 0, 1000, 48, kFALSE),
    fit3_cl(TString(name).Append("fit3_cl"), TString(title).Append("fit3_cl"), 1000, 0, 1, 48, kFALSE),
    fit3_Pull(TString(name).Append("fit3_Pull"), TString(title).Append("fit3_Pull"), 100, -5, 5, 24, 0, 24, kFALSE),
    fit4_im(TString(name).Append("fit4_im"), TString(title).Append("fit4_im"), 2000, 0, 2000, 48, kFALSE),
    fit4_cs(TString(name).Append("fit4_cs"), TString(title).Append("fit4_cs"), 1000, 0, 1000, 48, kFALSE),
    fit4_cl(TString(name).Append("fit4_cl"), TString(title).Append("fit4_cl"), 1000, 0, 1, 48, kFALSE),
    fit4_Pull(TString(name).Append("fit3_Pull"), TString(title).Append("fit3_Pull"), 100, -5, 5, 24, 0, 24, kFALSE),
    fit3Beam_im(TString(name).Append("fit3Beam_im"), TString(title).Append("fit3Beam_im"), 2000, 0, 2000, 48, kFALSE),
    fit3Beam_cs(TString(name).Append("fit3Beam_cs"), TString(title).Append("fit3Beam_cs"), 1000, 0, 1000, 48, kFALSE),
    fit3Beam_cl(TString(name).Append("fit3Beam_cl"), TString(title).Append("fit3Beam_cl"), 1000, 0, 1, 48, kFALSE),
    fit3Beam_Pull(TString(name).Append("fit3Beam_Pull"), TString(title).Append("fit3Beam_Pull"), 100, -5, 5, 28, 0, 28, kFALSE),
    fit4Beam_im(TString(name).Append("fit4Beam_im"), TString(title).Append("fit4Beam_im"), 2000, 0, 2000, 48, kFALSE),
    fit4Beam_cs(TString(name).Append("fit4Beam_cs"), TString(title).Append("fit4Beam_cs"), 1000, 0, 1000, 48, kFALSE),
    fit4Beam_cl(TString(name).Append("fit4Beam_cl"), TString(title).Append("fit4Beam_cl"), 1000, 0, 1, 48, kFALSE),
    fit4Beam_Pull(TString(name).Append("fit4Beam_Pull"), TString(title).Append("fit4Beam_Pull"), 100, -5, 5, 28, 0, 28, kFALSE),
    fit4Proton_im(TString(name).Append("fit4Proton_im"), TString(title).Append("fit4Proton_im"), 2000, 0, 2000, 48, kFALSE),
    fit4Proton_cs(TString(name).Append("fit4Proton_cs"), TString(title).Append("fit4Proton_cs"), 1000, 0, 1000, 48, kFALSE),
    fit4Proton_cl(TString(name).Append("fit4Proton_cl"), TString(title).Append("fit4Proton_cl"), 1000, 0, 1, 48, kFALSE),
    fit4Proton_Pull(TString(name).Append("fit4Proton_Pull"), TString(title).Append("fit4Proton_Pull"), 100, -5, 5, 28, 0, 28, kFALSE),
    fit3BeamProton_im(TString(name).Append("fit3BeamProton_im"), TString(title).Append("fit3BeamProton_im"), 2000, 0, 2000, 48, kFALSE),
    fit3BeamProton_cs(TString(name).Append("fit3BeamProton_cs"), TString(title).Append("fit3BeamProton_cs"), 1000, 0, 1000, 48, kFALSE),
    fit3BeamProton_cl(TString(name).Append("fit3BeamProton_cl"), TString(title).Append("fit3BeamProton_cl"), 1000, 0, 1, 48, kFALSE),
    fit3BeamProton_Pull(TString(name).Append("fit3BeamProton_Pull"), TString(title).Append("fit3BeamProton_Pull"), 100, -5, 5, 32, 0, 32, kFALSE),
    fit4BeamProton_im(TString(name).Append("fit4BeamProton_im"), TString(title).Append("fit4BeamProton_im"), 2000, 0, 2000, 48, kFALSE),
    fit4BeamProton_cs(TString(name).Append("fit4BeamProton_cs"), TString(title).Append("fit4BeamProton_cs"), 1000, 0, 1000, 48, kFALSE),
    fit4BeamProton_cl(TString(name).Append("fit4BeamProton_cl"), TString(title).Append("fit4BeamProton_cl"), 1000, 0, 1, 48, kFALSE),
    fit4BeamProton_Pull(TString(name).Append("fit4BeamProton_Pull"), TString(title).Append("fit4BeamProton_Pull"), 100, -5, 5, 32, 0, 32, kFALSE)
{

}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void   GAnalysis3MesonsProton::CalcResult()
{
    hist_raw.CalcResult();
    hist_SubImCut.CalcResult();
    hist_MMCut.CalcResult();
    fit3_im.CalcResult();
    fit3_cs.CalcResult();
    fit3_cl.CalcResult();
    fit3_Pull.CalcResult();
    fit4_im.CalcResult();
    fit4_cs.CalcResult();
    fit4_cl.CalcResult();
    fit4_Pull.CalcResult();
    fit3Beam_im.CalcResult();
    fit3Beam_cs.CalcResult();
    fit3Beam_cl.CalcResult();
    fit3Beam_Pull.CalcResult();
    fit4Beam_im.CalcResult();
    fit4Beam_cs.CalcResult();
    fit4Beam_cl.CalcResult();
    fit4Beam_Pull.CalcResult();
    fit4Proton_im.CalcResult();
    fit4Proton_cs.CalcResult();
    fit4Proton_cl.CalcResult();
    fit4Proton_Pull.CalcResult();
    fit3BeamProton_im.CalcResult();
    fit3BeamProton_cs.CalcResult();
    fit3BeamProton_cl.CalcResult();
    fit3BeamProton_Pull.CalcResult();
    fit4BeamProton_im.CalcResult();
    fit4BeamProton_cs.CalcResult();
    fit4BeamProton_cl.CalcResult();
    fit4BeamProton_Pull.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        if(CreateHistogramsForTaggerBinning==kTRUE)
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        else
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));
    }

    if((sub_im_0>500 && sub_im_0<580) &&
        (sub_im_1>100 && sub_im_1<170) &&
        (sub_im_2>100 && sub_im_2<170))
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();if(CreateHistogramsForTaggerBinning==kTRUE)
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            else
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

            if(mm>850 && mm<1025)
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                fit3.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5));

                if(fit3.Solve()==kTRUE)
                {
                    if(fit3.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit3_im.Fill(fit3.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit3_im.Fill(fit3.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit3_cs.Fill(fit3.GetChi2(), tagger.GetTagged_t(i));
                    fit3_cl.Fill(fit3.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<24; i++)
                        fit3_Pull.Fill(fit3.GetPull(i), i);
                }

                fit4.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit3Beam.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit4Beam.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit4Proton.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i), proton.Particle(0));
                fit3BeamProton.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i), proton.Particle(0));
                fit4BeamProton.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i), proton.Particle(0));

                if(fit4.Solve()>0)
                {
                    if(fit4.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit4_im.Fill(fit4.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit4_im.Fill(fit4.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit4_cs.Fill(fit4.GetChi2(), tagger.GetTagged_t(i));
                    fit4_cl.Fill(fit4.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<24; i++)
                        fit4_Pull.Fill(fit4.GetPull(i), i);
                }
                if(fit3Beam.Solve()>0)
                {
                    if(fit3Beam.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit3Beam_im.Fill(fit3Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit3Beam_im.Fill(fit3Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit3Beam_cs.Fill(fit3Beam.GetChi2(), tagger.GetTagged_t(i));
                    fit3Beam_cl.Fill(fit3Beam.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<28; i++)
                        fit3Beam_Pull.Fill(fit3Beam.GetPull(i), i);
                }
                if(fit4Beam.Solve()>0)
                {
                    if(fit4Beam.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit4Beam_im.Fill(fit4Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit4Beam_im.Fill(fit4Beam.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit4Beam_cs.Fill(fit4Beam.GetChi2(), tagger.GetTagged_t(i));
                    fit4Beam_cl.Fill(fit4Beam.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<28; i++)
                        fit4Beam_Pull.Fill(fit4Beam.GetPull(i), i);
                }
                if(fit4Proton.Solve()>0)
                {
                    if(fit4Proton.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit4Proton_im.Fill(fit4Proton.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit4Proton_im.Fill(fit4Proton.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit4Proton_cs.Fill(fit4Proton.GetChi2(), tagger.GetTagged_t(i));
                    fit4Proton_cl.Fill(fit4Proton.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<28; i++)
                        fit4Proton_Pull.Fill(fit4Proton.GetPull(i), i);
                }
                if(fit3BeamProton.Solve()>0)
                {
                    if(fit3BeamProton.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit3BeamProton_im.Fill(fit3BeamProton.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit3BeamProton_im.Fill(fit3BeamProton.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit3BeamProton_cs.Fill(fit3BeamProton.GetChi2(), tagger.GetTagged_t(i));
                    fit3BeamProton_cl.Fill(fit3BeamProton.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<32; i++)
                        fit3BeamProton_Pull.Fill(fit3BeamProton.GetPull(i), i);
                }
                if(fit4BeamProton.Solve()>0)
                {
                    if(fit4BeamProton.ConfidenceLevel()<0.1)
                        return;
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        fit4BeamProton_im.Fill(fit4BeamProton.GetTotalFitParticle().M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        fit4BeamProton_im.Fill(fit4BeamProton.GetTotalFitParticle().M(), tagger.GetTagged_t(i));
                    fit4BeamProton_cs.Fill(fit4BeamProton.GetChi2(), tagger.GetTagged_t(i));
                    fit4BeamProton_cl.Fill(fit4BeamProton.ConfidenceLevel(), tagger.GetTagged_t(i));
                    for(int i=0; i<32; i++)
                        fit4BeamProton_Pull.Fill(fit4BeamProton.GetPull(i), i);
                }
            }
        }
    }
}

void    GAnalysis3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* h  = arr->GetDirectory("WithProton");
    GHistWriteList* folder  = h->GetDirectory("Raw");
    hist_raw.PrepareWriteList(folder, TString(name).Append("_Raw").Data());
    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
    folder  = h->GetDirectory("fit3");
    fit3_im.PrepareWriteList(folder, TString(name).Append("_fit3_im").Data());
    fit3_cs.PrepareWriteList(folder, TString(name).Append("_fit3_cs").Data());
    fit3_cl.PrepareWriteList(folder, TString(name).Append("_fit3_cl").Data());
    fit3_Pull.PrepareWriteList(folder, TString(name).Append("_fit3_Pull").Data());
    folder  = h->GetDirectory("fit4");
    fit4_im.PrepareWriteList(folder, TString(name).Append("_fit4_im").Data());
    fit4_cs.PrepareWriteList(folder, TString(name).Append("_fit4_cs").Data());
    fit4_cl.PrepareWriteList(folder, TString(name).Append("_fit4_cl").Data());
    fit4_Pull.PrepareWriteList(folder, TString(name).Append("_fit4_Pull").Data());
    folder  = h->GetDirectory("fit3Beam");
    fit3Beam_im.PrepareWriteList(folder, TString(name).Append("_fit3Beam_im").Data());
    fit3Beam_cs.PrepareWriteList(folder, TString(name).Append("_fit3Beam_cs").Data());
    fit3Beam_cl.PrepareWriteList(folder, TString(name).Append("_fit3Beam_cl").Data());
    fit3Beam_Pull.PrepareWriteList(folder, TString(name).Append("_fit3Beam_Pull").Data());
    folder  = h->GetDirectory("fit4Beam");
    fit4Beam_im.PrepareWriteList(folder, TString(name).Append("_fit4Beam_im").Data());
    fit4Beam_cs.PrepareWriteList(folder, TString(name).Append("_fit4Beam_cs").Data());
    fit4Beam_cl.PrepareWriteList(folder, TString(name).Append("_fit4Beam_cl").Data());
    fit4Beam_Pull.PrepareWriteList(folder, TString(name).Append("_fit4Beam_Pull").Data());
    folder  = h->GetDirectory("fit4Proton");
    fit4Proton_im.PrepareWriteList(folder, TString(name).Append("_fit4Proton_im").Data());
    fit4Proton_cs.PrepareWriteList(folder, TString(name).Append("_fit4Proton_cs").Data());
    fit4Proton_cl.PrepareWriteList(folder, TString(name).Append("_fit4Proton_cl").Data());
    fit4Proton_Pull.PrepareWriteList(folder, TString(name).Append("_fit4Proton_Pull").Data());
    folder  = h->GetDirectory("fit3Beam");
    fit3BeamProton_im.PrepareWriteList(folder, TString(name).Append("_fit3BeamProton_im").Data());
    fit3BeamProton_cs.PrepareWriteList(folder, TString(name).Append("_fit3BeamProton_cs").Data());
    fit3BeamProton_cl.PrepareWriteList(folder, TString(name).Append("_fit3BeamProton_cl").Data());
    fit3BeamProton_Pull.PrepareWriteList(folder, TString(name).Append("_fit3BeamProton_Pull").Data());
    folder  = h->GetDirectory("fit4Beam");
    fit4BeamProton_im.PrepareWriteList(folder, TString(name).Append("_fit4BeamProton_im").Data());
    fit4BeamProton_cs.PrepareWriteList(folder, TString(name).Append("_fit4BeamProton_cs").Data());
    fit4BeamProton_cl.PrepareWriteList(folder, TString(name).Append("_fit4BeamProton_cl").Data());
    fit4BeamProton_Pull.PrepareWriteList(folder, TString(name).Append("_fit4BeamProton_Pull").Data());
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    hist_MMCut.Reset(option);
    fit3_im.Reset(option);
    fit3_cs.Reset(option);
    fit3_cl.Reset(option);
    fit3_Pull.Reset(option);
    fit4_im.Reset(option);
    fit4_cs.Reset(option);
    fit4_cl.Reset(option);
    fit4_Pull.Reset(option);
    fit3Beam_im.Reset(option);
    fit3Beam_cs.Reset(option);
    fit3Beam_cl.Reset(option);
    fit3Beam_Pull.Reset(option);
    fit4Beam_im.Reset(option);
    fit4Beam_cs.Reset(option);
    fit4Beam_cl.Reset(option);
    fit4Beam_Pull.Reset(option);
    fit4Proton_im.Reset(option);
    fit4Proton_cs.Reset(option);
    fit4Proton_cl.Reset(option);
    fit4Proton_Pull.Reset(option);
    fit3BeamProton_im.Reset(option);
    fit3BeamProton_cs.Reset(option);
    fit3BeamProton_cl.Reset(option);
    fit3BeamProton_Pull.Reset(option);
    fit4BeamProton_im.Reset(option);
    fit4BeamProton_cs.Reset(option);
    fit4BeamProton_cl.Reset(option);
    fit4BeamProton_Pull.Reset(option);
}

void    GAnalysis3MesonsProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MMCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3Beam_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Beam_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Proton_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Proton_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Proton_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4Proton_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3BeamProton_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3BeamProton_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3BeamProton_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3BeamProton_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4BeamProton_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4BeamProton_cs.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4BeamProton_cl.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4BeamProton_Pull.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

