#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const char* dirName, const Bool_t IsEtap) :
    isEtap(IsEtap),
    hist_raw(TString(name).Append("_Raw"), TString(title).Append(" Raw Data"), TString(dirName).Append("/Raw")),
    hist_SubImCut(TString(name).Append("_SubImCut"), TString(title).Append(" Sub inv. Mass Cut"), TString(dirName).Append("/IM_Cut"))
{
    if(IsEtap)
    {
        cutSubIM[0] = 500;
        cutSubIM[1] = 590;
    }
    else
    {
        cutSubIM[0] = 110;
        cutSubIM[1] = 155;
    }
    cutSubIM[2] = 110;
    cutSubIM[3] = 155;
    cutSubIM[4] = 110;
    cutSubIM[5] = 155;
}

GAnalysis3Mesons::~GAnalysis3Mesons()
{

}

void    GAnalysis3Mesons::Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    hist_raw.Fill(meson, tagger, CreateHistogramsForTaggerBinning);

    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
       (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
       (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        hist_SubImCut.Fill(meson, tagger, CreateHistogramsForTaggerBinning);
    }
}














GAnalysis3MesonsProton::GAnalysis3MesonsProton(const char* name, const char* title, const char* dirName, const Bool_t IsEtap) :
    hist_meson(name, title, TString(dirName).Append("/WithoutProton"), IsEtap),
    check_meson_proton(TString(name).Append("_checkProton").Data(), TString(title).Append(" Check Proton").Data(), TString(dirName).Append("/WithProton").Data()),
    hist_meson_proton(TString(name).Append("_proton").Data(), TString(title).Append(" Proton").Data(), TString(dirName).Append("/WithProton").Data(), IsEtap)
{

}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    if(proton.GetNParticles()>0)
    {
        if(check_meson_proton.Check(meson, proton, tagger, CreateHistogramsForTaggerBinning) == kTRUE)
        {
            hist_meson_proton.Fill(meson, tagger, CreateHistogramsForTaggerBinning);
            return;
        }
        hist_meson.Fill(meson, tagger, CreateHistogramsForTaggerBinning);
    }
}
