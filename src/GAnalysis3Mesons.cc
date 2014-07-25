#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const char* dirName) :
    hist_raw(name, title, dirName)
{

}

GAnalysis3Mesons::~GAnalysis3Mesons()
{

}

void    GAnalysis3Mesons::Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    hist_raw.Fill(meson, tagger, CreateHistogramsForTaggerBinning);
}












GAnalysis3MesonsProton::GAnalysis3MesonsProton(const char* name, const char* title, const char* dirName) :
    hist_meson(name, title, TString(dirName).Append("/WithoutProton")),
    check_meson_proton(TString(name).Append("_checkProton").Data(), TString(title).Append(" Check Proton").Data(), TString(dirName).Append("/WithProton").Data()),
    hist_meson_proton(TString(name).Append("_proton").Data(), TString(title).Append(" Proton").Data(), TString(dirName).Append("/WithProton").Data())
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
