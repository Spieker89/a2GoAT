#ifndef __GAnalysis3Mesons_h__
#define __GAnalysis3Mesons_h__


#include "GHistEvent.h"
#include "GCheckProton.h"


class GTreeTagger;
class GTreeParticle;
class GTreeMeson;


class   GAnalysis3Mesons
{
private:
    GHistEvent3Mesons   hist_raw;

protected:

public:
    GAnalysis3Mesons(const char* name, const char* title, const char* dirName);
    ~GAnalysis3Mesons();

    void    Fill(const GTreeMeson& meson, const GTreeTagger &tagger, const Bool_t CreateHistogramsForTaggerBinning);
};


class   GAnalysis3MesonsProton
{
private:
    GAnalysis3Mesons   hist_meson;
    GCheckProton       check_meson_proton;
    GAnalysis3Mesons   hist_meson_proton;

protected:

public:
    GAnalysis3MesonsProton(const char* name, const char* title, const char* dirName);
    ~GAnalysis3MesonsProton();

    void    Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
};

#endif
