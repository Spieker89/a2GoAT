#ifndef __GAnalysis3Mesons_h__
#define __GAnalysis3Mesons_h__


#include "GHistEvent.h"
#include "GCheckProton.h"
#include "GFit.h"


class GTreeTagger;
class GTreeParticle;
class GTreeMeson;


class   GAnalysis3Mesons
{
private:
    Bool_t              isEtap;
    GHistEvent3Mesons   hist_raw;

    Double_t            cutSubIM[6];
    GHistEvent3Mesons   hist_SubImCut;

    Double_t            cutMM[2];
    GHistEvent3Mesons   hist_MmCut;

protected:

public:
    GAnalysis3Mesons(const char* name, const char* title, const char* dirName, const Bool_t IsEtap);
    ~GAnalysis3Mesons();

    virtual Bool_t    Fill(const GTreeMeson& meson, const GTreeTagger &tagger, const Bool_t CreateHistogramsForTaggerBinning);
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);

    void    SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max);
    void    SetCutMM(const Double_t min, const Double_t max);
};


class   GAnalysis3MesonsProton
{
private:
    GAnalysis3Mesons   hist_meson;
    GFit               fit_meson;
    GCheckProton       check_meson_proton;
    GAnalysis3Mesons   hist_meson_proton;
    GFit               fit_meson_proton;

protected:

public:
    GAnalysis3MesonsProton(const char* name, const char* title, const char* dirName, const Bool_t IsEtap);
    ~GAnalysis3MesonsProton();

    virtual void    Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);

    void	SetHistMeson(const Double_t sub0_min, const Double_t sub0_max,
                         const Double_t sub1_min, const Double_t sub1_max,
                         const Double_t sub2_min, const Double_t sub2_max,
                         const Double_t mm_min, const Double_t mm_max);
    void	SetCheckProton(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity);
    void	SetHistMesonProton(const Double_t sub0_min, const Double_t sub0_max,
                               const Double_t sub1_min, const Double_t sub1_max,
                               const Double_t sub2_min, const Double_t sub2_max,
                               const Double_t mm_min, const Double_t mm_max);
};

#endif
