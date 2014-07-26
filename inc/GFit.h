#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GKinFitter.h"
#include "GKinFitterParticle.h"
#include "GHistEvent.h"


class   GTreeMeson;


class	GFit
{
private:
    Bool_t              isEtap;
    GKinFitter          fit3;
    GKinFitter          fit4;
    GH1                 fit3_ConfidenceLevel;
    GH1                 fit4_ConfidenceLevel;
    GH1                 fit3_ChiSq;
    GH1                 fit4_ChiSq;
    GH1                 im_fit3;
    GH1                 im_fit4;
    Double_t            cutConfidenceLevel;
    GH1                 im_fit3_cutCL;
    GH1                 im_fit4_cutCL;

    TFile*  GammaResFile;
    TH2F*   GammaEloss;
    TH2F*   GammaERes;
    TH2F*   GammaThetaRes;
    TH2F*   GammaPhiRes;

    void    Fit3(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    void    Fit4(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
protected:

public:
    GFit(const char* name, const char* title, const char* dirName, const Bool_t IsEtap);
    virtual ~GFit();

    void    Fit(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
};







#endif
