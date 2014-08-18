#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GKinFitter.h"
#include "GKinFitterParticle.h"
#include "GHistEvent.h"


class   GTreeMeson;

class  GFitPulls4Vector
{
private:
    GH1 Pull_Px;
    GH1 Pull_Py;
    GH1 Pull_Pz;
    GH1 Pull_E;

public:
    GFitPulls4Vector(const char* name, const char* title, const char* dirName);
    virtual ~GFitPulls4Vector();

    void    Fill(const Double_t px, const Double_t py, const Double_t pz, const Double_t e)  {Pull_Px.Fill(px); Pull_Py.Fill(py); Pull_Pz.Fill(pz); Pull_E.Fill(e);}
};

class  GFitPulls6Photons
{
private:
    GFitPulls4Vector g0;
    GFitPulls4Vector g1;
    GFitPulls4Vector g2;
    GFitPulls4Vector g3;
    GFitPulls4Vector g4;
    GFitPulls4Vector g5;

public:
    GFitPulls6Photons(const char* name, const char* title, const char* dirName);
    virtual ~GFitPulls6Photons();

    void    Fill(GKinFitter& fitter)
    {
        g0.Fill(fitter.Pull(0), fitter.Pull(1), fitter.Pull(2), fitter.Pull(3));
        g1.Fill(fitter.Pull(4), fitter.Pull(5), fitter.Pull(6), fitter.Pull(7));
        g2.Fill(fitter.Pull(8), fitter.Pull(9), fitter.Pull(10), fitter.Pull(11));
        g3.Fill(fitter.Pull(12), fitter.Pull(13), fitter.Pull(14), fitter.Pull(15));
        g4.Fill(fitter.Pull(16), fitter.Pull(17), fitter.Pull(18), fitter.Pull(19));
        g5.Fill(fitter.Pull(20), fitter.Pull(21), fitter.Pull(22), fitter.Pull(23));
    }
};


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
    GFitPulls6Photons   fit3_Pulls;
    GFitPulls6Photons   fit4_Pulls;
    GH1                 im_fit3;
    GH1                 im_fit4;
    Double_t            cutConfidenceLevel;
    GH1                 im_fit3_cutCL;
    GH1                 im_fit4_cutCL;
    GFitPulls6Photons   fit3_Pulls_cutCL;
    GFitPulls6Photons   fit4_Pulls_cutCL;

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
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
};







#endif
