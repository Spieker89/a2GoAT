#ifndef __GHistEvent_h__
#define __GHistEvent_h__

//#include <iostream>

#include <TLorentzVector.h>

#include "GH1.h"

class GTreeTagger;
class GTreeMeson;

class	GHistEvent
{
private:

protected:
    GH1     im;
    GH1     mm;

public:
    GHistEvent(const char* name, const char* title, const char* dirName);
    virtual ~GHistEvent();

    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime = 0, const Double_t taggerChannel = 0);
    virtual void    Fill(const TLorentzVector& part, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    Fill(const TLorentzVector& part, const TLorentzVector& rest, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
};




class	GHistEvent3Mesons   : public GHistEvent
{
private:
    GH1     sub0_im;
    GH1     sub1_im;
    GH1     sub2_im;

protected:

public:
    GHistEvent3Mesons(const char* name, const char* title, const char* dirName);
    virtual ~GHistEvent3Mesons();

    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime = 0, const Double_t taggerChannel = 0);
    virtual void    Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    Fill(const GTreeMeson& meson, const TLorentzVector& rest, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
};

#endif
