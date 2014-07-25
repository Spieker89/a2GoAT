#ifndef __GHistEvent_h__
#define __GHistEvent_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include <TLorentzVector.h>

#include "GH1.h"

class GTreeTagger;

class	GHistEvent
{
private:
    GH1     im;
    GH1     mm;

protected:

public:
    GHistEvent(const char* name, const char* title, const char* dirName);
    virtual ~GHistEvent();

    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime = 0, const Double_t taggerChannel = 0);
    virtual void    Fill(const TLorentzVector& part, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    Fill(const TLorentzVector& part, const TLorentzVector& rest, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
};

#endif
