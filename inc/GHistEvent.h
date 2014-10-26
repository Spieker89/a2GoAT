#ifndef __GHistEvent_h__
#define __GHistEvent_h__

#include <TLorentzVector.h>

#include "GH1.h"

class GTreeTagger;
class GTreeMeson;

class	GHistEvent  : public GHistLinked
{
private:

protected:
    GH1     im;
    GH1     mm;

public:
    GHistEvent(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistEvent();

    virtual void    CalcResult();
    virtual Int_t   Fill(Double_t x)    {}
    virtual void    Fill(const Double_t IM, const Double_t MM);
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime);
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime, const Double_t taggerChannel);
    virtual void    Fill(const TLorentzVector& part, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    Fill(const TLorentzVector& part, const TLorentzVector& rest, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};




class	GHistEvent3Mesons   : public GHistEvent
{
private:
    GH1     sub0_im;
    GH1     sub1_im;
    GH1     sub2_im;

protected:

public:
    GHistEvent3Mesons(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistEvent3Mesons();

    virtual void    CalcResult();
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM);
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime);
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime, const Double_t taggerChannel);
    virtual void    Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    Fill(const GTreeMeson& meson, const TLorentzVector& rest, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
};

#endif
