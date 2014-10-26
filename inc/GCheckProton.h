#ifndef __GCheckProton_h__
#define __GCheckProton_h__


#include "GH1.h"


class GTreeTagger;
class GTreeParticle;
class GTreeMeson;

class   GCheckProtonHist
{
private:
    GH1         protonAngeDiff;
    GH1         protonAngeDiffSmalest;
    GH1         protonCoplanarity;

public:
    GCheckProtonHist(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    ~GCheckProtonHist();

            void    CalcResult()                                                                                                                            {protonAngeDiff.CalcResult(); protonAngeDiffSmalest.CalcResult(); protonCoplanarity.CalcResult();}
            void    Fill(const Double_t _ProtonAngeDiff)                                            {protonAngeDiff.Fill(_ProtonAngeDiff);}
            void    Fill(const Double_t _ProtonAngeDiffSmalest, const Double_t _ProtonCoplanarity)  {protonAngeDiffSmalest.Fill(_ProtonAngeDiffSmalest); protonCoplanarity.Fill(_ProtonCoplanarity);}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "")                                                                                                            {protonAngeDiff.Reset(option); protonAngeDiffSmalest.Reset(option); protonCoplanarity.Reset(option);}
            //void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE)                       {protonAngeDiff.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); protonAngeDiffSmalest.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); protonCoplanarity.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);}
};


class   GCheckProton : public GHistLinked
{
private:
    GCheckProtonHist    raw;
    GCheckProtonHist    cutCoplanarity;
    GCheckProtonHist    cutProtonAngle;
    GCheckProtonHist    cutBoth;

    Double_t    CutProtonAngleDiff;
    Double_t    CutProtonCoplanarity[2];

protected:

public:
    GCheckProton(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    ~GCheckProton();

    virtual void    CalcResult();
            Bool_t  Check(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger);
    virtual Int_t   Fill(Double_t x)    {}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
            void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
            void    SetCuts(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};





#endif
