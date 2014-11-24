#ifndef __GCheckProton_h__
#define __GCheckProton_h__


#include <TLorentzVector.h>

#include "GH1.h"


class GTreeParticle;

class   GCheckProtonHist
{
private:
    GH1         protonAngeDiff;
    GH1         protonCoplanarity;

public:
    GCheckProtonHist(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    ~GCheckProtonHist();

            void    CalcResult()                                                            {protonAngeDiff.CalcResult(); protonCoplanarity.CalcResult();}
            void    Fill(const Double_t _ProtonAngeDiff, const Double_t _ProtonCoplanarity) {protonAngeDiff.Fill(_ProtonAngeDiff); protonCoplanarity.Fill(_ProtonCoplanarity);}
            void    Fill(const Double_t _ProtonAngeDiff, const Double_t _ProtonCoplanarity, const Double_t taggerTime) {protonAngeDiff.Fill(_ProtonAngeDiff, taggerTime); protonCoplanarity.Fill(_ProtonCoplanarity, taggerTime);}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "")                                            {protonAngeDiff.Reset(option); protonCoplanarity.Reset(option);}
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
            Bool_t  Check(const GTreeParticle& meson, const GTreeParticle& proton, const TLorentzVector& beamAndTarget, const Double_t taggerTime);
    virtual Int_t   Fill(Double_t x)    {}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
            void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
            void    SetCuts(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};





#endif
