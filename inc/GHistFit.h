#ifndef __GHistFit_h__
#define __GHistFit_h__


#include "GFit.h"
#include "GH1.h"


class   GTreeMeson;


class  GHistFitStruct    : GHistLinked
{
public:
    GH1    im;
    GH1    ChiSq;
    GH1    ConfidenceLevel;

    GHistFitStruct(const char *name, const char *title, const Bool_t linkHistogram = kTRUE);
    ~GHistFitStruct()  {}

    virtual void    CalcResult()        {im.CalcResult(); ChiSq.CalcResult(); ConfidenceLevel.CalcResult();}
    virtual Int_t   Fill(Double_t x)    {}
    virtual void    Fill(const GFitStruct& fit);
    virtual void    Fill(const GFitStruct& fit, const Double_t taggerTime);
    virtual void    Fill(const GFitStruct& fit, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "")        {im.Reset(option); ChiSq.Reset(option); ConfidenceLevel.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};



class	GHistFit    : GHistLinked
{
private:
    Bool_t      isEtap;

    GHistFitStruct    raw;
    GHistFitStruct    CutChiSq;
    GHistFitStruct    CutConfidenceLevel;
    GHistFitStruct    CutBoth;

    Double_t    ConfidenceLevelCut;
    Double_t    ChiSqCut;

protected:

public:
    GHistFit(const char *name, const char *title, const Bool_t IsEtap, const Bool_t linkHistogram = kTRUE);
    virtual ~GHistFit();

    virtual void        CalcResult()        {raw.CalcResult(); CutChiSq.CalcResult(); CutConfidenceLevel.CalcResult(); CutBoth.CalcResult();}
    virtual Int_t       Fill(Double_t x)    {}
    virtual void        Fill(const GFit3Constraints& fit);
    virtual void        Fill(const GFit3Constraints& fit, const Double_t taggerTime);
    virtual void        Fill(const GFit3Constraints& fit, const Double_t taggerTime, const Int_t taggerChannel);
            Bool_t      IsEtap()    const   {return isEtap;}
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "")        {raw.Reset(option); CutChiSq.Reset(option); CutConfidenceLevel.Reset(option); CutBoth.Reset(option);}
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


#endif
