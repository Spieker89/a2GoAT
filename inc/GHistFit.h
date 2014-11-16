#ifndef __GHistFit_h__
#define __GHistFit_h__


#include "GFit.h"
#include "GH1.h"


class   GTreeMeson;


class  GHistFitPullParticle    : GHistLinked
{
private:
    GH1    Pull_Px;
    GH1    Pull_Py;
    GH1    Pull_Pz;
    GH1    Pull_E;

public:
    GHistFitPullParticle(const char *name, const char *title, const Bool_t linkHistogram = kTRUE);
    ~GHistFitPullParticle()  {}

    virtual void    CalcResult()        {Pull_Px.CalcResult(); Pull_Py.CalcResult(); Pull_Pz.CalcResult(); Pull_E.CalcResult();}
    virtual Int_t   Fill(Double_t x)    {}
    virtual void    Fill(const Double_t* pulls);
    virtual void    Fill(const Double_t* pulls, const Double_t taggerTime);
    virtual void    Fill(const Double_t* pulls, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "")        {Pull_Px.Reset(option); Pull_Py.Reset(option); Pull_Pz.Reset(option); Pull_E.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


class  GHistFitPull6Photons    : GHistLinked
{
private:
    GHistFitPullParticle    p0;
    GHistFitPullParticle    p1;
    GHistFitPullParticle    p2;
    GHistFitPullParticle    p3;
    GHistFitPullParticle    p4;
    GHistFitPullParticle    p5;

public:
    GHistFitPull6Photons(const char *name, const char *title, const Bool_t linkHistogram = kTRUE);
    ~GHistFitPull6Photons()  {}

    virtual void    CalcResult()        {p0.CalcResult(); p1.CalcResult(); p2.CalcResult(); p3.CalcResult(); p4.CalcResult(); p5.CalcResult();}
    virtual Int_t   Fill(Double_t x)    {}
    virtual void    Fill(const GFitStruct& fit);
    virtual void    Fill(const GFitStruct& fit, const Double_t taggerTime);
    virtual void    Fill(const GFitStruct& fit, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "")        {p0.Reset(option); p1.Reset(option); p2.Reset(option); p3.Reset(option); p4.Reset(option); p5.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


class  GHistFitStruct    : GHistLinked
{
private:
    GH1                     im;
    GH1                     ChiSq;
    GH1                     ConfidenceLevel;

public:
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
