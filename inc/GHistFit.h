#ifndef __GHistFit_h__
#define __GHistFit_h__


#include "GFit.h"
#include "GH1.h"


class   GTreeMeson;


class  GHistFitStruct
{
public:
    GH1    im;
    GH1    ChiSq;
    GH1    ConfidenceLevel;

    GHistFitStruct(const char *name, const char *title, const Bool_t linkHistogram);
    ~GHistFitStruct()  {}

    void    Fill(const GFitStruct& fit);
    void    Fill(const GFitStruct& fit, const Double_t taggerTime);
    void    Fill(const GFitStruct& fit, const Double_t taggerTime, const Int_t taggerChannel);
};



class	GHistFit3Constraints    : GHistLinked
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
    GHistFit3Constraints(const char *name, const char *title, const Bool_t IsEtap, const Bool_t linkHistogram);
    virtual ~GHistFit3Constraints();

            Bool_t  IsEtap()    const   {return isEtap;}
    virtual void    Fill(const GFit3Constraints& fit);
    virtual void    Fill(const GFit3Constraints& fit, const Double_t taggerTime);
    virtual void    Fill(const GFit3Constraints& fit, const Double_t taggerTime, const Int_t taggerChannel);
};

/*
class	GHistFit4Constraints    :   public  GHistFit3Constraints
{
private:

protected:

public:
    GHistFit4Constraints(const Bool_t IsEtap);
    virtual ~GHistFit4Constraints();

    virtual Bool_t  Fit(const GTreeMeson& meson, const TLorentzVector &beamAndTarget);
};
*/


#endif
