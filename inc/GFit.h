#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GHistManager.h"
#include "GH1.h"
#include "GKinFitter.h"
#include "GKinFitterParticle.h"



class	GFit
{
private:

public:
    GFit()  {}
    ~GFit() {}

    virtual Double_t        ConfidenceLevel()           = 0;
    virtual TLorentzVector  GetTotalFitParticle()       = 0;
    virtual Double_t        GetChi2()                   = 0;
    virtual Double_t        GetPull(const Int_t index)  = 0;
};


class	GFit3Constraints    : public GFit
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit3Constraints(const Bool_t _IsEtap);
    ~GFit3Constraints();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};




class	GFit4Constraints    : public GFit
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4Constraints(const Bool_t _IsEtap);
    ~GFit4Constraints();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};








class	GFit4ConstraintsBeam    : public GFit
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsBeam(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeam();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};












class	GFit4ConstraintsProton    : public GFit
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsProton(const Bool_t _IsEtap);
    ~GFit4ConstraintsProton();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget,
                const TLorentzVector& proton);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};




class	GFit4ConstraintsProtonExact    : public GFit
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsProtonExact(const Bool_t _IsEtap);
    ~GFit4ConstraintsProtonExact();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget,
                const TLorentzVector& proton);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};








class	GFit4ConstraintsBeamProton    : public GFit
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsBeamProton(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeamProton();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget,
                const TLorentzVector& proton);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};









class   GHistFit    : public    GHistLinked
{
private:
    Int_t       nPulls;
    GH1         im;
    GH1         chiSq;
    GH1         confidenceLevel;
    GHistBGSub2 pulls;
public:
    GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram= kTRUE);
    ~GHistFit();

    virtual void        CalcResult()                    {im.CalcResult(); chiSq.CalcResult(); confidenceLevel.CalcResult(); pulls.CalcResult();}
    virtual Int_t       Fill(Double_t x)                {}
    virtual Int_t       Fill(GFit& fitter, const Double_t taggerTime);
    virtual Int_t       Fill(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "")    {im.Reset(option); chiSq.Reset(option); confidenceLevel.Reset(option); pulls.Reset(option);}
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE)   {im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); pulls.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);}
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


#endif
