#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GKinFitter.h"
#include "GKinFitterParticle.h"


class	GFit3Constraints
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit3Constraints(const Bool_t _IsEtap);
    ~GFit3Constraints();

    Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}
    Double_t        GetChi2()                       {return fitter.GetChi2();}
    Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};




class	GFit4Constraints
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4Constraints(const Bool_t _IsEtap);
    ~GFit4Constraints();

    Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}
    Double_t        GetChi2()                       {return fitter.GetChi2();}
    Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};








class	GFit4ConstraintsBeam
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsBeam(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeam();

    Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    TLorentzVector  GetTotalFitParticle();
    Double_t        GetChi2()                       {return fitter.GetChi2();}
    Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) return kTRUE; return kFALSE;}
};












class	GFit4ConstraintsProton
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsProton(const Bool_t _IsEtap);
    ~GFit4ConstraintsProton();

    Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    TLorentzVector  GetTotalFitParticle();
    Double_t        GetChi2()                       {return fitter.GetChi2();}
    Double_t        Pull(const Int_t index)         {return fitter.Pull(index);}
    Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
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




class	GFit4ConstraintsProtonExact
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsProtonExact(const Bool_t _IsEtap);
    ~GFit4ConstraintsProtonExact();

    Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    TLorentzVector  GetTotalFitParticle();
    Double_t        GetChi2()                       {return fitter.GetChi2();}
    Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
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








class	GFit4ConstraintsBeamProton
{
private:
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsBeamProton(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeamProton();

    Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    TLorentzVector  GetTotalFitParticle();
    Double_t        GetChi2()                       {return fitter.GetChi2();}
    Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
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





#endif
