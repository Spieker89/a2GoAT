#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GKinFitter.h"
#include "GKinFitterParticle.h"
#include "GHistEvent.h"


class   GTreeMeson;

class  GFitStruct
{
public:
    Double_t    im;
    Double_t    ChiSq;
    Double_t    ConfidenceLevel;
    Double_t    PullPhotons[24];
    Double_t    PullBeam[4];
    Double_t    PullProton[4];

    GFitStruct();
    //GFitStructEntry(const GFitStructEntry& c);
    ~GFitStruct()  {}

    //GFitStructEntry&    operator =(const GFitStructEntry& c);
};


class	GFit3Constraints
{
private:
    Bool_t      isEtap;

    GKinFitterParticle  photons[6];

    static  TFile*  GammaResFile;
    static  TH2F*   GammaEloss;
    static  TH2F*   GammaERes;
    static  TH2F*   GammaThetaRes;
    static  TH2F*   GammaPhiRes;

protected:
    GKinFitter  fitter;

    GFitStruct  result;


    GFit3Constraints(const Int_t npart, const Int_t ncon, const Bool_t IsEtap);

    virtual Bool_t InitFit(const GTreeMeson& meson);
            Bool_t SetPhotons(const GTreeMeson& meson);

public:
    GFit3Constraints(const Bool_t IsEtap);
    virtual ~GFit3Constraints();

    const   GFitStruct& GetResult() const   {return result;}
            Bool_t      IsEtap()    const   {return isEtap;}
    virtual Bool_t      Fit(const GTreeMeson& meson);
};


class	GFit4Constraints    :   public  GFit3Constraints
{
private:

protected:
    GFit4Constraints(const Int_t npart, const Int_t ncon, const Bool_t IsEtap);

    virtual Bool_t  InitFit(const GTreeMeson& meson, const TLorentzVector &beamAndTarget);

public:
    GFit4Constraints(const Bool_t IsEtap);
    virtual ~GFit4Constraints();

    virtual Bool_t  Fit(const GTreeMeson& meson, const TLorentzVector &beamAndTarget);
};


class	GFit3ConstraintsBeam    :   public  GFit3Constraints
{
private:

protected:
    GFit3ConstraintsBeam(const Int_t npart, const Int_t ncon, const Bool_t IsEtap);

    virtual Bool_t  InitFit(const GTreeMeson& meson, const TLorentzVector &beamAndTarget);

public:
    GFit3ConstraintsBeam(const Bool_t IsEtap);
    virtual ~GFit3ConstraintsBeam();

    virtual Bool_t  Fit(const GTreeMeson& meson, const TLorentzVector &beamAndTarget);
};



class	GFit4ConstraintsBeam    :   public  GFit4Constraints
{
private:

protected:
    GFit4ConstraintsBeam(const Int_t npart, const Int_t ncon, const Bool_t IsEtap);

    virtual Bool_t  InitFit(const GTreeMeson& meson, const TLorentzVector &beamAndTarget);

public:
    GFit4ConstraintsBeam(const Bool_t IsEtap);
    virtual ~GFit4ConstraintsBeam();

    virtual Bool_t  Fit(const GTreeMeson& meson, const TLorentzVector &beamAndTarget);
};



#endif
