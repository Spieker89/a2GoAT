#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GKinFitter.h"
#include "GKinFitterParticle.h"
#include "GHistEvent.h"


class   GTreeMeson;
/*
class  GFitPulls4Vector
{
private:
    Double_t    Pull_Px;
    Double_t    Pull_Py;
    Double_t    Pull_Pz;
    Double_t    Pull_E;

public:
    GFitPulls4Vector(const Double_t px, const Double_t py, const Double_t pz, const Double_t e)   : Pull_Px(px), Pull_Py(py), Pull_Pz(pz), Pull_E(e)  {}
    ~GFitPulls4Vector();

    void    Fill(const Double_t px, const Double_t py, const Double_t pz, const Double_t e)  {Pull_Px = px; Pull_Py = py); Pull_Pz(pz); Pull_E(e);}
};

class  GFitPulls6Photons
{
private:
    GFitPulls4Vector g0;
    GFitPulls4Vector g1;
    GFitPulls4Vector g2;
    GFitPulls4Vector g3;
    GFitPulls4Vector g4;
    GFitPulls4Vector g5;

public:
    GFitPulls6Photons(const char* name, const char* title);
    virtual ~GFitPulls6Photons();

    virtual void    CalcResult()
    {
        g0.CalcResult();
        g1.CalcResult();
        g2.CalcResult();
        g3.CalcResult();
        g4.CalcResult();
        g5.CalcResult();
    }
    void    Fill(GKinFitter& fitter)
    {
        g0.Fill(fitter.Pull(0), fitter.Pull(1), fitter.Pull(2), fitter.Pull(3));
        g1.Fill(fitter.Pull(4), fitter.Pull(5), fitter.Pull(6), fitter.Pull(7));
        g2.Fill(fitter.Pull(8), fitter.Pull(9), fitter.Pull(10), fitter.Pull(11));
        g3.Fill(fitter.Pull(12), fitter.Pull(13), fitter.Pull(14), fitter.Pull(15));
        g4.Fill(fitter.Pull(16), fitter.Pull(17), fitter.Pull(18), fitter.Pull(19));
        g5.Fill(fitter.Pull(20), fitter.Pull(21), fitter.Pull(22), fitter.Pull(23));
    }
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0)
    {
        if(!arr)
            return;

        g0.PrepareWriteList(arr, TString(name).Append("_g0"));
        g1.PrepareWriteList(arr, TString(name).Append("_g1"));
        g2.PrepareWriteList(arr, TString(name).Append("_g2"));
        g3.PrepareWriteList(arr, TString(name).Append("_g3"));
        g4.PrepareWriteList(arr, TString(name).Append("_g4"));
        g5.PrepareWriteList(arr, TString(name).Append("_g5"));
    }
    virtual void    Reset(Option_t* option = "")
    {
        g0.Reset(option);
        g1.Reset(option);
        g2.Reset(option);
        g3.Reset(option);
        g4.Reset(option);
        g5.Reset(option);
    }
};
*/

class  GFitStructEntry
{
public:
    Double_t    im;
    Double_t    ChiSq;
    Double_t    ConfidenceLevel;
    Double_t    PullPhotons[24];
    Double_t    PullBeam[4];
    Double_t    PullProton[4];

    GFitStructEntry();
    GFitStructEntry(const GFitStructEntry& c);
    ~GFitStructEntry()  {}

    GFitStructEntry&    operator =(const GFitStructEntry& c);
};

struct  GFitStruct
{
    GFitStructEntry    raw;
    GFitStructEntry    CutChiSq;
    GFitStructEntry    CutConfidenceLevel;
    GFitStructEntry    CutBoth;
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
    Double_t    cutConfidenceLevel;
    Double_t    cutChiSq;

    GFitStruct  result;


    GFit3Constraints(const Int_t npart, const Int_t ncon, const Bool_t IsEtap);

    virtual void    InitFit(const GTreeMeson& meson);
            void    SetPhotons(const GTreeMeson& meson);

public:
    GFit3Constraints(const Bool_t IsEtap);
    virtual ~GFit3Constraints();

            Bool_t  IsEtap()    const   {return isEtap;}
    virtual void    Fit(const GTreeMeson& meson);
};


class	GFit4Constraints    :   public  GFit3Constraints
{
private:

protected:
    //GFit3Constraints(const Int_t npart, const Int_t ncon, const Bool_t IsEtap);

    virtual void    InitFit(const GTreeMeson& meson);
            void    SetPhotons(const GTreeMeson& meson);

public:
    GFit4Constraints(const Bool_t IsEtap);
    virtual ~GFit4Constraints();

    virtual void    Fit(const GTreeMeson& meson);
};



#endif
