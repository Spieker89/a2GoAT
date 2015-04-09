#ifndef __GFit_h__
#define __GFit_h__


#include <TLorentzVector.h>

#include "GH1.h"
#include "APLCON.hpp"



class	GFit    : public GHistLinked
{
protected:
    // lightweight structure for linking to fitter
    struct FitParticle {
        void SetFromVector(const TLorentzVector& p_) {
            Ek = p_.E()-p_.M();
            Theta = p_.Theta();
            Phi = p_.Phi();
            Ek_Sigma = 0.02*Ek*pow(Ek,-0.36);
            Theta_Sigma = 2.5*TMath::DegToRad();
            if(Theta>20*TMath::DegToRad() && Theta<160*TMath::DegToRad()) {
                Phi_Sigma = Theta_Sigma/sin(Theta);
            }
            else {
                Phi_Sigma = 1*TMath::DegToRad();
            }
        }
        static TLorentzVector Make(const std::vector<double>& EkThetaPhi,
                                   const Double_t m){
            const double E = EkThetaPhi[0] + m;
            const Double_t p = sqrt( E*E - m*m );
            TVector3 pv(1,0,0);
            pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
            TLorentzVector l(pv, E);
            return l;
        }
        static TLorentzVector Make(const FitParticle& p,
                                   const Double_t m) {
            return Make(std::vector<double>{p.Ek, p.Theta, p.Phi}, m);
        }
        std::vector<double*> Link() {
            return {std::addressof(Ek),
                        std::addressof(Theta),
                        std::addressof(Phi)};
        }
        std::vector<double*> LinkSigma() {
            return {std::addressof(Ek_Sigma),
                        std::addressof(Theta_Sigma),
                        std::addressof(Phi_Sigma)};
        }


        double Ek;
        double Ek_Sigma;
        double Theta;
        double Theta_Sigma;
        double Phi;
        double Phi_Sigma;

    };

private:
    std::vector<FitParticle>    aplconPhotons;

    GH1         im;
    GH1         zVertex;
    GH1         theta;
    GH1         phi;
    GH1         chiSq;
    GH1         confidenceLevel;

protected:
    std::string name;
    APLCON      fitter;
    APLCON::Result_t result;

    virtual void    AddConstraints() = 0;
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5)
    {
        aplconPhotons[0].SetFromVector(p0);
        aplconPhotons[1].SetFromVector(p1);
        aplconPhotons[2].SetFromVector(p2);
        aplconPhotons[3].SetFromVector(p3);
        aplconPhotons[4].SetFromVector(p4);
        aplconPhotons[5].SetFromVector(p5);
    }


    GFit(const char* name);

public:
    ~GFit()                                 {}

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)    {return 0;}
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* _Name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}
    virtual void Solve(const double time, const int channel);
};










class	GFit1Constraints    : public GFit
{
private:
    TLorentzVector  beamAndTarget;

    GHistBGSub2     pulls;

protected:
    GFit1Constraints(const char* name)  :
        GFit(name),
        pulls(TString(name).Append("_pulls"), TString(name).Append("_pulls"), 100, -5, 5, 18, 0, 18, kFALSE)
    {AddConstraints();}

public:
    GFit1Constraints()                  :
        GFit("fit1"),
        pulls(TString("fit1").Append("_pulls"), TString("fit1").Append("_pulls"), 100, -5, 5, 18, 0, 18, kFALSE)
    {AddConstraints();}
    ~GFit1Constraints() {}

    virtual void        AddConstraints();
    virtual void        CalcResult();
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* _Name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        Solve(const double time, const int channel);
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& _BeamAndTarget)
    {
        beamAndTarget   = _BeamAndTarget;
        GFit::Set(p0, p1, p2, p3, p4, p5);
    }
};









class	GFit3Constraints    : public GFit
{
protected:
    GFit3Constraints(const char* name)  :   GFit(name)      {AddConstraints();}

public:
    GFit3Constraints()                  :   GFit("fit3")    {AddConstraints();}
    ~GFit3Constraints()                                     {}

    virtual void    AddConstraints();
    /*virtual TLorentzVector  GetMeson()                      {return fitter.GetTotalFitParticle();}
    virtual TLorentzVector  GetTotal()                      {return GetMeson();}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}*/
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5)
    {
        GFit::Set(p0, p1, p2, p3, p4, p5);
    }
};









class	GFit4Constraints    : public GFit1Constraints
{
public:
    GFit4Constraints()                  :   GFit1Constraints("fit4")    {AddConstraints();}
    ~GFit4Constraints()                                                 {}

    virtual void    AddConstraints();
    /*virtual TLorentzVector  GetMeson()                      {return fitter.GetTotalFitParticle();}
    virtual TLorentzVector  GetTotal()                      {return GetMeson();}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}*/
};











class	GFitProton    : public GFit
{
private:
    FitParticle     proton;
    TLorentzVector  beamAndTarget;

public:
    GFitProton();
    ~GFitProton()                   {}

    virtual void    AddConstraints()    {}
    /*virtual TLorentzVector  GetMeson();
    virtual TLorentzVector  GetTotal()                      {return GetMeson() + fitter.GetParticle(6);}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return fitter.GetParticle(6);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}*/
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& _Proton,
                const TLorentzVector& _BeamAndTarget)
    {
        beamAndTarget   = _BeamAndTarget;
        proton.SetFromVector(_Proton);
        GFit::Set(p0, p1, p2, p3, p4, p5);
    }
};










/*
class   GHistFit    : public    GHistLinked
{
private:
    Int_t       nPulls;
    GH1         im;
    GH1         sub0im;
    GH1         sub1im;
    GH1         sub2im;
    GH1         theta;
    GH1         phi;
    GH1         Pim;
    GH1         Ptheta;
    GH1         Pphi;
    GH1         chiSq;
    GH1         VchiSq;
    GH1         CchiSq;
    GH1         confidenceLevel;
    GH1         VconfidenceLevel;
    GH1         CconfidenceLevel;
    GHistBGSub2 pulls;
public:
    GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram= kTRUE);
    ~GHistFit();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {return 0;}
    virtual Int_t       Fill(GFit& fitter, const Double_t taggerTime);
    virtual Int_t       Fill(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}
};


class   GHistIterativeFit    : public    GHistLinked
{
private:
    GHistBGSub2 im;
    GHistBGSub2 sub0im;
    GHistBGSub2 sub1im;
    GHistBGSub2 sub2im;
    GHistBGSub2 mm;
    GHistBGSub2 totE;
    GHistBGSub2 totPx;
    GHistBGSub2 totPy;
    GHistBGSub2 totPz;
    GHistBGSub2 VchiSq;
    GHistBGSub2 VconfidenceLevel;
    GHistBGSub2 CchiSq;
    GHistBGSub2 CconfidenceLevel;
    GHistBGSub2 chiSq;
    GHistBGSub2 confidenceLevel;
    GHistFit    final;

public:
    GHistIterativeFit(const char* name, const char* title, const Int_t _NPulls, const Int_t _NSteps, Bool_t linkHistogram= kTRUE);
    ~GHistIterativeFit();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {return 0;}
    virtual Int_t       Fill(GFit& fitter);
    virtual Int_t       FillFinal(GFit& fitter, const Double_t taggerTime)                              {return final.Fill(fitter, taggerTime);}
    virtual Int_t       FillFinal(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel)   {return final.Fill(fitter, taggerTime, taggerChannel);}
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}
};
*/

#endif
