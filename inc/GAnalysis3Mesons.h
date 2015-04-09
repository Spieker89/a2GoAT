#ifndef __GAnalysis3Mesons_h__
#define __GAnalysis3Mesons_h__


#include "GHistEvent.h"
#include "GCheckProton.h"
#include "GFit.h"
//#include "GHistFit.h"
#include "APLCON.hpp"


class GTreeTagger;
class GTreeParticle;
class GTreeMeson;


class   GAnalysis3Mesons  : public GHistLinked
{
private:
    Bool_t              isEtap;

    GHistEvent3Mesons   hist_raw;
    GHistEvent3Mesons   hist_SubImCut;

    GHistEvent3Mesons   hist_MMCut;
//    GHistIterativeFit   hist_fit1;
//    GHistIterativeFit   hist_fit3;
//    GHistIterativeFit   hist_fit4;

//    GFit1Constraints        fit1;
//    GFit3Constraints        fit3;
//    GFit4Constraints        fit4;

    APLCON  fitter;

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

    std::vector<FitParticle>    aplconPhotons;


    GH1     fitResultIM;
    GH1     fitZVertex;

protected:

public:
    Double_t            cutSubIM[6];
    Double_t            cutMM[2];

    GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram = kTRUE);
    ~GAnalysis3Mesons();

    virtual void    CalcResult();
            Bool_t  IsEtap()    const   {return isEtap;}
    virtual Int_t   Fill(Double_t x)    {return 0;}
    virtual void    Fill(const GTreeMeson& meson, const GTreeParticle& photons, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}

    void    SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max);
    void    SetCutMM(const Double_t min, const Double_t max);
};


class   GAnalysis3MesonsProton  : public GHistLinked
{
private:
    Bool_t              isEtap;

    GCheckProton        checkProton;

    GHistEvent3Mesons   hist_raw;
    GHistBGSub2         hist_raw_TOF;

    GHistEvent3Mesons   hist_SubImCut;
    GHistBGSub2         hist_SubImCut_TOF;

    GHistEvent3Mesons   hist_MMCut;
    GHistBGSub2         hist_TOF;
//    GHistIterativeFit   hist_fit1;
//    GHistIterativeFit   hist_fit3;
//    GHistIterativeFit   hist_fit4;
//    GHistIterativeFit   hist_fit7Proton;

//    GFit1Constraints            fit1;
//    GFit3Constraints            fit3;
//    GFit4Constraints            fit4;
//    GFit7ConstraintsProton      fit7Proton;

protected:

public:
    GAnalysis3MesonsProton(const char* name, const char* title, const Bool_t IsEtap, Bool_t linkHistogram = kTRUE);
    ~GAnalysis3MesonsProton();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)    {return 0;}
    virtual void        Fill(const GTreeMeson& meson, const GTreeParticle& photons, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}

    void	SetHistMeson(const Double_t sub0_min, const Double_t sub0_max,
                         const Double_t sub1_min, const Double_t sub1_max,
                         const Double_t sub2_min, const Double_t sub2_max,
                         const Double_t mm_min, const Double_t mm_max);
    void    SetFitMeson(const Double_t fit3_CutConfidenceLevel, const Double_t fit4_CutConfidenceLevel);
    void	SetCheckProton(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity);
    void	SetHistMesonProton(const Double_t sub0_min, const Double_t sub0_max,
                               const Double_t sub1_min, const Double_t sub1_max,
                               const Double_t sub2_min, const Double_t sub2_max,
                               const Double_t mm_min, const Double_t mm_max);
    void    SetFitMesonProton(const Double_t fit3_CutConfidenceLevel, const Double_t fit4_CutConfidenceLevel);
};

#endif
