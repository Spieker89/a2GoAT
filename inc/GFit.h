#ifndef     ___GFIT___
#define	    ___GFIT___

#include <vector>
#include <TLorentzVector.h>
#include "/disk/nobackup/001/spieker/Mainz1/APLCON/src/APLCON.hpp"
#include "TH2.h"
#include "TFile.h"
#include "GTreeParticle.h"
using namespace std;

class	GFit
{

private:


    double Ek;
    double Ek_Sigma;
    double              Theta;
    double              Phi;
    double              Theta_Sigma;
    double              Phi_Sigma;

    double              aplconProtonTheta;
    double              aplconProtonPhi;
    double              aplconProtonThetaSigma;
    double              aplconProtonPhiSigma;

    APLCON::Variable_Settings_t Ek_Setting;
    APLCON::Variable_Settings_t Theta_Setting;
    APLCON::Variable_Settings_t Phi_Setting;


struct FitParticle {

    //get the uncertainties(cite Adlarson in case)
    TFile *photonfile = new TFile("/disk/nobackup/001/spieker/Mainz1/a2GoAT/configfiles/uncertainties/photon_uncertainties.root","OPEN");
    TFile *protonfile = new TFile("/disk/nobackup/001/spieker/Mainz1/a2GoAT/configfiles/uncertainties/proton_uncertainties.root","OPEN");
	
    TH2F *photon_CB_energy = (TH2F*)photonfile->Get("CB_E");
    TH2F *photon_CB_theta = (TH2F*)photonfile->Get("CB_theta");
    TH2F *photon_CB_phi = (TH2F*)photonfile->Get("CB_phi");
    TH2F *photon_TAPS_energy = (TH2F*)photonfile->Get("TAPS_E"); 
    TH2F *photon_TAPS_theta = (TH2F*)photonfile->Get("TAPS_theta");
    TH2F *photon_TAPS_phi = (TH2F*)photonfile->Get("TAPS_phi");
	
    TH2F *proton_TAPS_energy = (TH2F*)protonfile->Get("TAPS_E");
    TH2F *proton_TAPS_theta = (TH2F*)protonfile->Get("TAPS_theta");
    TH2F *proton_TAPS_phi = (TH2F*)protonfile->Get("TAPS_phi");



	void SetFromBeam(const Double_t beam) {

		Ek = beam;
		Theta = 0;
		Phi = 0;
		Ek_Sigma = 1.9/TMath::Sqrt(3);
		Theta_Sigma = 1e-04;	
		Phi_Sigma = 1e-04;

	}

  	void	SetFromProton(const GTreeParticle *proton1, const Int_t particle1){

		aplconProtonTheta = proton1->Particle(particle1).Theta();
        	aplconProtonPhi = proton1->Particle(particle1).Phi();

		if(proton1->HasCB(particle1)){

			aplconProtonThetaSigma=4.35*TMath::DegToRad();
			aplconProtonPhiSigma=aplconProtonThetaSigma*0.7/TMath::Sin(aplconProtonTheta);

		}
		if(proton1->HasTAPS(particle1)){

			aplconProtonThetaSigma=1.15*TMath::DegToRad()*proton_TAPS_theta->GetBinContent(proton_TAPS_theta->FindBin(proton1->Particle(particle1).E()-proton1->Particle(particle1).M(),aplconProtonTheta*TMath::RadToDeg()));
			aplconProtonPhiSigma=0.765*TMath::DegToRad()*proton_TAPS_phi->GetBinContent(proton_TAPS_phi->FindBin(proton1->Particle(particle1).E()-proton1->Particle(particle1).M(),aplconProtonTheta*TMath::RadToDeg()));

			//if histo is not finding something
			if(aplconProtonThetaSigma==0) aplconProtonThetaSigma = 3.7*TMath::DegToRad();
	      	       	if(aplconProtonPhiSigma==0) aplconProtonPhiSigma = 1.2*TMath::DegToRad();

     	   	}

    	}
     	void	SetFromVector(const GTreeParticle *photon, const Int_t particle){

		Ek = photon->Particle(particle).E()-photon->Particle(particle).M();
		Theta = photon->Particle(particle).Theta();
		Phi = photon->Particle(particle).Phi();

		if(photon->HasCB(particle)){

			Ek_Sigma = Ek*1.5*photon_CB_energy->GetBinContent(photon_CB_energy->FindBin(Ek,Theta*TMath::RadToDeg()));
			Phi_Sigma = 0.83*TMath::DegToRad()*photon_CB_phi->GetBinContent(photon_CB_phi->FindBin(Ek,Theta*TMath::RadToDeg()));
			Theta_Sigma = 0.68*TMath::DegToRad()*photon_CB_theta->GetBinContent(photon_CB_theta->FindBin(Ek,Theta*TMath::RadToDeg()));

			//if histo is not finding something
			if(Ek_Sigma==0) Ek_Sigma=(0.02*(Ek/1.0e3)*pow((Ek/1.0e3),-0.36))*1.0e3*2.0;
			if(Theta_Sigma==0) Theta_Sigma = 5.2*TMath::DegToRad();
 			if(Phi_Sigma==0) Phi_Sigma = 0.7*Theta_Sigma/TMath::Sin(Theta);

		}

		if(photon->HasTAPS(particle)){

			Ek_Sigma = Ek*1*photon_TAPS_energy->GetBinContent(photon_TAPS_energy->FindBin(Ek,Theta*TMath::RadToDeg()));
			Phi_Sigma = 1*TMath::DegToRad();
			Theta_Sigma = 2.5*TMath::DegToRad();

			//if histo is not finding something
			if(Ek_Sigma==0) Ek_Sigma = ((0.018 + 0.008*TMath::Sqrt(Ek/1.0e3))*(Ek/1.0e3))*1.0e3*2.0;
			if(Theta_Sigma==0) Theta_Sigma = 2*TMath::DegToRad();
 			if(Phi_Sigma==0) Phi_Sigma = 2*TMath::DegToRad();

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

        static TLorentzVector Make(const Double_t Ek,
                                   const std::vector<double>& ThetaPhi,
                                   const Double_t m){

		const double E = Ek + m;
		const Double_t p = sqrt( E*E - m*m );
		TVector3 pv(1,0,0);
		pv.SetMagThetaPhi(p, ThetaPhi[0], ThetaPhi[1]);
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

	std::vector<APLCON::Variable_Settings_t> LinkSettings()
	{
        	return{Ek_Setting, Theta_Setting, Phi_Setting};
	}

    	std::vector<double*>    GetLink()       {return {&aplconProtonTheta, &aplconProtonPhi};}
	std::vector<double*>    GetLinkSigma()  {return {&aplconProtonThetaSigma, &aplconProtonPhiSigma};}


        double Ek;
        double Ek_Sigma;
        double Theta;
        double Theta_Sigma;
        double Phi;
        double Phi_Sigma;
 
        double              aplconProtonTheta;
        double              aplconProtonPhi;
        double              aplconProtonThetaSigma;
        double              aplconProtonPhiSigma;

   	APLCON::Variable_Settings_t Ek_Setting;
    	APLCON::Variable_Settings_t Theta_Setting;
    	APLCON::Variable_Settings_t Phi_Setting;


};
    	APLCON                      fitter;
	std::vector<FitParticle>    gamma;
	FitParticle		    proton;
// 	Double_t		    beam;
	FitParticle		    beam;

	APLCON::Result_t	    result;

public:
	GFit();
	~GFit();

	void	AddConstraintsTotMomentum();
	void	AddConstraintsTotEnergy();	
	void	AddConstraintsMM();
	void	AddConstraintsIM();
	bool    Solve();

	void	Set(const GTreeParticle *p0, const Int_t particle1, const Int_t particle2, const Int_t particle3,const Int_t particle4)
{
	gamma[0].SetFromVector(p0,particle1);
	gamma[1].SetFromVector(p0,particle2);
	gamma[2].SetFromVector(p0,particle3);
	gamma[3].SetFromVector(p0,particle4);
}

	void 	SetProton(const GTreeParticle *pr, const Int_t particle)	{proton.SetFromProton(pr,particle);}

//   	void    SetBeam(const Double_t _Beam)    {beam=_Beam;}

 	void    SetBeam(const double _Beam)    {beam.SetFromBeam(_Beam);}

const   GFit::FitParticle&   GetFittedProton() const
{
	    const APLCON::Result_Variable_t& pe = result.Variables.at("PE");
	    static FitParticle p;
	    p.Ek            = pe.Value.After;
	    p.Ek_Sigma      = pe.Sigma.After;
	    p.Theta         = aplconProtonTheta;
	    p.Theta_Sigma   = aplconProtonThetaSigma;
	    p.Phi         = aplconProtonPhi;
	    p.Phi_Sigma   = aplconProtonPhiSigma;
	    return p;
}


	TLorentzVector		GetFittedGamma(const int index)	{return FitParticle::Make(gamma[index], 0);}
	Double_t		GetFittedProtonTheta()		{return aplconProtonTheta;}
	Double_t		GetFittedProtonPhi()		{return aplconProtonPhi;}
	TLorentzVector		GetFittedBeam()		{return FitParticle::Make(beam, 0);}
	APLCON::Result_t	GetFittedResult()		{return result;}
};




#endif
