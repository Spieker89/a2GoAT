#ifndef     ___GFIT___
#define	    ___GFIT___


#include <vector>
#include <TLorentzVector.h>
#include "/disk/nobackup/001/spieker/Mainz1/APLCON/src/APLCON.hpp"



class	GFit
{
private:

    double              Theta;
    double              Phi;
    double              Theta_Sigma;
    double              Phi_Sigma;

    double              aplconProtonTheta;
    double              aplconProtonPhi;
    double              aplconProtonThetaSigma;
    double              aplconProtonPhiSigma;

	struct FitParticle {

	void SetFromBeam(const Double_t beam) {
            Ek = beam;
            Theta = 0;
            Phi = 0;
            Ek_Sigma = 1.5;
            Theta_Sigma = 2*TMath::DegToRad();
            Phi_Sigma = 6*TMath::DegToRad();
        }

     	void    SetFromProton(const TLorentzVector proton){
        	aplconProtonTheta = proton.Theta();
        	aplconProtonPhi = proton.Phi();
        	Theta_Sigma = 5*TMath::DegToRad();
     	   if(aplconProtonTheta>20*TMath::DegToRad() && aplconProtonTheta<160*TMath::DegToRad()) {
     	       aplconProtonPhiSigma = aplconProtonThetaSigma/sin(aplconProtonTheta);
     	   }
     	   else {
     	       aplconProtonPhiSigma = 4*TMath::DegToRad();
     	   }
    	}

        void SetFromVector(const TLorentzVector& p_) {
            Ek = p_.E()-p_.M();
            Theta = p_.Theta();
            Phi = p_.Phi();
            Ek_Sigma = Ek*(0.018+0.008*pow(Ek,-0.8));
            Theta_Sigma = 4*TMath::DegToRad();
            if(Theta>20*TMath::DegToRad() && Theta<160*TMath::DegToRad()) {
                Phi_Sigma = Theta_Sigma/sin(Theta);
           	 Ek_Sigma = 0.02*Ek*pow(Ek,-0.5);
            }
            else {
                Phi_Sigma = 5*TMath::DegToRad();
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

        std::vector<double*> LinkProton() {
            return {
                        std::addressof(aplconProtonTheta),
                        std::addressof(aplconProtonPhi)};
        }
        std::vector<double*> LinkSigmaProton() {
            return {
                        std::addressof(aplconProtonThetaSigma),
                        std::addressof(aplconProtonPhiSigma)};
        }


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


    };

    	APLCON                      fitter;
	std::vector<FitParticle>    gamma;
	FitParticle		    proton;
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



	void	Set(const TLorentzVector& p0, const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3)
{
	gamma[0].SetFromVector(p0);
	gamma[1].SetFromVector(p1);
	gamma[2].SetFromVector(p2);
	gamma[3].SetFromVector(p3);
}
	void 	SetProton(const TLorentzVector &pr)	{proton.SetFromProton(pr);}

 	void    SetBeam(const double _Beam)    {beam.SetFromBeam(_Beam);}

	TLorentzVector		GetFittedGamma(const int index)	{return FitParticle::Make(gamma[index], 0);}
	Double_t		GetFittedProtonTheta()		{return aplconProtonTheta;}
	Double_t		GetFittedProtonPhi()		{return aplconProtonPhi;}
	TLorentzVector		GetFittedBeam()		{return FitParticle::Make(beam, 0);}
	APLCON::Result_t	GetFittedResult()		{return result;}
};




#endif
