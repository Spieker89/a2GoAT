#ifndef     ___GFIT___
#define	    ___GFIT___


#include <vector>
#include <TLorentzVector.h>
#include "/disk/nobackup/001/spieker/Mainz1/APLCON/src/APLCON.hpp"



class	GFit
{
private:
	struct FitParticle {
        void SetFromVector(const TLorentzVector& p_) {
            Ek = p_.E()-p_.M();
            Theta = p_.Theta();
            Phi = p_.Phi();
            Ek_Sigma = 0.01*Ek*pow(Ek,-0.36);
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

	void	Set(const TLorentzVector& p0, const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3, const TLorentzVector& pr, const TLorentzVector& _Beam)
{
	gamma[0].SetFromVector(p0);
	gamma[1].SetFromVector(p1);
	gamma[2].SetFromVector(p2);
	gamma[3].SetFromVector(p3);
	proton.SetFromVector(pr);
	beam.SetFromVector(_Beam);
}

	TLorentzVector		GetFittedGamma(const int index)	{return FitParticle::Make(gamma[index], 0);}
	TLorentzVector		GetFittedProton()		{return FitParticle::Make(proton, 0);}
	TLorentzVector		GetFittedBeam()		{return FitParticle::Make(beam, 0);}
	APLCON::Result_t	GetFittedResult()		{return result;}
};




#endif
