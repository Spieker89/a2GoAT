#include "GFit.h"



using namespace std;




GFit::GFit()	:
	fitter("fitter")
{
	gamma.resize(4);

	for(size_t i=0;i<gamma.size();i++) {
        	std::stringstream s;
        	s << "Ph" << i;
        	fitter.LinkVariable(s.str(), gamma[i].Link(), gamma[i].LinkSigma());
    	}

	fitter.LinkVariable("Pr", proton.Link(), proton.LinkSigma());
	fitter.LinkVariable("Be", beam.Link(), beam.LinkSigma());
    	APLCON::Fit_Settings_t settings = fitter.GetSettings();
    	settings.MaxIterations = 100;
    	fitter.SetSettings(settings);
}

GFit::~GFit()
{
}





void    GFit::AddConstraintsTotMomentum()
{
    fitter.AddConstraint("px", {"Be","Ph0", "Ph1", "Ph2", "Ph3", "Pr"},
                        [&] (vector<double>& be,vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& pr) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector vp = FitParticle::Make(pr, 938);

        return vb.Px()-v0.Px()-v1.Px()-v2.Px()-v3.Px()-vp.Px();
        }
    );



    fitter.AddConstraint("py", {"Be","Ph0", "Ph1", "Ph2", "Ph3", "Pr"},
                        [&] (vector<double>& be,vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& pr) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector vp = FitParticle::Make(pr, 938);

        return vb.Py()-v0.Py()-v1.Py()-v2.Py()-v3.Py()-vp.Py();
        }
    );


    fitter.AddConstraint("pz", {"Be","Ph0", "Ph1", "Ph2", "Ph3", "Pr"},
                        [&] (vector<double>& be,vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& pr) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector vp = FitParticle::Make(pr, 938);

        return vb.Pz()-v0.Pz()-v1.Pz()-v2.Pz()-v3.Pz()-vp.Pz();
        }
    );
}

void    GFit::AddConstraintsTotEnergy()
{
    fitter.AddConstraint("e", {"Be","Ph0", "Ph1", "Ph2", "Ph3", "Pr"},
                        [&] (vector<double>& be,vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& pr) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector vp = FitParticle::Make(pr, 938);

        return vb.E()+938.272046-v0.E()-v1.E()-v2.E()-v3.E()-vp.E();
        }
    );
}

void    GFit::AddConstraintsMM()
{
    fitter.AddConstraint("mm", {"Be","Ph0", "Ph1", "Ph2", "Ph3"},
                        [&] (vector<double>& be,vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);

	TLorentzVector proton_target(0.,0.,0.,938.272046);

        return (vb + proton_target-v0-v1-v2-v3).M()-938.272046;
        }
    );
}

void    GFit::AddConstraintsIM()
{
    fitter.AddConstraint("im0", {"Ph0", "Ph1"},
                        [&] (vector<double>& p0, vector<double>& p1) {

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);


        return (v0+v1).M()-134.9766;
        }
    );

    fitter.AddConstraint("im1", {"Ph2", "Ph3"},
                        [&] (vector<double>& p0, vector<double>& p1) {

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);


        return (v0+v1).M()-134.9766;
        }
    );
}


bool GFit::Solve()
{
    result = fitter.DoFit();
    if(result.Status == APLCON::Result_Status_t::Success)
    {
        Double_t    cl  = TMath::Prob(result.ChiSquare, 15-result.NDoF);
	
        return true;
    }
    //cout << name << " not working" << endl;
    return false;
}
