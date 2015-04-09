#include "GFit.h"
#include "GTreeMeson.h"



GFit::GFit(const char* name)   :
    GHistLinked(kFALSE),
    im(TString(name).Append("_im"), TString(name).Append("_im"), 500, 700, 1200, 48, kFALSE),
    zVertex(TString(name).Append("_zVertex"), TString(name).Append("_zVertex"), 200, -10, 10, 48, kFALSE),
    theta(TString(name).Append("_theta"), TString(name).Append("_theta"), 180, 0, 180, 48, kFALSE),
    phi(TString(name).Append("_phi"), TString(name).Append("_phi"), 360, -180, 180, 48, kFALSE),
    chiSq(TString(name).Append("_chiSq"), TString(name).Append("_chiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(name).Append("_confidenceLevel"), TString(name).Append("_confidenceLevel"), 1000, 0, 1, 48, kFALSE),
    name(name),
    fitter(name)
{
    aplconPhotons.resize(6);

    for(size_t i=0;i<aplconPhotons.size();i++) {
        std::stringstream s;
        s << "Photon" << i;
        fitter.LinkVariable(s.str(), aplconPhotons[i].Link(), aplconPhotons[i].LinkSigma());
    }

    fitter.AddUnmeasuredVariable("zVertex");

    APLCON::Fit_Settings_t settings = fitter.GetSettings();
    settings.MaxIterations = 100;
    fitter.SetSettings(settings);
}

void    GFit::CalcResult()
{
    im.CalcResult();
    zVertex.CalcResult();
    theta.CalcResult();
    phi.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();
}

void    GFit::PrepareWriteList(GHistWriteList* arr, const char* _Name)
{
    if(!arr)
        return;

    im.PrepareWriteList(arr, "IM");
    zVertex.PrepareWriteList(arr, "zVertex");
    theta.PrepareWriteList(arr, "theta");
    phi.PrepareWriteList(arr, "phi");
    chiSq.PrepareWriteList(arr, "chiSq");
    confidenceLevel.PrepareWriteList(arr, "confidenceLevel");
}

void    GFit::Reset(Option_t* option)
{
    im.Reset(option);
    zVertex.Reset(option);
    theta.Reset(option);
    phi.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);
}

void GFit::Solve(const double time, const int channel)
{
    result = fitter.DoFit();
    if(result.Status == APLCON::Result_Status_t::Success)
    {
        TLorentzVector etap(0,0,0,0);
        for(size_t t=0;t<aplconPhotons.size();t++)
            etap += FitParticle::Make(aplconPhotons[t], 0);
        im.Fill(etap.M(), time, channel);
        const APLCON::Result_Variable_t& var = result.Variables.at("zVertex");
        zVertex.Fill(var.Value.After, time, channel);
        theta.Fill(etap.Theta()*TMath::RadToDeg(), time, channel);
        phi.Fill(etap.Phi()*TMath::RadToDeg(), time);
        chiSq.Fill(result.ChiSquare, time);
        confidenceLevel.Fill(TMath::Prob(result.ChiSquare, result.NDoF), time);
        //cout << result << endl;
        return;
    }
    //cout << "fit " << name << " not working" << endl;
}











void    GFit1Constraints::AddConstraints()
{
    fitter.AddConstraint("mm", {"Photon0", "Photon1", "Photon2", "Photon3", "Photon4", "Photon5", "zVertex"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        p0[1]    = std::atan2( R*sin(p0[1]), R*cos(p0[1]) - v_z);
        p1[1]    = std::atan2( R*sin(p1[1]), R*cos(p1[1]) - v_z);
        p2[1]    = std::atan2( R*sin(p2[1]), R*cos(p2[1]) - v_z);
        p3[1]    = std::atan2( R*sin(p3[1]), R*cos(p3[1]) - v_z);
        p4[1]    = std::atan2( R*sin(p4[1]), R*cos(p4[1]) - v_z);
        p5[1]    = std::atan2( R*sin(p5[1]), R*cos(p5[1]) - v_z);

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);

        return (beamAndTarget-v0-v1-v2-v3-v4-v5).M() - MASS_PROTON;
        }
    );
}

void GFit1Constraints::Solve(const double time, const int channel)
{
    GFit::Solve(time, channel);
    for(int i=0; i<6; i++)
    {
        for(int p=0; p<3; p++)
        {
            std::stringstream s;
            s << "Photon" << i;
            s << "[" << p << "]";
            //cout << result << endl;
            //cout << s.str() << endl;
            const APLCON::Result_Variable_t& var = result.Variables.at(s.str());
            pulls.Fill(var.Pull, (3*i)+p);
        }
    }
}

void    GFit1Constraints::CalcResult()
{
    GFit::CalcResult();
    pulls.CalcResult();
}

void    GFit1Constraints::PrepareWriteList(GHistWriteList* arr, const char* _Name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory(name);
    GFit::PrepareWriteList(folder);
    pulls.PrepareWriteList(folder, "pulls");
}

void    GFit1Constraints::Reset(Option_t* option)
{
    GFit::Reset(option);
    pulls.Reset(option);
}






void    GFit3Constraints::AddConstraints()
{
    fitter.AddConstraint("im0", {"Photon0", "Photon1", "zVertex"},
                        [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_ETA;
    }
    );

    fitter.AddConstraint("im1", {"Photon2", "Photon3", "zVertex"},
                        [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );

    fitter.AddConstraint("im2", {"Photon4", "Photon5", "zVertex"},
                        [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );
}









void    GFit4Constraints::AddConstraints()
{
    fitter.AddConstraint("im0", {"Photon0", "Photon1", "zVertex"},
    [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_ETA;
    }
    );

    fitter.AddConstraint("im1", {"Photon2", "Photon3", "zVertex"},
    [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );

    fitter.AddConstraint("im2", {"Photon4", "Photon5", "zVertex"},
    [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );
}




GFitProton::GFitProton() :
    GFit("fit7")
{
    fitter.LinkVariable(std::string("Proton"), proton.Link(), proton.LinkSigma());

    fitter.AddConstraint("im0", {"Photon0", "Photon1", "zVertex"},
    [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_ETA;
    }
    );

    fitter.AddConstraint("im1", {"Photon2", "Photon3", "zVertex"},
    [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );

    fitter.AddConstraint("im2", {"Photon4", "Photon5", "zVertex"},
    [] (vector<double>& p0, vector<double>& p1, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        const double theta0 = p0[1]; // second element is theta
        const double theta0_p = std::atan2( R*sin(theta0), R*cos(theta0) - v_z);
        const double theta1 = p1[1]; // second element is theta
        const double theta1_p = std::atan2( R*sin(theta1), R*cos(theta1) - v_z);
        p0[1] = theta0_p;
        p1[1] = theta1_p;

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );

    fitter.AddConstraint("Px", {"Photon0", "Photon1", "Photon2", "Photon3", "Photon4", "Photon5", "Proton", "zVertex"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pr, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        p0[1]    = std::atan2( R*sin(p0[1]), R*cos(p0[1]) - v_z);
        p1[1]    = std::atan2( R*sin(p1[1]), R*cos(p1[1]) - v_z);
        p2[1]    = std::atan2( R*sin(p2[1]), R*cos(p2[1]) - v_z);
        p3[1]    = std::atan2( R*sin(p3[1]), R*cos(p3[1]) - v_z);
        p4[1]    = std::atan2( R*sin(p4[1]), R*cos(p4[1]) - v_z);
        p5[1]    = std::atan2( R*sin(p5[1]), R*cos(p5[1]) - v_z);
        pr[1]    = std::atan2( R*sin(pr[1]), R*cos(pr[1]) - v_z);

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector pro= FitParticle::Make(pr, 0);

        return beamAndTarget.Px()-v0.Px()-v1.Px()-v2.Px()-v3.Px()-v4.Px()-v5.Px()-pro.Px();
        }
    );

    fitter.AddConstraint("Py", {"Photon0", "Photon1", "Photon2", "Photon3", "Photon4", "Photon5", "Proton", "zVertex"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pr, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        p0[1]    = std::atan2( R*sin(p0[1]), R*cos(p0[1]) - v_z);
        p1[1]    = std::atan2( R*sin(p1[1]), R*cos(p1[1]) - v_z);
        p2[1]    = std::atan2( R*sin(p2[1]), R*cos(p2[1]) - v_z);
        p3[1]    = std::atan2( R*sin(p3[1]), R*cos(p3[1]) - v_z);
        p4[1]    = std::atan2( R*sin(p4[1]), R*cos(p4[1]) - v_z);
        p5[1]    = std::atan2( R*sin(p5[1]), R*cos(p5[1]) - v_z);

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector pro= FitParticle::Make(pr, 0);

        return beamAndTarget.Py()-v0.Py()-v1.Py()-v2.Py()-v3.Py()-v4.Py()-v5.Py()-pro.Py();
        }
    );

    fitter.AddConstraint("Pz", {"Photon0", "Photon1", "Photon2", "Photon3", "Photon4", "Photon5", "Proton", "zVertex"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pr, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        p0[1]    = std::atan2( R*sin(p0[1]), R*cos(p0[1]) - v_z);
        p1[1]    = std::atan2( R*sin(p1[1]), R*cos(p1[1]) - v_z);
        p2[1]    = std::atan2( R*sin(p2[1]), R*cos(p2[1]) - v_z);
        p3[1]    = std::atan2( R*sin(p3[1]), R*cos(p3[1]) - v_z);
        p4[1]    = std::atan2( R*sin(p4[1]), R*cos(p4[1]) - v_z);
        p5[1]    = std::atan2( R*sin(p5[1]), R*cos(p5[1]) - v_z);

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector pro= FitParticle::Make(pr, 0);

        return beamAndTarget.Pz()-v0.Pz()-v1.Pz()-v2.Pz()-v3.Pz()-v4.Pz()-v5.Pz()-pro.Pz();
        }
    );

    fitter.AddConstraint("E", {"Photon0", "Photon1", "Photon2", "Photon3", "Photon4", "Photon5", "Proton", "zVertex"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pr, const vector<double>& zv) {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = zv[0];
        p0[1]    = std::atan2( R*sin(p0[1]), R*cos(p0[1]) - v_z);
        p1[1]    = std::atan2( R*sin(p1[1]), R*cos(p1[1]) - v_z);
        p2[1]    = std::atan2( R*sin(p2[1]), R*cos(p2[1]) - v_z);
        p3[1]    = std::atan2( R*sin(p3[1]), R*cos(p3[1]) - v_z);
        p4[1]    = std::atan2( R*sin(p4[1]), R*cos(p4[1]) - v_z);
        p5[1]    = std::atan2( R*sin(p5[1]), R*cos(p5[1]) - v_z);

        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector pro= FitParticle::Make(pr, 0);

        return beamAndTarget.E()-v0.E()-v1.E()-v2.E()-v3.E()-v4.E()-v5.E()-pro.E();
        }
    );
}

/*
GHistFit::GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    nPulls(_NPulls),
    im(TString(name).Append("im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, 48, kFALSE),
    sub0im(TString(name).Append("sub0im"), TString(title).Append(" sub0 inv. Mass"), 800, 0, 800, 48, kFALSE),
    sub1im(TString(name).Append("sub1im"), TString(title).Append(" sub1 inv. Mass"), 400, 0, 400, 48, kFALSE),
    sub2im(TString(name).Append("sub2im"), TString(title).Append(" sub2 inv. Mass"), 400, 0, 400, 48, kFALSE),
    theta(TString(name).Append("theta"), TString(title).Append(" theta"), 180, 0, 180, 48, kFALSE),
    phi(TString(name).Append("phi"), TString(title).Append(" phi"), 360, -180, 180, 48, kFALSE),
    Pim(TString(name).Append("Pim"), TString(title).Append(" Proton Mass"), 2000, 0, 2000, 48, kFALSE),
    Ptheta(TString(name).Append("Ptheta"), TString(title).Append(" Proton theta"), 180, 0, 180, 48, kFALSE),
    Pphi(TString(name).Append("Pphi"), TString(title).Append(" Proton phi"), 360, -180, 180, 48, kFALSE),
    chiSq(TString(name).Append("ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, 48, kFALSE),
    VchiSq(TString(name).Append("VChiSq"), TString(title).Append(" VChiSq"), 1000, 0, 100, 48, kFALSE),
    CchiSq(TString(name).Append("CChiSq"), TString(title).Append(" CChiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(name).Append("ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE),
    VconfidenceLevel(TString(name).Append("VConfLev"), TString(title).Append(" VConfLev"), 1000, 0, 1, 48, kFALSE),
    CconfidenceLevel(TString(name).Append("CConfLev"), TString(title).Append(" CConfLev"), 1000, 0, 1, 48, kFALSE),
    pulls(TString(name).Append("Pulls"), TString(title).Append(" Pulls"), 100, -5, 5, nPulls, 0, nPulls, kFALSE)
{

}

GHistFit::~GHistFit()
{

}

void        GHistFit::CalcResult()
{
    im.CalcResult();
    sub0im.CalcResult();
    sub1im.CalcResult();
    sub2im.CalcResult();
    theta.CalcResult();
    phi.CalcResult();
    Pim.CalcResult();
    Ptheta.CalcResult();
    Pphi.CalcResult();
    chiSq.CalcResult();
    VchiSq.CalcResult();
    CchiSq.CalcResult();
    confidenceLevel.CalcResult();
    VconfidenceLevel.CalcResult();
    CconfidenceLevel.CalcResult();
    pulls.CalcResult();
}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime)
{
    TLorentzVector  etap(fitter.GetMeson());
    Int_t   ret = im.Fill(etap.M(), taggerTime);
    ret += sub0im.Fill(fitter.GetSub(0).M(), taggerTime);
    ret += sub1im.Fill(fitter.GetSub(1).M(), taggerTime);
    ret += sub2im.Fill(fitter.GetSub(2).M(), taggerTime);
    ret += theta.Fill(etap.Theta()*TMath::RadToDeg(), taggerTime);
    ret += phi.Fill(etap.Phi()*TMath::RadToDeg(), taggerTime);
    ret += Pim.Fill(fitter.GetRecoil().M(), taggerTime);
    ret += Ptheta.Fill(fitter.GetRecoil().Theta()*TMath::RadToDeg(), taggerTime);
    ret += Pphi.Fill(fitter.GetRecoil().Phi()*TMath::RadToDeg(), taggerTime);
    ret += chiSq.Fill(fitter.GetChi2(), taggerTime);
    ret += VchiSq.Fill(fitter.GetVariablesChi2(), taggerTime);
    ret += CchiSq.Fill(fitter.GetConstraintsChi2(), taggerTime);
    ret += confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime);
    ret += VconfidenceLevel.Fill(fitter.VariablesConfidenceLevel(), taggerTime);
    ret += CconfidenceLevel.Fill(fitter.ConstraintsConfidenceLevel(), taggerTime);
    for(int i=0; i<nPulls; i++)
        ret += pulls.Fill(fitter.GetPull(i), i);
    return ret;
}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel)
{
    TLorentzVector  etap(fitter.GetMeson());
    Int_t   ret = im.Fill(etap.M(), taggerTime, taggerChannel);
    ret += sub0im.Fill(fitter.GetSub(0).M(), taggerTime, taggerChannel);
    ret += sub1im.Fill(fitter.GetSub(1).M(), taggerTime, taggerChannel);
    ret += sub2im.Fill(fitter.GetSub(2).M(), taggerTime, taggerChannel);
    ret += theta.Fill(etap.Theta()*TMath::RadToDeg(), taggerTime, taggerChannel);
    ret += phi.Fill(etap.Phi()*TMath::RadToDeg(), taggerTime, taggerChannel);
    ret += Pim.Fill(fitter.GetRecoil().M(), taggerTime, taggerChannel);
    ret += Ptheta.Fill(fitter.GetRecoil().Theta()*TMath::RadToDeg(), taggerTime, taggerChannel);
    ret += Pphi.Fill(fitter.GetRecoil().Phi()*TMath::RadToDeg(), taggerTime, taggerChannel);
    ret += chiSq.Fill(fitter.GetChi2(), taggerTime, taggerChannel);
    ret += VchiSq.Fill(fitter.GetVariablesChi2(), taggerTime, taggerChannel);
    ret += CchiSq.Fill(fitter.GetConstraintsChi2(), taggerTime, taggerChannel);
    ret += confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime, taggerChannel);
    ret += VconfidenceLevel.Fill(fitter.VariablesConfidenceLevel(), taggerTime, taggerChannel);
    ret += CconfidenceLevel.Fill(fitter.ConstraintsConfidenceLevel(), taggerTime, taggerChannel);
    for(int i=0; i<nPulls; i++)
        ret += pulls.Fill(fitter.GetPull(i), i);
    return ret;
}

void    GHistFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        sub0im.PrepareWriteList(arr, TString(name).Append("_Sub0IM").Data());
        sub1im.PrepareWriteList(arr, TString(name).Append("_Sub1IM").Data());
        sub2im.PrepareWriteList(arr, TString(name).Append("_Sub2IM").Data());
        theta.PrepareWriteList(arr, TString(name).Append("_Theta").Data());
        phi.PrepareWriteList(arr, TString(name).Append("_Phi").Data());
        Pim.PrepareWriteList(arr, TString(name).Append("_PrIM").Data());
        Ptheta.PrepareWriteList(arr, TString(name).Append("_PrTheta").Data());
        Pphi.PrepareWriteList(arr, TString(name).Append("_PrPhi").Data());
        chiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
        VchiSq.PrepareWriteList(arr, TString(name).Append("_VChiSq").Data());
        CchiSq.PrepareWriteList(arr, TString(name).Append("_CChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        VconfidenceLevel.PrepareWriteList(arr, TString(name).Append("_VConfLev").Data());
        CconfidenceLevel.PrepareWriteList(arr, TString(name).Append("_CConfLev").Data());
        pulls.PrepareWriteList(arr, TString(name).Append("_Pulls").Data());
    }
    else
    {
        im.PrepareWriteList(arr);
        sub0im.PrepareWriteList(arr);
        sub1im.PrepareWriteList(arr);
        sub2im.PrepareWriteList(arr);
        theta.PrepareWriteList(arr);
        phi.PrepareWriteList(arr);
        Pim.PrepareWriteList(arr);
        Ptheta.PrepareWriteList(arr);
        Pphi.PrepareWriteList(arr);
        chiSq.PrepareWriteList(arr);
        VchiSq.PrepareWriteList(arr);
        CchiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        VconfidenceLevel.PrepareWriteList(arr);
        CconfidenceLevel.PrepareWriteList(arr);
        pulls.PrepareWriteList(arr);
    }
}

void        GHistFit::Reset(Option_t* option)
{
    im.Reset(option);
    sub0im.Reset(option);
    sub1im.Reset(option);
    sub2im.Reset(option);
    theta.Reset(option);
    phi.Reset(option);
    Pim.Reset(option);
    Ptheta.Reset(option);
    Pphi.Reset(option);
    chiSq.Reset(option);
    VchiSq.Reset(option);
    CchiSq.Reset(option);
    confidenceLevel.Reset(option);
    VconfidenceLevel.Reset(option);
    CconfidenceLevel.Reset(option);
    pulls.Reset(option);
}

void        GHistFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub0im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub1im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub2im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    theta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    phi.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Pim.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Ptheta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Pphi.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    VchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    VconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}









GHistIterativeFit::GHistIterativeFit(const char* name, const char* title, const Int_t _NPulls, const Int_t _NSteps, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    im(TString(name).Append("im"), TString(title).Append(" inv. Mass"), 750, 500, 1250, _NSteps, 0, _NSteps, kFALSE),
    sub0im(TString(name).Append("sub0im"), TString(title).Append(" sub0 inv. Mass"), 1000, 500, 600, _NSteps, 0, _NSteps, kFALSE),
    sub1im(TString(name).Append("sub1im"), TString(title).Append(" sub1 inv. Mass"), 1000, 85, 185, _NSteps, 0, _NSteps, kFALSE),
    sub2im(TString(name).Append("sub2im"), TString(title).Append(" sub2 inv. Mass"), 1000, 85, 185, _NSteps, 0, _NSteps, kFALSE),
    mm(TString(name).Append("mm"), TString(title).Append(" mm"), 10500, -50, 1000, _NSteps, 0, _NSteps, kFALSE),
    totE(TString(name).Append("totE"), TString(title).Append(" totE"), 1000, -100, 1900, _NSteps, 0, _NSteps, kFALSE),
    totPx(TString(name).Append("totPx"), TString(title).Append(" totPx"), 1000, -500, 500, _NSteps, 0, _NSteps, kFALSE),
    totPy(TString(name).Append("totPy"), TString(title).Append(" totPy"), 1000, -500, 500, _NSteps, 0, _NSteps, kFALSE),
    totPz(TString(name).Append("totPz"), TString(title).Append(" totPz"), 1000, -100, 1900, _NSteps, 0, _NSteps, kFALSE),
    VchiSq(TString(name).Append("VChiSq"), TString(title).Append(" VChiSq"), 1000, 0, 100, _NSteps, 0, _NSteps, kFALSE),
    VconfidenceLevel(TString(name).Append("VConfLev"), TString(title).Append(" VConfLev"), 1000, 0, 1, _NSteps, 0, _NSteps, kFALSE),
    CchiSq(TString(name).Append("CChiSq"), TString(title).Append(" CChiSq"), 1000, 0, 100, _NSteps, 0, _NSteps, kFALSE),
    CconfidenceLevel(TString(name).Append("CConfLev"), TString(title).Append(" CConfLev"), 1000, 0, 1, _NSteps, 0, _NSteps, kFALSE),
    chiSq(TString(name).Append("ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, _NSteps, 0, _NSteps, kFALSE),
    confidenceLevel(TString(name).Append("ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, _NSteps, 0, _NSteps, kFALSE),
    final(TString(name).Append("Final"), TString(title).Append(" Final"), _NPulls, kFALSE)
{

}

GHistIterativeFit::~GHistIterativeFit()
{

}


void        GHistIterativeFit::CalcResult()
{
    //std::cout << im.GetName() << std::endl;
    im.CalcResult();
    sub0im.CalcResult();
    sub1im.CalcResult();
    sub2im.CalcResult();
    mm.CalcResult();
    totE.CalcResult();
    totPx.CalcResult();
    totPy.CalcResult();
    totPz.CalcResult();
    VchiSq.CalcResult();
    VconfidenceLevel.CalcResult();
    CchiSq.CalcResult();
    CconfidenceLevel.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();

    final.CalcResult();
}

Int_t       GHistIterativeFit::Fill(GFit& fitter)
{
    TLorentzVector  etap(fitter.GetMeson());
    TLorentzVector  tot(fitter.GetTotal());
    Int_t   ret = im.Fill(etap.M(), fitter.GetIterations());
    ret += sub0im.Fill(fitter.GetSub(0).M(), fitter.GetIterations());
    ret += sub1im.Fill(fitter.GetSub(1).M(), fitter.GetIterations());
    ret += sub2im.Fill(fitter.GetSub(2).M(), fitter.GetIterations());
    ret += mm.Fill(tot.M(), fitter.GetIterations());
    ret += totE.Fill(tot.E(), fitter.GetIterations());
    ret += totPx.Fill(tot.Px(), fitter.GetIterations());
    ret += totPy.Fill(tot.Py(), fitter.GetIterations());
    ret += totPz.Fill(tot.Pz(), fitter.GetIterations());
    ret += VchiSq.Fill(fitter.GetVariablesChi2(), fitter.GetIterations());
    ret += VconfidenceLevel.Fill(fitter.VariablesConfidenceLevel(), fitter.GetIterations());
    ret += CchiSq.Fill(fitter.GetConstraintsChi2(), fitter.GetIterations());
    ret += CconfidenceLevel.Fill(fitter.ConstraintsConfidenceLevel(), fitter.GetIterations());
    ret += chiSq.Fill(fitter.GetChi2(), fitter.GetIterations());
    ret += confidenceLevel.Fill(fitter.ConfidenceLevel(), fitter.GetIterations());
    return ret;
}

void    GHistIterativeFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

        im.PrepareWriteList(arr, "IM");
        sub0im.PrepareWriteList(arr, "_Sub0IM");
        sub1im.PrepareWriteList(arr, "_Sub1IM");
        sub2im.PrepareWriteList(arr, "_Sub2IM");
        mm.PrepareWriteList(arr, "_MM");
        totE.PrepareWriteList(arr, "_TotE");
        totPx.PrepareWriteList(arr, "_TotPx");
        totPy.PrepareWriteList(arr, "_TotPy");
        totPz.PrepareWriteList(arr, "_TotPz");
        VchiSq.PrepareWriteList(arr, "_VChiSq");
        VconfidenceLevel.PrepareWriteList(arr, "_VConfLev");
        CchiSq.PrepareWriteList(arr, "_CChiSq");
        CconfidenceLevel.PrepareWriteList(arr, "_CConfLev");
        chiSq.PrepareWriteList(arr, "_ChiSq");
        confidenceLevel.PrepareWriteList(arr, "_ConfLev");

    GHistWriteList* finalDir    = arr->GetDirectory("Final");
    final.PrepareWriteList(finalDir, "Final");
}

void        GHistIterativeFit::Reset(Option_t* option)
{
    im.Reset(option);
    sub0im.Reset(option);
    sub1im.Reset(option);
    sub2im.Reset(option);
    mm.Reset(option);
    totE.Reset(option);
    totPx.Reset(option);
    totPy.Reset(option);
    totPz.Reset(option);
    VchiSq.Reset(option);
    VconfidenceLevel.Reset(option);
    CchiSq.Reset(option);
    CconfidenceLevel.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);

    final.Reset(option);
}

void        GHistIterativeFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub0im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub1im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub2im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    mm.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totE.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totPx.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totPy.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totPz.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    VchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    VconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);

    final.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}*/
