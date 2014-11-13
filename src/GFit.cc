#include "GFit.h"
#include "GTreeMeson.h"


GFitStruct::GFitStruct()   :
    im(0),
    ChiSq(0),
    ConfidenceLevel(0)
{
    for(int i=0; i<24; i++)
        PullPhotons[i]     = 0;
    for(int i=0; i<4; i++)
        PullBeam[i]        = 0;
    for(int i=0; i<4; i++)
        PullProton[i]      = 0;
}

/*GFitPulls4Vector::GFitPulls4Vector(const char* name, const char* title) :
    Pull_Px(TString(name).Append("_Px"), TString(title).Append(" Px"), 2000, -10, 10, 48, kFALSE),
    Pull_Py(TString(name).Append("_Py"), TString(title).Append(" Py"), 2000, -10, 10, 48, kFALSE),
    Pull_Pz(TString(name).Append("_Pz"), TString(title).Append(" Pz"), 2000, -10, 10, 48, kFALSE),
    Pull_E(TString(name).Append("_E"), TString(title).Append(" E"), 2000, -10, 10, 48, kFALSE)
{

}

GFitPulls4Vector::~GFitPulls4Vector()
{

}


GFitPulls6Photons::GFitPulls6Photons(const char* name, const char* title) :
    g0(TString(name).Append("_g0"), TString(title).Append(" Photon 0")),
    g1(TString(name).Append("_g1"), TString(title).Append(" Photon 1")),
    g2(TString(name).Append("_g2"), TString(title).Append(" Photon 2")),
    g3(TString(name).Append("_g3"), TString(title).Append(" Photon 3")),
    g4(TString(name).Append("_g4"), TString(title).Append(" Photon 4")),
    g5(TString(name).Append("_g5"), TString(title).Append(" Photon 5"))
{

}

GFitPulls6Photons::~GFitPulls6Photons()
{

}
*/


TFile*  GFit3Constraints::GammaResFile  = 0;
TH2F*   GFit3Constraints::GammaEloss    = 0;
TH2F*   GFit3Constraints::GammaERes     = 0;
TH2F*   GFit3Constraints::GammaThetaRes = 0;
TH2F*   GFit3Constraints::GammaPhiRes   = 0;




GFit3Constraints::GFit3Constraints(const Int_t npart, const Int_t ncon, const Bool_t IsEtap) :
    isEtap(IsEtap),
    fitter(npart, ncon, 0)
{
    if(!GammaResFile)   GammaResFile   = new TFile("~/GammaRes.root");
    if(!GammaEloss)     GammaEloss     = (TH2F*)GammaResFile->Get("Eloss");
    if(!GammaERes)      GammaERes      = (TH2F*)GammaResFile->Get("EResIter");
    if(!GammaThetaRes)  GammaThetaRes  = (TH2F*)GammaResFile->Get("ThetaRes;1");
    if(!GammaPhiRes)    GammaPhiRes    = (TH2F*)GammaResFile->Get("PhiRes;1");
}

GFit3Constraints::GFit3Constraints(const Bool_t IsEtap) :
    isEtap(IsEtap),
    fitter(6, 3, 0)
{
    if(!GammaResFile)   GammaResFile   = new TFile("~/GammaRes.root");
    if(!GammaEloss)     GammaEloss     = (TH2F*)GammaResFile->Get("Eloss");
    if(!GammaERes)      GammaERes      = (TH2F*)GammaResFile->Get("EResIter");
    if(!GammaThetaRes)  GammaThetaRes  = (TH2F*)GammaResFile->Get("ThetaRes;1");
    if(!GammaPhiRes)    GammaPhiRes    = (TH2F*)GammaResFile->Get("PhiRes;1");
}

GFit3Constraints::~GFit3Constraints()
{

}

void    GFit3Constraints::InitFit(const GTreeMeson& meson)
{
    SetPhotons(meson);

    fitter.Reset();
    fitter.AddPosKFParticle(photons[0]);
    fitter.AddPosKFParticle(photons[1]);
    fitter.AddPosKFParticle(photons[2]);
    fitter.AddPosKFParticle(photons[3]);
    fitter.AddPosKFParticle(photons[4]);
    fitter.AddPosKFParticle(photons[5]);

    Int_t	sub[6] = {0, 1, 2, 3, 4, 5};
    if(isEtap)
        fitter.AddSubInvMassConstraint(2, &sub[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &sub[0], MASS_PI0);

    fitter.AddSubInvMassConstraint(2, &sub[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &sub[4], MASS_PI0);
}

Bool_t  GFit3Constraints::Fit(const GTreeMeson& meson)
{
    InitFit(meson);

    if(fitter.Solve()>=0)
    {
        result.im               = fitter.GetTotalFitParticle().Get4Vector().M();
        result.ChiSq            = fitter.GetChi2();
        result.ConfidenceLevel  = fitter.ConfidenceLevel();
        for(int i=0; i<24; i++)
            result.PullPhotons[i]   = fitter.Pull(i);
        return kTRUE;
    }
    return kFALSE;
}

void    GFit3Constraints::SetPhotons(const GTreeMeson& meson)
{
    Int_t Ebin  = 0;
    Int_t Thbin = 0;
    Float_t resth = 0;
    Float_t resph = 0;
    Float_t resE  = 0;

    for(int i=0; i<6; i++)
    {
        Ebin  = GammaEloss->GetXaxis()->FindFixBin(meson.SubPhotons(0, i).E());
        Thbin = GammaEloss->GetYaxis()->FindFixBin(meson.SubPhotons(0, i).Theta()*TMath::RadToDeg());
        // Get resolutions
        resth = GammaThetaRes->GetBinContent(Ebin, Thbin);
        resph = GammaPhiRes->GetBinContent(Ebin, Thbin);
        resE  = GammaERes->GetBinContent(Ebin, Thbin);
        if(resth==0 || resph==0 || resE==0 ) return; // If energy or angle is out of calibrated  range!
        // Now set particle parameters
        //                     LorentzVector
        photons[i].Set4Vector(meson.SubPhotons(0, i));
        //std::cout << "Res: " << resth << ", " << resph << ", " << photons->Particle(i).E()*resE << std::endl;
        photons[i].SetResolutions(resth, resph, 2 *meson.SubPhotons(0, i).E()*resE);
    }
}










GFit4Constraints::GFit4Constraints(const Bool_t IsEtap) :
    GFit3Constraints(6, 4, 0)
{
}

GFit4Constraints::~GFit4Constraints()
{
}

void    GFit4Constraints::InitFit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
   GFit3Constraints::InitFit(meson);

   Int_t    index[6] = {0,1,2,3,4,5};
   fitter.AddSubMissMassConstraint(beamAndTarget, 6, index, MASS_PROTON);
}

Bool_t    GFit4Constraints::Fit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
    InitFit(meson, beamAndTarget);

    return GFit3Constraints::Fit(meson);
}
