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

Bool_t  GFit3Constraints::InitFit(const GTreeMeson& meson)
{
    if(SetPhotons(meson)==kFALSE)
        return kFALSE;

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

    return kTRUE;
}

Bool_t  GFit3Constraints::Fit(const GTreeMeson& meson)
{
    if(InitFit(meson)==kFALSE)
        return kFALSE;

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

Bool_t  GFit3Constraints::SetPhotons(const GTreeMeson& meson)
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
        if(resth==0 || resph==0 || resE==0 ) return kFALSE; // If energy or angle is out of calibrated  range!
        // Now set particle parameters
        //                     LorentzVector
        photons[i].Set4Vector(meson.SubPhotons(0, i));
        //std::cout << "Res: " << resth << ", " << resph << ", " << photons->Particle(i).E()*resE << std::endl;
        photons[i].SetResolutions(resth, resph, 2 *meson.SubPhotons(0, i).E()*resE);
    }
    return kTRUE;
}










GFit4Constraints::GFit4Constraints(const Int_t npart, const Int_t ncon, const Bool_t IsEtap) :
    GFit3Constraints(npart, ncon, 0)
{
}

GFit4Constraints::GFit4Constraints(const Bool_t IsEtap) :
    GFit3Constraints(6, 4, 0)
{
}

GFit4Constraints::~GFit4Constraints()
{
}

Bool_t  GFit4Constraints::InitFit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
   if(GFit3Constraints::InitFit(meson)==kFALSE)
       return kFALSE;

   Int_t    index[6] = {0,1,2,3,4,5};
   fitter.AddSubMissMassConstraint(beamAndTarget, 6, index, MASS_PROTON);

   return kTRUE;
}

Bool_t  GFit4Constraints::Fit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
    if(InitFit(meson, beamAndTarget)==kFALSE)
        return kFALSE;

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







GFit3ConstraintsBeam::GFit3ConstraintsBeam(const Bool_t IsEtap) :
    GFit3Constraints(7, 3, 0)
{
}

GFit3ConstraintsBeam::~GFit3ConstraintsBeam()
{
}

Bool_t  GFit3ConstraintsBeam::InitFit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
   if(GFit3Constraints::InitFit(meson)==kFALSE)
       return kFALSE;

   //printf("%lf, %lf, %lf, %lf, %lf\n", (beamAndTarget-meson.Particle(0)).Px(), (beamAndTarget-meson.Particle(0)).Py(), (beamAndTarget-meson.Particle(0)).Pz(), (beamAndTarget-meson.Particle(0)).E(), (beamAndTarget-meson.Particle(0)).M());
   GKinFitterParticle   initial(beamAndTarget, 1, 1, 1);
   fitter.AddNegKFParticle(initial);

   return kTRUE;
}

Bool_t  GFit3ConstraintsBeam::Fit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
    if(InitFit(meson, beamAndTarget)==kFALSE)
        return kFALSE;

    if(fitter.Solve()>=0)
    {
        result.im               = (fitter.GetParticle(0).Get4Vector()+fitter.GetParticle(1).Get4Vector()+fitter.GetParticle(2).Get4Vector()+fitter.GetParticle(3).Get4Vector()+fitter.GetParticle(4).Get4Vector()+fitter.GetParticle(5).Get4Vector()).M();
        result.ChiSq            = fitter.GetChi2();
        result.ConfidenceLevel  = fitter.ConfidenceLevel();
        for(int i=0; i<24; i++)
            result.PullPhotons[i]   = fitter.Pull(i);
        for(int i=0; i<4; i++)
            result.PullBeam[i]   = fitter.Pull(i+24);
        return kTRUE;
    }
    return kFALSE;
}






GFit4ConstraintsBeam::GFit4ConstraintsBeam(const Bool_t IsEtap) :
    GFit4Constraints(7, 4, 0)
{
}

GFit4ConstraintsBeam::~GFit4ConstraintsBeam()
{
}

Bool_t  GFit4ConstraintsBeam::InitFit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
   if(GFit4Constraints::InitFit(meson, beamAndTarget)==kFALSE)
       return kFALSE;

   GKinFitterParticle   initial(beamAndTarget, 1, 1, 1);
   fitter.AddNegKFParticle(initial);

   return kTRUE;
}

Bool_t  GFit4ConstraintsBeam::Fit(const GTreeMeson& meson, const TLorentzVector& beamAndTarget)
{
    if(InitFit(meson, beamAndTarget)==kFALSE)
        return kFALSE;

    if(fitter.Solve()>=0)
    {
        result.im               = (fitter.GetParticle(0).Get4Vector()+fitter.GetParticle(1).Get4Vector()+fitter.GetParticle(2).Get4Vector()+fitter.GetParticle(3).Get4Vector()+fitter.GetParticle(4).Get4Vector()+fitter.GetParticle(5).Get4Vector()).M();
        result.ChiSq            = fitter.GetChi2();
        result.ConfidenceLevel  = fitter.ConfidenceLevel();
        for(int i=0; i<24; i++)
            result.PullPhotons[i]   = fitter.Pull(i);
        for(int i=0; i<4; i++)
            result.PullBeam[i]   = fitter.Pull(i+24);
        return kTRUE;
    }
    return kFALSE;
}
