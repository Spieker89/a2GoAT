#include "GFit.h"
#include "GTreeMeson.h"


GFit3Constraints::GFit3Constraints(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    fitter(6, 3, 0)
{

}

GFit3Constraints::~GFit3Constraints()
{

}

void    GFit3Constraints::Set(const TLorentzVector& p0,
                              const TLorentzVector& p1,
                              const TLorentzVector& p2,
                              const TLorentzVector& p3,
                              const TLorentzVector& p4,
                              const TLorentzVector& p5)
{
    fitter.Reset();

    GKinFitterParticle  photons[6];

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
}














GFit4Constraints::GFit4Constraints(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    fitter(6, 4, 0)
{

}

GFit4Constraints::~GFit4Constraints()
{

}

void    GFit4Constraints::Set(const TLorentzVector& p0,
                              const TLorentzVector& p1,
                              const TLorentzVector& p2,
                              const TLorentzVector& p3,
                              const TLorentzVector& p4,
                              const TLorentzVector& p5,
                              const TLorentzVector& beamAndTarget)
{
    fitter.Reset();

    GKinFitterParticle  photons[6];

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddSubMissMassConstraint(beamAndTarget, 6, &index[0], MASS_PROTON);
}



















GFit4ConstraintsBeam::GFit4ConstraintsBeam(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    fitter(7, 4, 0)
{

}

GFit4ConstraintsBeam::~GFit4ConstraintsBeam()
{

}

TLorentzVector  GFit4ConstraintsBeam::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit4ConstraintsBeam::Set(const TLorentzVector& p0,
                              const TLorentzVector& p1,
                              const TLorentzVector& p2,
                              const TLorentzVector& p3,
                              const TLorentzVector& p4,
                              const TLorentzVector& p5,
                              const TLorentzVector& beamAndTarget)
{
    fitter.Reset();

    GKinFitterParticle  photons[6];
    GKinFitterParticle  beam;

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }
    beam.Set4Vector(beamAndTarget);
    beam.SetResolutions(1, 1, 2);
    fitter.AddNegKFParticle(beam);

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddInvMassConstraint(MASS_PROTON);
}














GFit4ConstraintsProton::GFit4ConstraintsProton(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    fitter(7, 4, 0)
{

}

GFit4ConstraintsProton::~GFit4ConstraintsProton()
{

}

TLorentzVector  GFit4ConstraintsProton::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit4ConstraintsProton::Set(const TLorentzVector& p0,
                                    const TLorentzVector& p1,
                                    const TLorentzVector& p2,
                                    const TLorentzVector& p3,
                                    const TLorentzVector& p4,
                                    const TLorentzVector& p5,
                                    const TLorentzVector& beamAndTarget,
                                    const TLorentzVector& proton)
{
    fitter.Reset();

    GKinFitterParticle  photons[6];
    GKinFitterParticle  pro;

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }
    TLorentzVector  help(proton);
    help.SetE(beamAndTarget.E()-p0.E()-p1.E()-p2.E()-p3.E()-p4.E()-p5.E());
    help.SetVect(proton.Vect().Unit()*(TMath::Sqrt(help.E()*help.E() - MASS_PROTON*MASS_PROTON)));
    pro.Set4Vector(help);
    pro.SetResolutions(3, 3, 0.1);
    fitter.AddPosKFParticle(pro);

    Int_t   index[7]    = {0, 1, 2, 3, 4, 5, 6};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddSubMissMassConstraint(beamAndTarget, 7, &index[0], 0);
}
























GFit4ConstraintsBeamProton::GFit4ConstraintsBeamProton(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    fitter(8, 5, 0)
{

}

GFit4ConstraintsBeamProton::~GFit4ConstraintsBeamProton()
{

}

TLorentzVector  GFit4ConstraintsBeamProton::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit4ConstraintsBeamProton::Set(const TLorentzVector& p0,
                                        const TLorentzVector& p1,
                                        const TLorentzVector& p2,
                                        const TLorentzVector& p3,
                                        const TLorentzVector& p4,
                                        const TLorentzVector& p5,
                                        const TLorentzVector& beamAndTarget,
                                        const TLorentzVector& proton)
{
    fitter.Reset();

    GKinFitterParticle  photons[6];
    GKinFitterParticle  pro;
    GKinFitterParticle  beam;

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    };
    beam.Set4Vector(beamAndTarget);
    beam.SetResolutions(1, 1, 2);
    fitter.AddNegKFParticle(beam);
    TLorentzVector  help(proton);
    help.SetE(beamAndTarget.E()-p0.E()-p1.E()-p2.E()-p3.E()-p4.E()-p5.E());
    help.SetVect(proton.Vect().Unit()*(TMath::Sqrt(help.E()*help.E() - MASS_PROTON*MASS_PROTON)));
    //help.Print();
    //std::cout << help.M() << std::endl;
    pro.Set4Vector(help);
    pro.SetResolutions(3, 3, 0.1);
    fitter.AddPosKFParticle(pro);

    Int_t   index[8]    = {0, 1, 2, 3, 4, 5, 6, 7};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    //fitter.AddSubInvMassConstraint(1, &index[7], MASS_PROTON);
    fitter.AddTotEnergyConstraint(0);
    //fitter.AddTotMomentumConstraint(TVector3(0, 0, 0));
    fitter.AddInvMassConstraint(0);
}












GHistFit::GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    nPulls(_NPulls),
    im(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, 48, kFALSE),
    chiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE),
    pulls(TString(name).Append("_Pulls"), TString(title).Append(" Pulls"), 100, -5, 5, nPulls, 0, nPulls, kFALSE)
{

}

GHistFit::~GHistFit()
{

}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime)
{
    im.Fill(fitter.GetTotalFitParticle().M(), taggerTime);
    chiSq.Fill(fitter.GetChi2(), taggerTime);
    confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime);
    for(int i=0; i<nPulls; i++)
        pulls.Fill(fitter.GetPull(i), i);
}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel)
{
    im.Fill(fitter.GetTotalFitParticle().M(), taggerTime, taggerChannel);
    chiSq.Fill(fitter.GetChi2(), taggerTime);
    confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime);
    for(int i=0; i<nPulls; i++)
        pulls.Fill(fitter.GetPull(i), i);
}

void    GHistFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        chiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        pulls.PrepareWriteList(arr, TString(name).Append("_Pulls").Data());
    }
    else
    {
        im.PrepareWriteList(arr);
        chiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        pulls.PrepareWriteList(arr);
    }
}
