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













GFit3ConstraintsBeam::GFit3ConstraintsBeam(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    fitter(7, 3, 0)
{

}

GFit3ConstraintsBeam::~GFit3ConstraintsBeam()
{

}

TLorentzVector  GFit3ConstraintsBeam::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit3ConstraintsBeam::Set(const TLorentzVector& p0,
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
    fitter(7, 5, 0)
{

}

GFit4ConstraintsProton::~GFit4ConstraintsProton()
{

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
    fitter.AddSubInvMassConstraint(1, &index[6], MASS_PROTON);
    fitter.AddSubMissMassConstraint(beamAndTarget, 6, &index[0], MASS_PROTON);
}













GFit3ConstraintsBeamProton::GFit3ConstraintsBeamProton(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    fitter(8, 4, 0)
{

}

GFit3ConstraintsBeamProton::~GFit3ConstraintsBeamProton()
{

}

TLorentzVector  GFit3ConstraintsBeamProton::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit3ConstraintsBeamProton::Set(const TLorentzVector& p0,
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
    }
    beam.Set4Vector(beamAndTarget);
    beam.SetResolutions(1, 1, 2);
    fitter.AddNegKFParticle(beam);
    TLorentzVector  help(proton);
    help.SetE(beamAndTarget.E()-p0.E()-p1.E()-p2.E()-p3.E()-p4.E()-p5.E());
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
    fitter.AddSubInvMassConstraint(1, &index[7], MASS_PROTON);
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
    fitter.AddSubInvMassConstraint(1, &index[7], MASS_PROTON);
    fitter.AddInvMassConstraint(0);
}
