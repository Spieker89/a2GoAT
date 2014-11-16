#ifndef __MyPhysics_h__
#define __MyPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GHistEvent.h"
#include "GFit.h"

class	MyPhysics  : public GTreeManager
{
private:
    GHistEvent3Mesons   hist_raw;
    GHistEvent3Mesons   hist_SubImCut;
    GHistEvent3Mesons   hist_MMCut;

    GFit3Constraints    fit3;
    GFit4Constraints    fit4;
    GFit3ConstraintsBeam    fit3Beam;
    GFit4ConstraintsBeam    fit4Beam;
    GFit4ConstraintsProton  fit4Proton;
    GFit3ConstraintsBeamProton    fit3BeamProton;
    GFit4ConstraintsBeamProton    fit4BeamProton;

    GH1                 fit3_im;
    GH1                 fit3_cs;
    GH1                 fit3_cl;
    GH1                 fit4_im;
    GH1                 fit4_cs;
    GH1                 fit4_cl;
    GH1                 fit3Beam_im;
    GH1                 fit3Beam_cs;
    GH1                 fit3Beam_cl;
    GH1                 fit4Beam_im;
    GH1                 fit4Beam_cs;
    GH1                 fit4Beam_cl;
    GH1                 fit4Proton_im;
    GH1                 fit4Proton_cs;
    GH1                 fit4Proton_cl;
    GH1                 fit3BeamProton_im;
    GH1                 fit3BeamProton_cs;
    GH1                 fit3BeamProton_cl;
    GH1                 fit4BeamProton_im;
    GH1                 fit4BeamProton_cs;
    GH1                 fit4BeamProton_cl;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    MyPhysics();
    virtual ~MyPhysics();

    virtual Bool_t	Init(const char* configfile);

};
#endif
