#ifndef __GFitBeamProton_h__
#define __GFitBeamProton_h__


#include "GFitBeam.h"



class	GFitBeamProton    : public GFitBeam
{
private:
    FitParticle     aplconProton;

    GH1                 protonEnergy;
    GH1                 protonTheta;
    GHistBGSub          protonPhi;

    virtual int GetNDOF()   {return 24;}

public:
    GFitBeamProton(const char* _Name, const Bool_t linkHistogram);
    virtual ~GFitBeamProton()               {}

    virtual void    AddConstraintsTotMomentum();
    virtual void    AddConstraintsTotEnergy();
    virtual void    CalcResult();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    SetProton(const TLorentzVector proton)    {aplconProton.SetFromVector(proton); aplconProton.SetToProtonResolution();}
    virtual bool    Solve(const double time, const int channel);
};





class	GFitBeamProtonVertex    : public GFitBeamProton
{
private:
    GH1     vertex;

    void    ConvertTheta(std::vector<double>& p, std::vector<double>& v)  {p[1] = std::atan2( 25.4*sin(p[1]), 25.4*cos(p[1]) - v[0]);}

    virtual int GetNDOF()   {return 24;}

public:
    GFitBeamProtonVertex(const char* _Name, const Bool_t linkHistogram);
    ~GFitBeamProtonVertex()             {}

    virtual void    AddConstraintsIM();
    virtual void    AddConstraintMM();
    virtual void    AddConstraintsTotMomentum();
    virtual void    AddConstraintsTotEnergy();
    virtual void    CalcResult();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual bool    Solve(const double time, const int channel);
};


#endif
