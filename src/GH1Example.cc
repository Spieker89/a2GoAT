#include "GH1Example.h"



GH1Example::GH1Example()    :
    hist_eta("eta", "eta", "eta"),
    hist_etap("etap", "etap", "etap")
{ 
        GHistBGSub::InitCuts(-20, 20, -55, -35);
        GHistBGSub::AddRandCut(35, 55);
}

GH1Example::~GH1Example()
{

}

Bool_t	GH1Example::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();


    TraverseValidEvents();


	return kTRUE;
}

void	GH1Example::ProcessEvent()
{
    if(eta->GetNParticles()>0)
        hist_eta.Fill(eta->Particle(0), *tagger, kTRUE);
    if(etap->GetNParticles()>0)
        hist_etap.Fill(etap->Particle(0), *tagger, kTRUE);
}

void	GH1Example::ProcessScalerRead()
{
    hist_eta.ScalerReadCorrection(1/0.65);
    hist_etap.ScalerReadCorrection(1/0.65);
}
