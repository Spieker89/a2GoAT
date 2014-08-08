#include "MyPhysics.h"



MyPhysics::MyPhysics()    :
    hist_eta("eta", "eta", "eta", kFALSE),
    hist_etap("etap", "etap", "etap", kTRUE)
{ 
        GHistBGSub::InitCuts(-20, 20, -55, -35);
        GHistBGSub::AddRandCut(35, 55);
}

MyPhysics::~MyPhysics()
{

}

Bool_t	MyPhysics::Start()
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

void	MyPhysics::ProcessEvent()
{
    if(eta->GetNParticles()>0)
        hist_eta.Fill(*eta, *protons, *tagger, kTRUE);
    if(etap->GetNParticles()>0)
        hist_etap.Fill(*etap, *protons, *tagger, kTRUE);
}

void	MyPhysics::ProcessScalerRead()
{
    //hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    //hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}
