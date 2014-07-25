#include "GH1Example.h"



GH1Example::GH1Example()    :
    hist_eta("eta", "eta", "eta"),
    check_eta_proton("check_eta_proton", "check_eta_proton", "eta"),
    hist_etap("etap", "etap", "etap"),
    check_etap_proton("check_etap_proton", "check_etap_proton", "etap")
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
    {
        if(protons->GetNParticles()>0)
            check_eta_proton.Check(*eta, *protons, *tagger);
        else
            hist_eta.Fill(*eta, *tagger, kTRUE);
    }
    if(etap->GetNParticles()>0)
    {
        if(protons->GetNParticles()>0)
            check_etap_proton.Check(*etap, *protons, *tagger);
        else
            hist_etap.Fill(*etap, *tagger, kTRUE);
    }
}

void	GH1Example::ProcessScalerRead()
{
    //hist_eta.ScalerReadCorrection(1/0.65);
    //hist_etap.ScalerReadCorrection(1/0.65);
}
