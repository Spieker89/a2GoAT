#include "GH1Example.h"



GH1Example::GH1Example()    :
    test("test", "test", "test")
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
        test.Fill(eta->Particle(0), *tagger, kTRUE);
}

void	GH1Example::ProcessScalerRead()
{
    //test.ScalerReadCorrection(5);
}
