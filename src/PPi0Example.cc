#include "PPi0Example.h"

PPi0Example::PPi0Example()
{ 
    IM 		= new GH1("IM", 	"IM", 		400,   0, 400);
    MM		= new GH1("MM", 	"MM", 	 	400,   800, 1200);
    IM1 		= new TH1F("IM1", 	"IM1", 		400,   0, 400);
    IM1_side 		= new TH1F("IM1_side", 	"IM1_side", 		400,   0, 400);

    MM1		= new TH1F("MM1", 	"MM1", 	 	400,   800, 1200);

GHistBGSub::InitCuts(-8,8,-300,-100);
GHistBGSub::AddRandCut(100,300);

}

PPi0Example::~PPi0Example()
{
}

Bool_t	PPi0Example::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;

	return kTRUE;
}

Bool_t	PPi0Example::Start()
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

void	PPi0Example::ProcessEvent()
{

	if(GetPhotons()->GetNParticles()==2){
		TLorentzVector pi0_4vect = GetPhotons()->Particle(0)+GetPhotons()->Particle(1);
		//invMass
		Double_t inv=pi0_4vect.M();
	
		for(Int_t j=0; j < GetTagger()->GetNTagged();j++)
		{

			if(GetPhotons()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;
			if(GetPhotons()->GetDetectors(1)==GTreeTrack::DETECTOR_PbWO4) continue;
			if(GetRootinos()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;

			Double_t time= GetTagger()->GetTaggedTime(j) - 0.5*(GetPhotons()->GetTime(1)+GetPhotons()->GetTime(0));	

			//get target, beam and missng particle information
			TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
			TLorentzVector  beam_4vect = TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(j),GetTagger()->GetTaggedEnergy(j));
			TLorentzVector  missingp_4vect = beam_4vect + protonvektor_target - pi0_4vect;
			Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
			Double_t beamphoton1E=GetTagger()->GetTaggedEnergy(j);

			//due to kinematic not possible -> kick them out
			if(anglethetaproton_rek > 90){continue;}

			//Missing Mass
			Double_t missingmass=missingp_4vect.M();
			if(time > -8 && time < 8){
				IM1->Fill(inv,0.5);

			}

			if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
				IM1->Fill(inv,-0.5*0.04);

				IM1_side->Fill(inv,0.04*0.5);
			}
			IM->Fillweighted(inv,time,0.5); //WIE MACHE ICH MIT GH1 WICHTUNGEN?????
		}
	}



}

void	PPi0Example::ProcessScalerRead()
{
	// Fill Tagger Scalers
}

Bool_t	PPi0Example::Write()
{

	IM1->Write();
    // Write all GH1's and TObjects defined in this class
    return GTreeManager::Write();
}
