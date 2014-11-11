#include "PPi0Analyse.h"


PPi0Analyse::PPi0Analyse()
{ 
    GHistBGSub::InitCuts(-20, 15, -100, -40);
   GHistBGSub::AddRandCut(35, 95);

  	SetTarget(938); 
        
time_prompt = new TH1F("time_prompt", 	"time_prompt", 	1400, -700, 700);
time1= new TH1F("time1", 	"time1", 	1400, -700, 700);
time_side 	= new TH1F("time_side", 	"time_side", 	1400, -700, 700);
   
IM 		= new TH1F("IM", 	"IM", 		1000,   0, 1000);
MM		= new TH1F("MM", 	"MM", 	 	400,   800, 1200);

IM_proton 		= new TH1F("IM_proton", 	"IM_proton", 		1000,   0, 1000);
MM_proton		= new TH1F("MM_proton", 	"MM_proton", 	 	400,   800, 1200);
theta_proton 	= new TH1F("theta_proton", 	"theta_proton", 		400,   -200, 200);
coplanarity_proton	= new TH1F("coplanarity_proton","coplanarity_proton",400,0,360);

IM_all          = new TH1F("IM_all",         "IM_all",           1000,   0, 1000);
MM_all          = new TH1F("MM_all",         "MM_all",           400,   800, 1200);


IM_all_proton          = new TH1F("IM_all_proton",         "IM_all_proton",           1000,   0, 1000);
MM_all_proton          = new TH1F("MM_all_proton",         "MM_all_proton",           400,   800, 1200);
theta_all_proton       = new TH1F("theta_all_proton",      "theta_all_proton",                400,   -200, 200);
coplanarity_all_proton = new TH1F("coplanarity_all_proton","coplanarity_all_proton",400,0,360);
 
poltable_energy          = new TH1F("poltable_energy",         "poltable_energy",           1400,   0, 1400);
poltable_energy_weight          = new TH1F("poltable_energy_weight",         "poltable_energy_weight",           1400,   0, 1400);

test = new TH1F("test","test",200,-1,1);

//kinematic variables in dependence of energy
IM_energy_kplustplus = new TH3F("IM_energy_kplustplus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kplustminus = new TH3F("IM_energy_kplustminus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustplus = new TH3F("IM_energy_kminustplus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustminus = new TH3F("IM_energy_kminustminus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);

MM_energy_kplustplus = new TH3F("MM_energy_kplustplus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kplustminus = new TH3F("MM_energy_kplustminus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustplus = new TH3F("MM_energy_kminustplus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustminus = new TH3F("MM_energy_kminustminus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);

// //kinematic variables in dependence of energy
IM_energy_kplustplus_proton = new TH3F("IM_energy_kplustplus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kplustminus_proton = new TH3F("IM_energy_kplustminus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustplus_proton = new TH3F("IM_energy_kminustplus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustminus_proton = new TH3F("IM_energy_kminustminus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);

MM_energy_kplustplus_proton = new TH3F("MM_energy_kplustplus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kplustminus_proton = new TH3F("MM_energy_kplustminus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustplus_proton = new TH3F("MM_energy_kminustplus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustminus_proton = new TH3F("MM_energy_kminustminus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);

coplanarityverteilung_collerated = new TH2F("coplanarityverteilung_collerated", "Polar-Verteilung;#phi_{#pi}-#phi_{p}[deg]", 400,0,360,200,200,800);
coplanarityverteilung_collerated->Sumw2();

missingmassverteilung_collerated = new TH2F("missingmassverteilung_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,200,200,800);
missingmassverteilung_collerated->Sumw2();

missingmassverteilung_collerated_proton = new TH2F("missingmassverteilung_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,200,200,800);
missingmassverteilung_collerated_proton->Sumw2();

//kinematic variables in dependence of energy
cosverteilung_collerated = new TH2F("cosverteilung_collerated", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosverteilung_collerated->Sumw2();

// //kinematic variables in dependence of energy
cosverteilung_collerated_proton = new TH2F("cosverteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosverteilung_collerated_proton->Sumw2();
// 
// 
// //Histogramme für 2pi0 System mit Cos-Verteilung___________________________________________________
// 
kristallminus_targetplus_collerated_proton = new TH3F("kristallminus_targetplus_collerated_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_proton->Sumw2();

kristallminus_targetminus_collerated_proton = new TH3F("kristallminus_targetminus_collerated_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_proton->Sumw2();

kristallplus_targetplus_collerated_proton = new TH3F("kristallplus_targetplus_collerated_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_proton->Sumw2();

kristallplus_targetminus_collerated_proton = new TH3F("kristallplus_targetminus_collerated_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_proton->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb_proton = new TH3F("kristallminus_targetplus_collerated_pb_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pb_proton->Sumw2();

kristallminus_targetminus_collerated_pb_proton = new TH3F("kristallminus_targetminus_collerated_pb_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pb_proton->Sumw2();

kristallplus_targetplus_collerated_pb_proton = new TH3F("kristallplus_targetplus_collerated_pb_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pb_proton->Sumw2();

kristallplus_targetminus_collerated_pb_proton = new TH3F("kristallplus_targetminus_collerated_pb_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pb_proton->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt_proton = new TH3F("kristallminus_targetplus_collerated_pt_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pt_proton->Sumw2();

kristallminus_targetminus_collerated_pt_proton = new TH3F("kristallminus_targetminus_collerated_pt_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pt_proton->Sumw2();

kristallplus_targetplus_collerated_pt_proton = new TH3F("kristallplus_targetplus_collerated_pt_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pt_proton->Sumw2();

kristallplus_targetminus_collerated_pt_proton = new TH3F("kristallplus_targetminus_collerated_pt_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pt_proton->Sumw2();

//Histogramme für 2pi0 System mit Cos-Verteilung___________________________________________________

kristallminus_targetplus_collerated = new TH3F("kristallminus_targetplus_collerated","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated->Sumw2();

kristallminus_targetminus_collerated = new TH3F("kristallminus_targetminus_collerated","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated->Sumw2();

kristallplus_targetplus_collerated = new TH3F("kristallplus_targetplus_collerated","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated->Sumw2();

kristallplus_targetminus_collerated = new TH3F("kristallplus_targetminus_collerated","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb = new TH3F("kristallminus_targetplus_collerated_pb","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pb->Sumw2();

kristallminus_targetminus_collerated_pb = new TH3F("kristallminus_targetminus_collerated_pb","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pb->Sumw2();

kristallplus_targetplus_collerated_pb = new TH3F("kristallplus_targetplus_collerated_pb","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pb->Sumw2();

kristallplus_targetminus_collerated_pb = new TH3F("kristallplus_targetminus_collerated_pb","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pb->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt = new TH3F("kristallminus_targetplus_collerated_pt","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pt->Sumw2();

kristallminus_targetminus_collerated_pt = new TH3F("kristallminus_targetminus_collerated_pt","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pt->Sumw2();

kristallplus_targetplus_collerated_pt = new TH3F("kristallplus_targetplus_collerated_pt","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pt->Sumw2();

kristallplus_targetminus_collerated_pt = new TH3F("kristallplus_targetminus_collerated_pt","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pt->Sumw2();

pionmasse = 134.9766;
unten_mass = -67.0;
oben_mass = 67.0;
unten_copl = -25.;
oben_copl = 25.;
unten_inv = -23.;
oben_inv = 23.;
unten_theta=-10.;
oben_theta=10.;
timebackground = -1*35./120.;

}

PPi0Analyse::~PPi0Analyse()
{
}

Bool_t	PPi0Analyse::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();
TDirectory* curDir1  = file_out->mkdir(Form("Butanol%i",linpol->GetEdgeSetting()));
poltable_energy->Write();
poltable_energy_weight->Write();
TDirectory* curDir3  = curDir1->mkdir("Selektion");
curDir3->cd();
time_prompt->Write();
time1->Write();
time_side->Write();
IM_all->Write();
MM_all->Write();  
IM->Write();
MM->Write();
test->Write();
IM_energy_kplustplus->Write();
IM_energy_kplustminus->Write();
IM_energy_kminustplus->Write();
IM_energy_kminustminus->Write();
MM_energy_kplustplus->Write();
MM_energy_kplustminus->Write();
MM_energy_kminustplus->Write();
MM_energy_kminustminus->Write();

curDir1->cd();
TDirectory* curDir5  = curDir1->mkdir("Selektion_withProton");
curDir5->cd();
theta_all_proton->Write();
coplanarity_all_proton->Write();
IM_all_proton->Write();
MM_all_proton->Write();  
IM_proton->Write();
MM_proton->Write();
theta_proton->Write();
coplanarity_proton->Write();
IM_energy_kplustplus_proton->Write();
IM_energy_kplustminus_proton->Write();
IM_energy_kminustplus_proton->Write();
IM_energy_kminustminus_proton->Write();

MM_energy_kplustplus_proton->Write();
MM_energy_kplustminus_proton->Write();
MM_energy_kminustplus_proton->Write();
MM_energy_kminustminus_proton->Write();

curDir1->cd();
TDirectory* curDir2  = curDir1->mkdir("Allgemein");
curDir2->cd();
kristallminus_targetplus_collerated->Write();
kristallminus_targetminus_collerated->Write();
kristallplus_targetplus_collerated->Write();
kristallplus_targetminus_collerated->Write();

kristallminus_targetplus_collerated_pb->Write();
kristallminus_targetminus_collerated_pb->Write();
kristallplus_targetplus_collerated_pb->Write();
kristallplus_targetminus_collerated_pb->Write();

kristallminus_targetplus_collerated_pt->Write();
kristallminus_targetminus_collerated_pt->Write();
kristallplus_targetplus_collerated_pt->Write();
kristallplus_targetminus_collerated_pt->Write();
missingmassverteilung_collerated->Write();
cosverteilung_collerated->Write();
curDir1->cd();
TDirectory* curDir4  = curDir1->mkdir("Allgemein_withProton");
curDir4->cd();
kristallminus_targetplus_collerated_proton->Write();
kristallminus_targetminus_collerated_proton->Write();
kristallplus_targetplus_collerated_proton->Write();
kristallplus_targetminus_collerated_proton->Write();

kristallminus_targetplus_collerated_pb_proton->Write();
kristallminus_targetminus_collerated_pb_proton->Write();
kristallplus_targetplus_collerated_pb_proton->Write();
kristallplus_targetminus_collerated_pb_proton->Write();

kristallminus_targetplus_collerated_pt_proton->Write();
kristallminus_targetminus_collerated_pt_proton->Write();
kristallplus_targetplus_collerated_pt_proton->Write();
kristallplus_targetminus_collerated_pt_proton->Write();
cosverteilung_collerated_proton->Write();
coplanarityverteilung_collerated->Write();
missingmassverteilung_collerated_proton->Write();
cout << file_in->GetName() << endl;
file_out->Close();

	return kTRUE;
}

void	PPi0Analyse::ProcessEvent()	
{
// 			cout << pt << "\t" << planesetting << "\t" << polsetting << endl; 

	
	TLorentzVector pi0_4vect = photons->Particle(0)+photons->Particle(1);
	//invMass
	Double_t inv=pi0_4vect.M();


	test->Fill(pt);


		for(Int_t j=0; j < tagger->GetNTagged();j++)
		{
			Double_t time= tagger->GetTagged_t(j) - 0.5*(photons->GetTime(1)+photons->GetTime(0));
			time1->Fill(time);

			if(GHistBGSub::IsPrompt(time)){
				time_prompt->Fill(time);
			}

				if(GHistBGSub::IsRandom(time)){
					time_side->Fill(time);
				}

			poltable_energy->Fill(tagger->GetPhotonBeam_E(j));
			poltable_energy_weight->Fill(tagger->GetPhotonBeam_E(j),linpol->GetPolDegree(tagger->GetTagged_ch(j)));

			
// 				cout << linpol->GetPolPlane() << endl;
				cout << trigger->GetErrCode() << endl;
		
				TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
				TLorentzVector  beam_4vect = TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(j),tagger->GetPhotonBeam_E(j));
				TLorentzVector  missingp_4vect = beam_4vect + protonvektor_target - pi0_4vect;
				Double_t beamphoton1E=tagger->GetPhotonBeam_E(j);
				Double_t pb=linpol->GetPolDegree(tagger->GetTagged_ch(j));
				//Missing Mass
				Double_t missingmass=missingp_4vect.M();

				//Boost-System
				TLorentzVector pi0_4vect_boost = CMVector(pi0_4vect, beam_4vect, protonvektor_target);
				TVector3 pi0_3vektor_boost = pi0_4vect_boost.Vect();
				Double_t cospi0 = pi0_3vektor_boost.CosTheta();
				Double_t phimeson = 360*(pi0_4vect.Vect().Phi())/(2*TMath::Pi());

				//without any cuts!
				if(GHistBGSub::IsPrompt(time)){
					MM->Fill(missingmass);
					IM->Fill(inv);
				}
	
				if(GHistBGSub::IsRandom(time)){
					MM->Fill(missingmass,timebackground);
					IM->Fill(inv,timebackground);
				}


				//with cuts
				if(GHistBGSub::IsPrompt(time)){
		
					if(((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){

						if(linpol->GetPolPlane()==0){
							if(pt>0 || pt==5){
								MM_energy_kminustplus->Fill(missingmass,cospi0,beamphoton1E);
							}
	
							if(pt<0 || pt==5){
								MM_energy_kminustminus->Fill(missingmass,cospi0,beamphoton1E);
							}
						}

						if(linpol->GetPolPlane()==1){
							if(pt>0 || pt==5){
								MM_energy_kplustplus->Fill(missingmass,cospi0,beamphoton1E);
							}
				
							if(pt<0 || pt==5){
								MM_energy_kplustminus->Fill(missingmass,cospi0,beamphoton1E);
							}
						}

						MM_all->Fill(missingmass);
						missingmassverteilung_collerated->Fill(missingmass,beamphoton1E);

					}
	
					if(((missingmass > (938.+unten_mass) && missingmass < (938.+oben_mass)))){

						if(linpol->GetPolPlane()==0){
							if(pt>0 || pt==5){
								IM_energy_kminustplus->Fill(inv,cospi0,beamphoton1E);
							}
							if(pt<0 || pt==5){
								IM_energy_kminustminus->Fill(inv,cospi0,beamphoton1E);
							}
						}

						if(linpol->GetPolPlane()==1){
							if(pt>0 || pt==5){
								IM_energy_kplustplus->Fill(inv,cospi0,beamphoton1E);
							}
							if(pt<0 || pt==5){
								IM_energy_kplustminus->Fill(inv,cospi0,beamphoton1E);
							}
						}
						IM_all->Fill(inv);
					}
				}

				if(GHistBGSub::IsRandom(time)){
	
					if(((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){

						if(linpol->GetPolPlane()==0){
							if(pt>0 || pt==5){
								MM_energy_kminustplus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
							}
							if(pt<0 || pt==5){
								MM_energy_kminustminus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
							}
						}

						if(linpol->GetPolPlane()==1){
			
							if(pt>0 || pt==5){
								MM_energy_kplustplus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
							}		
							if(pt<0 || pt==5){
								MM_energy_kplustminus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
							}
						}						
						MM_all->Fill(missingmass,timebackground);
						missingmassverteilung_collerated->Fill(missingmass,beamphoton1E,timebackground);

					}
	
					if(((missingmass > (938.+unten_mass) && missingmass < (938.+oben_mass)))){
						if(linpol->GetPolPlane()==0){
							if(pt>0 || pt==5){
								IM_energy_kminustplus->Fill(inv,cospi0,beamphoton1E,timebackground);
							}
							if(pt<0 || pt==5){
								IM_energy_kminustminus->Fill(inv,cospi0,beamphoton1E,timebackground);
							}
						}

						if(linpol->GetPolPlane()==1){
							if(pt>0 || pt==5){
								IM_energy_kplustplus->Fill(inv,cospi0,beamphoton1E,timebackground);
							}
							if(pt<0 || pt==5){
								IM_energy_kplustminus->Fill(inv,cospi0,beamphoton1E,timebackground);
							}
						}
						IM_all->Fill(inv,timebackground);
					}
				}
// 				cout << rootinos->GetNParticles() << endl;
				if(protons->GetNParticles()==1 || electrons->GetNParticles()==1 || chargedPi->GetNParticles()==1 || rootinos->GetNParticles()==1){

					//Phi-Difference
					TLorentzVector proton_4vect_meas;
					if(electrons->GetNParticles()==1){proton_4vect_meas=electrons->Particle(0);}
					if(protons->GetNParticles()==1){proton_4vect_meas=protons->Particle(0);}
					if(chargedPi->GetNParticles()==1){proton_4vect_meas=chargedPi->Particle(0);}
					if(rootinos->GetNParticles()==1){proton_4vect_meas=rootinos->Particle(0);}

					Double_t anglephimeson = pi0_4vect.Vect().Phi();
					Double_t anglephiproton = proton_4vect_meas.Vect().Phi();
					Double_t phi = 360*(anglephimeson - anglephiproton)/(2*TMath::Pi());
					if(phi < 0){phi = phi + 360;}
			
					//theta differenz
					Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
					Double_t anglethetaproton_meas = TMath::RadToDeg()*proton_4vect_meas.Vect().Theta();
					Double_t thetadiff=anglethetaproton_rek-anglethetaproton_meas;

					//without any cuts!
					if(GHistBGSub::IsPrompt(time)){
						MM_proton->Fill(missingmass);
						theta_proton->Fill(thetadiff);
						coplanarity_proton->Fill(phi);
						IM_proton->Fill(inv);
					}
		
					if(GHistBGSub::IsRandom(time)){
						MM_proton->Fill(missingmass,timebackground);
						IM_proton->Fill(inv,timebackground);
						theta_proton->Fill(thetadiff,timebackground);
						coplanarity_proton->Fill(phi,timebackground);
					}


					if(GHistBGSub::IsPrompt(time)){
						if((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl))) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){
							theta_all_proton->Fill(thetadiff);
						}
	
						if((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass)) && (thetadiff > (unten_theta) && thetadiff < (oben_theta)) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){
							coplanarity_all_proton->Fill(phi);
							coplanarityverteilung_collerated->Fill(phi, beamphoton1E);
						}

						if((thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl))) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){
		
							if(linpol->GetPolPlane()==0){
								if(pt>0 || pt==5){
									MM_energy_kminustplus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
		
								if(pt<0 || pt==5){
									MM_energy_kminustminus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
							}	

							if(linpol->GetPolPlane()==1){
								if(pt>0 || pt==5){
									MM_energy_kplustplus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
					
								if(pt<0 || pt==5){
									MM_energy_kplustminus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
							}	
							MM_all_proton->Fill(missingmass);
							missingmassverteilung_collerated_proton->Fill(missingmass,beamphoton1E);

						}
	
						if(((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl)))){
							if(linpol->GetPolPlane()==0){
								if(pt>0 || pt==5){
									IM_energy_kminustplus_proton->Fill(inv,cospi0,beamphoton1E);
								}
								if(pt<0 || pt==5){
									IM_energy_kminustminus_proton->Fill(inv,cospi0,beamphoton1E);
								}
							}

							if(linpol->GetPolPlane()==1){
								if(pt>0 || pt==5){
									IM_energy_kplustplus_proton->Fill(inv,cospi0,beamphoton1E);
								}
								if(pt<0 || pt==5){
									IM_energy_kplustminus_proton->Fill(inv,cospi0,beamphoton1E);
								}
							}							
							IM_all_proton->Fill(inv);
						}

					}

					if(GHistBGSub::IsRandom(time)){
						if((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl))) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){
							theta_all_proton->Fill(thetadiff,timebackground);
						}
	
						if((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass)) && (thetadiff > (unten_theta) && thetadiff < (oben_theta)) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){
							coplanarity_all_proton->Fill(phi,timebackground);
							coplanarityverteilung_collerated->Fill(phi, beamphoton1E,timebackground);
						}

						if((thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl))) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){	
							if(linpol->GetPolPlane()==0){
								if(pt>0 || pt==5){
									MM_energy_kminustplus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
								}
								if(pt<0 || pt==5){
									MM_energy_kminustminus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
								}
							}

							if(linpol->GetPolPlane()==1){
				
								if(pt>0 || pt==5){
									MM_energy_kplustplus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
								}		
								if(pt<0 || pt==5){
									MM_energy_kplustminus->Fill(missingmass,cospi0,beamphoton1E,timebackground);
								}
							}								
							MM_all_proton->Fill(missingmass,timebackground);
							missingmassverteilung_collerated_proton->Fill(missingmass,beamphoton1E,timebackground);

						}
	
						if(((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl)))){
							if(linpol->GetPolPlane()==0){
								if(pt>0 || pt==5){
									IM_energy_kminustplus->Fill(inv,cospi0,beamphoton1E,timebackground);
								}
								if(pt<0 || pt==5){
									IM_energy_kminustminus->Fill(inv,cospi0,beamphoton1E,timebackground);
								}
							}

							if(linpol->GetPolPlane()==1){
								if(pt>0 || pt==5){
									IM_energy_kplustplus->Fill(inv,cospi0,beamphoton1E,timebackground);
								}
								if(pt<0 || pt==5){
									IM_energy_kplustminus->Fill(inv,cospi0,beamphoton1E,timebackground);
								}
							}							
							IM_all_proton->Fill(inv,timebackground);
						}
					}

					if(GHistBGSub::IsPrompt(time)){
	
						if(((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl))) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){
	
						cosverteilung_collerated_proton->Fill(cospi0,beamphoton1E);
					
							if(linpol->GetPolPlane()==0){
			
								if(pt>0 || pt==5){
									kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0, beamphoton1E);
									kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0, beamphoton1E,pb);
									kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0, beamphoton1E,pt);
								}
	
				
								if(pt<0 || pt==5){
									kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E);
									kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,pb);
									kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,pt);
								}
	
							}
	
							if(linpol->GetPolPlane()==1){
			
								if(pt>0 || pt==5){
									kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E);
									kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,pb);
									kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,pt);;
								}
	
				
								if(pt<0 || pt==5){
									kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E);
									kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,pb);
									kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,pt);
								}
	
							}						
	
						}
	
					}
	
					if(GHistBGSub::IsRandom(time)){
	
						if(((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((180.+unten_copl)) && phi < ((180.+oben_copl))) && ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){
	
						cosverteilung_collerated_proton->Fill(cospi0,beamphoton1E,timebackground);
	
					
							if(linpol->GetPolPlane()==0){
			
								if(pt>0 || pt==5){
									kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0, beamphoton1E,timebackground);
									kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0, beamphoton1E,timebackground*pb);
									kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0, beamphoton1E,timebackground*pt);
								}
	
				
								if(pt<0 || pt==5){
									kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground);
									kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground*pb);
									kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground*pt);
								}
	
							}
	
							if(linpol->GetPolPlane()==1){
			
								if(pt>0 || pt==5){
									kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground);
									kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground*pb);
									kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground*pt);;
								}
	
				
								if(pt<0 || pt==5){
									kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground);
									kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground*pb);
									kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,timebackground*pt);
								}
	
							}						
	
						}
	
					}
// 				
				}
// 				//all events with proton identified


				if(GHistBGSub::IsPrompt(time)){

					if((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass))&& ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){

					cosverteilung_collerated->Fill(cospi0,beamphoton1E);
				
						if(linpol->GetPolPlane()==0){
		
							if(pt>0 || pt==5){
								kristallminus_targetplus_collerated->Fill(phimeson,cospi0, beamphoton1E);
								kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0, beamphoton1E,pb);
								kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0, beamphoton1E,pt);
							}

			
							if(pt<0 || pt==5){
								kristallminus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E);
 								kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,pb);
 								kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,pt);
							}

						}

						if(linpol->GetPolPlane()==1){
		
							if(pt>0 || pt==5){
								kristallplus_targetplus_collerated->Fill(phimeson,cospi0,beamphoton1E);
								kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,pb);
								kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,pt);;
							}

			
							if(pt<0 || pt==5){
								kristallplus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E);
 								kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,pb);
 								kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,pt);
							}

						}						

					}

				}

				if(GHistBGSub::IsRandom(time)){

					if((missingmass>(938.+unten_mass) && missingmass<(938.+oben_mass))&& ((inv > (pionmasse + unten_inv) && inv < (pionmasse + oben_inv)))){

					cosverteilung_collerated->Fill(cospi0,beamphoton1E,timebackground);

				
						if(linpol->GetPolPlane()==0){
		
							if(pt>0  || pt==5){
								kristallminus_targetplus_collerated->Fill(phimeson,cospi0, beamphoton1E,timebackground);
								kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0, beamphoton1E,timebackground*pb);
								kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0, beamphoton1E,timebackground*pt);
							}

			
							if(pt<0  || pt==5){
								kristallminus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E,timebackground);
 								kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,timebackground*pb);
 								kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,timebackground*pt);
							}

						}

						if(linpol->GetPolPlane()==1){
		
							if(pt>0  || pt==5){
								kristallplus_targetplus_collerated->Fill(phimeson,cospi0,beamphoton1E,timebackground);
								kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,timebackground*pb);
								kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,timebackground*pt);;
							}

			
							if(pt<0  || pt==5){
								kristallplus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E,timebackground);
 								kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,timebackground*pb);
 								kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,timebackground*pt);
							}

						}						

					}

				}

// 			}//is proton?
		}//tagger
}


void PPi0Analyse::fOnEndProcessing() {
}

void PPi0Analyse::fOnBeforeEventProcessing() {

}


void	PPi0Analyse::ProcessScalerRead()
{
pt=targetpol(file_in);

    //time.ScalerReadCorrection(5);
}
