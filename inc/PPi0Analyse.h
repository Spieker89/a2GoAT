#ifndef __PPi0Analyse_h__
#define __PPi0Analyse_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GTreeLinPol.h"
#include "GTreeTrigger.h"
#include "PPhysics.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCutG.h"
#include <sstream>


class	PPi0Analyse  : public PPhysics
{
protected:
    TH1F*	time1;
    TH1F*	time_prompt;
    TH1F*	time_side;

    TH1F*        IM_all;
    TH1F*        MM_all;  
    TH1F*        IM;
    TH1F*        MM;

    TH1F*        theta_all_proton;
    TH1F*        coplanarity_all_proton;
    TH1F*        IM_all_proton;
    TH1F*        MM_all_proton;  

    TH1F*        IM_proton;
    TH1F*        MM_proton;
     TH1F*        theta_proton;
     TH1F*        coplanarity_proton;
   TCutG *cuttaps = new TCutG("CUT_TAPS_dE_E",10);
   TCutG *cutcb = new TCutG("CUT_CB_dE_E",12);

Int_t kTAPS;
Int_t kCB;
Double_t cutPromptMin,cutPromptMax,cutSideMax,cutSideMin,backgroundSubstractionFactor;
Double_t unten_mass, oben_mass, unten_copl, oben_copl, oben_inv, unten_inv, oben_theta, unten_theta;
Double_t timebackground;
Double_t pionmasse;
Double_t pt;
Double_t pt1;
Int_t runNumber;
TNamed filenumber;
TString* path11;

//kinematic variables in dependence of energy
TH3F *IM_energy_kplustplus;
TH3F *IM_energy_kplustminus;
TH3F *IM_energy_kminustplus;
TH3F *IM_energy_kminustminus;

TH3F *MM_energy_kplustplus;
TH3F *MM_energy_kplustminus;
TH3F *MM_energy_kminustplus;
TH3F *MM_energy_kminustminus;

TH2F *invmassverteilung_collerated;
TH2F *invmassverteilung_collerated_proton;

// //kinematic variables in dependence of energy
TH3F *IM_energy_kplustplus_proton;
TH3F *IM_energy_kplustminus_proton;
TH3F *IM_energy_kminustplus_proton;
TH3F *IM_energy_kminustminus_proton;

TH3F *MM_energy_kplustplus_proton;
TH3F *MM_energy_kplustminus_proton;
TH3F *MM_energy_kminustplus_proton;
TH3F *MM_energy_kminustminus_proton;

TH3F *events_all;
TH3F *events_witherror;
TH2F *Check_CBdE_E;
TH2F *Check_CBdE_E_nocuts;
TH2F *Check_TAPSdE_E;
TH2F *Check_TAPSdE_E_nocuts;
TH2F *Check_TAPS_TOF_proton;
TH2F *Check_TAPS_TOF_photon_3ped;
TH2F *Check_TAPS_TOF_photon_2_3ped;


string polsetting;
string planesetting;
TH1F* test;

TH2F *cosverteilung_collerated;

TH2F *coplanarityverteilung_collerated;
TH2F *thetaverteilung_collerated;
TH2F *missingmassverteilung_collerated;
TH2F *missingmassverteilung_collerated_proton;
TH3F *kristallminus_targetplus_collerated;
TH3F *kristallminus_targetminus_collerated;
TH3F *kristallplus_targetplus_collerated;
TH3F *kristallplus_targetminus_collerated;

TH3F *kristallminus_targetplus_collerated_pb;
TH3F *kristallminus_targetminus_collerated_pb;
TH3F *kristallplus_targetplus_collerated_pb;
TH3F *kristallplus_targetminus_collerated_pb;

TH3F *kristallminus_targetplus_collerated_pt;
TH3F *kristallminus_targetminus_collerated_pt;
TH3F *kristallplus_targetplus_collerated_pt;
TH3F *kristallplus_targetminus_collerated_pt;

TH2F *cosverteilung_collerated_proton;

TH3F *kristallminus_targetplus_collerated_proton;
TH3F *kristallminus_targetminus_collerated_proton;
TH3F *kristallplus_targetplus_collerated_proton;
TH3F *kristallplus_targetminus_collerated_proton;

TH3F *kristallminus_targetplus_collerated_pb_proton;
TH3F *kristallminus_targetminus_collerated_pb_proton;
TH3F *kristallplus_targetplus_collerated_pb_proton;
TH3F *kristallplus_targetminus_collerated_pb_proton;

TH3F *kristallminus_targetplus_collerated_pt_proton;
TH3F *kristallminus_targetminus_collerated_pt_proton;
TH3F *kristallplus_targetplus_collerated_pt_proton;
TH3F *kristallplus_targetminus_collerated_pt_proton;

TH2F *cos_beam_hel1_targetplus;
TH2F *cos_beam_hel0_targetplus;
TH2F *cos_beam_hel1_targetminus;
TH2F *cos_beam_hel0_targetminus;

TH2F *cos_beam_hel1_targetplus_pc;
TH2F *cos_beam_hel0_targetplus_pc;
TH2F *cos_beam_hel1_targetminus_pc;
TH2F *cos_beam_hel0_targetminus_pc;

TH2F *cos_beam_hel1_targetplus_pt;
TH2F *cos_beam_hel0_targetplus_pt;
TH2F *cos_beam_hel1_targetminus_pt;
TH2F *cos_beam_hel0_targetminus_pt;

TH2F *cos_beam_hel1_targetplus_proton;
TH2F *cos_beam_hel0_targetplus_proton;
TH2F *cos_beam_hel1_targetminus_proton;
TH2F *cos_beam_hel0_targetminus_proton;

TH2F *cos_beam_hel1_targetplus_proton_pc;
TH2F *cos_beam_hel0_targetplus_proton_pc;
TH2F *cos_beam_hel1_targetminus_proton_pc;
TH2F *cos_beam_hel0_targetminus_proton_pc;

TH2F *cos_beam_hel1_targetplus_proton_pt;
TH2F *cos_beam_hel0_targetplus_proton_pt;
TH2F *cos_beam_hel1_targetminus_proton_pt;
TH2F *cos_beam_hel0_targetminus_proton_pt;

TH1F *triggertest;
TH1F *poltable_energy;
TH1F *poltable_energy_weight;

TH2F *openingangle_ptopi0_energy_allcuts;
TH2F *openingangle_gammatogamma_energy_allcuts;
TH2F *openingangle_gammatogamma_energy_allcuts_proton;
TH2F *openingangle_ptopi0_energy_nocuts;
TH2F *openingangle_gammatogamma_energy_nocuts;

TH3F *coplanarityverteilung_cospi0_collerated;
TH3F *thetaverteilung_cospi0_collerated;
TH3F *missingmassverteilung_cospi0_collerated;
TH3F *missingmassverteilung_cospi0_collerated_proton;
TH3F *invmassverteilung_cospi0_collerated;
TH3F *invmassverteilung_cospi0_collerated_proton;
TH1F *h_energy_sum_pi0_3ped;
TH1F *h_energy_sum_pi0_2_3ped;
TH1F *beamenergy_gen;
TH2F *cosmeson_beamenergy_monte;
TH2F *cosmeson_beamenergy_rek_3ped;
TH2F *cosmeson_beamenergy_rek_2_3ped;

TH3F *cosmeson_beamenergy_energysum_rek;
TH3F *cosmeson_beamenergy_energysum_rek_proton;
TH3F *cosmeson_beamenergy_energysum_monte;

TH3F *clustersize_cospi0_energy_collerated;
TH3F *clustersize_cospi0_energy_collerated_proton;
TH3F *thetaproton_cospi0_energy_collerated;

// 	Double_t test1;
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual void fOnEndProcessing();
    virtual void fOnBeforeEventProcessing();
    virtual Bool_t    Write();
//     virtual Double_t targetpol(TFile *f);


Double_t GetCircPol(Double_t bph_E, Int_t runnumber1){

	Double_t E_0=1557.0;
	Double_t p_e=0.77;
	Double_t x=bph_E/E_0;

	//for November 2013 beam time
	if(runnumber1 > 200 && runnumber1< 380) p_e=0.7582;
	if(runnumber1 > 379 && runnumber1< 500) p_e=0.7497;
	if(runnumber1 > 499 && runnumber1< 610) p_e=0.7649;
	if(runnumber1 > 609 && runnumber1< 665) p_e=0.7564;
	if(runnumber1 > 664 && runnumber1< 715) p_e=0.7563;
	if(runnumber1 > 714 && runnumber1< 775) p_e=0.7453;
	if(runnumber1 > 774 && runnumber1< 850) p_e=0.7485;
	if(runnumber1 > 849 && runnumber1< 1000) p_e=0.7689;
	if(runnumber1 > 999 && runnumber1< 1200) p_e=0.7856;
	if(runnumber1 > 1199 && runnumber1< 1500) p_e=0.7858;

	//for May 2014 beam time
	if(runnumber1 > 4180 && runnumber1< 4200) p_e=0.7901;
	if(runnumber1 > 4199 && runnumber1< 4244) p_e=0.7864;
	if(runnumber1 > 4243 && runnumber1< 4280) p_e=0.7874;
	if(runnumber1 > 4279 && runnumber1< 4310) p_e=0.7846;
	if(runnumber1 > 4309 && runnumber1< 4359) p_e=0.7870;
	if(runnumber1 > 4358 && runnumber1< 4384) p_e=0.7772;
	if(runnumber1 > 4383 && runnumber1< 4509) p_e=0.7423;

	
	//Double_t pol_deg = x * (3+1-x)/(3-2*(1-x)+3*(1-x)*(1-x)) * p_e;
	Double_t pol_deg = (4*x-x*x)/(4-4*x+3*x*x) * p_e;

	return pol_deg;
}

Double_t targetpol(TFile *f){

	TString* filename = new TString(f->GetPath());
	TString* path1 = new TString(filename->Tokenize("_")->At(filename->Tokenize("_")->GetEntries()-1)->GetName());
	path1->Resize(path1->Length()-7);
	stringstream ss(path1->Data());
	Int_t runnumber1;
	ss >> runnumber1;
	string plane;
	string poledge;

	ifstream in2;
  	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
//	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Diamond_May14_linpol.txt");
// 	in2.open("/hadron/spieker/Mainz_Analyse/Carbon_Diamond_Apr14_linpol.txt");


	Double_t runnumber;
	Double_t pol;
	Double_t pol1;
	Int_t N = 0;

	if (in2.is_open()) {
		while (1) {
		
		in2 >> runnumber >> plane >> poledge >> pol;
	//  	cout << runnumber << endl;
		if (!in2.good()){pol1=500.;break;}
		if(runnumber==runnumber1){pol1=pol;break;}
		
		}
		N++;
	}
	else{cout << "kacke: Datei oeffnet sich nicht!" << endl;}
return pol1/100;
}



string poledge(TFile *f){

	TString* filename = new TString(f->GetPath());
	TString* path1 = new TString(filename->Tokenize("_")->At(filename->Tokenize("_")->GetEntries()-1)->GetName());
	path1->Resize(path1->Length()-7);
	stringstream ss(path1->Data());
	Int_t runnumber1;
	ss >> runnumber1;
	ifstream in2;
// 	in2.open("/hadron/spieker/Mainz_Analyse/Carbon_Diamond_Apr14_linpol.txt");
  	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
//	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Diamond_May14_linpol.txt");



	Double_t runnumber;
	string poledge;
	string poledge1;
	Int_t N = 0;
	string plane;
	Double_t pol;
	if (in2.is_open()) {
		while (1) {
		
		in2 >> runnumber >> plane >> poledge >> pol;
	//  	cout << runnumber << endl;
		if (!in2.good()){poledge1="nichtzuordbar";break;}
		if(runnumber==runnumber1){poledge1=poledge;break;}
		
		}
		N++;
	}
	else{cout << "kacke: Datei oeffnet sich nicht!" << endl;}

return poledge1;
}

string polplane(TFile *f){

	TString* filename = new TString(f->GetPath());
	TString* path1 = new TString(filename->Tokenize("_")->At(filename->Tokenize("_")->GetEntries()-1)->GetName());
	path1->Resize(path1->Length()-7);
	stringstream ss(path1->Data());
	Int_t runnumber1;
	ss >> runnumber1;
	ifstream in2;
// 	in2.open("/hadron/spieker/Mainz_Analyse/Carbon_Diamond_Apr14_linpol.txt");
  	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
//	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Diamond_May14_linpol.txt");


	Double_t runnumber;
	string poledge;
	string poledge1;
	Int_t N = 0;
	string plane;
	string plane1;
	Double_t pol;

	if (in2.is_open()) {
		while (1) {
		
		in2 >> runnumber >> plane >> poledge >> pol;
// 	  	cout << runnumber << "\t" << poledge << "\t" << plane << endl;
		if (!in2.good()){plane1="nichtzuordbar";break;}
		if(runnumber==runnumber1){plane1=plane;break;}
		
		}
		N++;
	}
	else{cout << "kacke: Datei oeffnet sich nicht!" << endl;}

return plane1;
}

TLorentzVector CMVector(TLorentzVector vec,TLorentzVector vec_t,TLorentzVector vec_b)
{
TLorentzVector vec_cm=vec;
TVector3 beta=((vec_t+vec_b).BoostVector())*(-1.0);
vec_cm.Boost(beta);
return vec_cm;
}

			
public:
    PPi0Analyse();
    virtual ~PPi0Analyse();

    //virtual Bool_t	Init(const char* configfile);

};
#endif

