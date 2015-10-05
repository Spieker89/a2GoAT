#ifndef __NPiPlusAnalyse_h__
#define __NPiPlusAnalyse_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GTreeLinPol.h"
#include "GTreeTrigger.h"
#include "PPhysics.h"
#include "HistoManu.h"
#include "HistoManu2.h"
#include "HistoManu3.h"
#include "TCutG.h"
#include <sstream>
#include "HistoManu.h"
#include "HistoManu2.h"
#include "HistoManu3.h"

class	NPiPlusAnalyse  : public PPhysics
{
protected:
    TH1F*	time1;
    TH1F*	time_prompt;
    TH1F*	time_side;


   TCutG *cuttaps = new TCutG("CUT_TAPS_dE_E",10);
   TCutG *cutcb = new TCutG("CUT_CB_dE_E",12);

Int_t kTAPS;
Int_t kCB_all;
Int_t kCB_CBandPID;
Int_t kCB_CBandMWPC;
Int_t kPSA_neutron;
Int_t kTOF_neutron;

Double_t cutPromptMin,cutPromptMax,cutSideMax,cutSideMin,backgroundSubstractionFactor;
Double_t unten_mass, oben_mass, unten_copl, oben_copl, oben_inv, unten_inv, oben_theta, unten_theta;
Double_t timebackground;
Double_t pionmasse;
Double_t pt;
Double_t pt1;
Int_t runNumber;
TNamed filenumber;
TString* path11;


TH3F *events_all;
TH3F *events_witherror;

HistoManu *coplanarity_withoutcuts;

//dE cuts
HistoManu2 *Check_CBdE_E_nocuts_all;
HistoManu2 *Check_CBdE_E_nocuts_CBandPID;
HistoManu2 *Check_CBdE_E_nocuts_CBandMWPC;
HistoManu2 *Check_TAPSdE_E_nocuts;

//TOF
HistoManu2 *Check_TAPS_TOF_piplus;
HistoManu2 *Check_TAPS_TOF_neutron;

HistoManu2 *Check_TAPS_energy_time_neutron;
HistoManu2 *Check_TAPS_energy_time_piplus;

HistoManu2 *Check_PSA_together;

//PSA
HistoManu2 *Check_PSA_piplus;
HistoManu2 *Check_PSA_neutron;

HistoManu *theta_neutron;
HistoManu *theta_all_neutron;

//Clustersize cuts
HistoManu3 *clustersize_cospiplus_energy_collerated_piplus_CB;
HistoManu3 *clustersize_cospiplus_energy_collerated_piplus_TAPS;
HistoManu3 *clustersize_cospiplus_energy_collerated_neutron_CB;
HistoManu3 *clustersize_cospiplus_energy_collerated_neutron_TAPS;

string polsetting;
string planesetting;
TH1F* test;

HistoManu2 *cosverteilung_collerated;
HistoManu2 *cosverteilung_forE_collerated;

HistoManu3 *kristallminus_targetplus_collerated;
HistoManu3 *kristallminus_targetminus_collerated;
HistoManu3 *kristallplus_targetplus_collerated;
HistoManu3 *kristallplus_targetminus_collerated;

HistoManu3 *kohlenstoffasymmetrie_collerated;
HistoManu3 *kohlenstoffasymmetrie_collerated_neutron;

HistoManu3 *kristallminus_targetplus_collerated_pb;
HistoManu3 *kristallminus_targetminus_collerated_pb;
HistoManu3 *kristallplus_targetplus_collerated_pb;
HistoManu3 *kristallplus_targetminus_collerated_pb;

HistoManu3 *kristallminus_targetplus_collerated_pt;
HistoManu3 *kristallminus_targetminus_collerated_pt;
HistoManu3 *kristallplus_targetplus_collerated_pt;
HistoManu3 *kristallplus_targetminus_collerated_pt;

HistoManu2 *cosverteilung_collerated_neutron;
HistoManu2 *cosverteilung_forE_collerated_neutron;
HistoManu3 *kristallminus_targetplus_collerated_neutron;
HistoManu3 *kristallminus_targetminus_collerated_neutron;
HistoManu3 *kristallplus_targetplus_collerated_neutron;
HistoManu3 *kristallplus_targetminus_collerated_neutron;

HistoManu3 *kristallminus_targetplus_collerated_pb_neutron;
HistoManu3 *kristallminus_targetminus_collerated_pb_neutron;
HistoManu3 *kristallplus_targetplus_collerated_pb_neutron;
HistoManu3 *kristallplus_targetminus_collerated_pb_neutron;

HistoManu3 *kristallminus_targetplus_collerated_pt_neutron;
HistoManu3 *kristallminus_targetminus_collerated_pt_neutron;
HistoManu3 *kristallplus_targetplus_collerated_pt_neutron;
HistoManu3 *kristallplus_targetminus_collerated_pt_neutron;

HistoManu2 *cos_beam_hel1_targetplus;
HistoManu2 *cos_beam_hel0_targetplus;
HistoManu2 *cos_beam_hel1_targetminus;
HistoManu2 *cos_beam_hel0_targetminus;

HistoManu2 *cos_beam_hel1_targetplus_pc;
HistoManu2 *cos_beam_hel0_targetplus_pc;
HistoManu2 *cos_beam_hel1_targetminus_pc;
HistoManu2 *cos_beam_hel0_targetminus_pc;

HistoManu2 *cos_beam_hel1_targetplus_pt;
HistoManu2 *cos_beam_hel0_targetplus_pt;
HistoManu2 *cos_beam_hel1_targetminus_pt;
HistoManu2 *cos_beam_hel0_targetminus_pt;

HistoManu2 *cos_beam_hel1_targetplus_neutron;
HistoManu2 *cos_beam_hel0_targetplus_neutron;
HistoManu2 *cos_beam_hel1_targetminus_neutron;
HistoManu2 *cos_beam_hel0_targetminus_neutron;

HistoManu2 *cos_beam_hel1_targetplus_neutron_pc;
HistoManu2 *cos_beam_hel0_targetplus_neutron_pc;
HistoManu2 *cos_beam_hel1_targetminus_neutron_pc;
HistoManu2 *cos_beam_hel0_targetminus_neutron_pc;

HistoManu2 *cos_beam_hel1_targetplus_neutron_pt;
HistoManu2 *cos_beam_hel0_targetplus_neutron_pt;
HistoManu2 *cos_beam_hel1_targetminus_neutron_pt;
HistoManu2 *cos_beam_hel0_targetminus_neutron_pt;

TH1F *triggertest;
TH1F *poltable_energy;
TH1F *poltable_energy_weight;

HistoManu *h_energy_sum_piplus;
HistoManu *h_energy_sum_piplus_2_3ped;

HistoManu *beamenergy_gen;
HistoManu2 *cosmeson_beamenergy_monte;
HistoManu2 *cosmeson_beamenergy_rek_3ped;
HistoManu2 *cosmeson_beamenergy_rek_2_3ped;

HistoManu3 *cosmeson_beamenergy_energysum_rek;
HistoManu3 *cosmeson_beamenergy_energysum_rek_neutron;
HistoManu3 *cosmeson_beamenergy_energysum_monte;


HistoManu3 *thetaneutron_cospiplus_energy_collerated;
HistoManu3 *thetaneutron_cospiplus_energy_collerated_withoutcuts;

TH1F* h_anzahl_piplus;


// 	Double_t test1;
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual void fOnEndProcessing();

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
  	in2.open("/hadron/spieker/Mainz_Analyse/allinfos_recentbeamtimes.txt");
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
  	in2.open("/hadron/spieker/Mainz_Analyse/allinfos_recentbeamtimes.txt");

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
  	in2.open("/hadron/spieker/Mainz_Analyse/allinfos_recentbeamtimes.txt");

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
    NPiPlusAnalyse();
    virtual ~NPiPlusAnalyse();

    //virtual Bool_t	Init(const char* configfile);

};
#endif

