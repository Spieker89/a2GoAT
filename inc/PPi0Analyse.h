#ifndef __PPi0Analyse_h__
#define __PPi0Analyse_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GTreeLinPol.h"
#include "PPhysics.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

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

//kinematic variables in dependence of energy
TH3F *IM_energy_kplustplus;
TH3F *IM_energy_kplustminus;
TH3F *IM_energy_kminustplus;
TH3F *IM_energy_kminustminus;

TH3F *MM_energy_kplustplus;
TH3F *MM_energy_kplustminus;
TH3F *MM_energy_kminustplus;
TH3F *MM_energy_kminustminus;

// //kinematic variables in dependence of energy
TH3F *IM_energy_kplustplus_proton;
TH3F *IM_energy_kplustminus_proton;
TH3F *IM_energy_kminustplus_proton;
TH3F *IM_energy_kminustminus_proton;

TH3F *MM_energy_kplustplus_proton;
TH3F *MM_energy_kplustminus_proton;
TH3F *MM_energy_kminustplus_proton;
TH3F *MM_energy_kminustminus_proton;

Double_t unten_mass, oben_mass, unten_copl, oben_copl, oben_inv, unten_inv, oben_theta, unten_theta;
Double_t timebackground;
Double_t pionmasse;
Double_t pt;
string polsetting;
string planesetting;
TH1F* test;

TH2F *cosverteilung_collerated;

TH2F *coplanarityverteilung_collerated;
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

TH1F *poltable_energy;
TH1F *poltable_energy_weight;
// 	Double_t test1;
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual void fOnEndProcessing();
    virtual void fOnBeforeEventProcessing();
//     virtual Double_t targetpol(TFile *f);

Double_t targetpol(TFile *f){

	TString* filename = new TString(f->GetPath());
	TString* path1 = new TString(filename->Tokenize("_")->At(filename->Tokenize("_")->GetEntries()-1)->GetName());
	path1->Resize(path1->Length()-7);
	stringstream ss(path1->Data());
	Int_t runnumber1;
	ss >> runnumber1;
	ifstream in2;
	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_2013_11.txt");
	Double_t runnumber;
	Double_t pol;
	Double_t runnumber_array[1000];
	Double_t pol_array[1000];
	Double_t pol1;
	Int_t N = 0;

	if (in2.is_open()) {
		while (1) {
		
		in2 >> runnumber >> pol;
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
	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
	Double_t runnumber;
	string poledge;
	string poledge1;
	Int_t N = 0;
	string plane;
	if (in2.is_open()) {
		while (1) {
		
		in2 >> runnumber >> plane >> poledge;
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
	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
	Double_t runnumber;
	string poledge;
	string poledge1;
	Int_t N = 0;
	string plane;
	string plane1;

	if (in2.is_open()) {
		while (1) {
		
		in2 >> runnumber >> plane >> poledge;
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

