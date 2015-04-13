#ifndef __P2Pi0Analyse_h__
#define __P2Pi0Analyse_h__

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

#include <sstream>

class	P2Pi0Analyse  : public PPhysics
{
protected:
    TH1F*	time1;
    TH1F*	time_prompt;
    TH1F*	time_side;

TH3F *events_all;
TH3F *events_witherror;
TH2F *Check_CBdE_E;
TH2F *Check_TAPSdE_E;
Double_t unten_mass, oben_mass, unten_copl, oben_copl, oben_inv, unten_inv, oben_theta, unten_theta;
Double_t timebackground;
Double_t pionmasse;
Double_t pt1;

string polsetting;
string planesetting;
Double_t pt;
TH1F* test;

TH2F *thetaverteilung_collerated;
TH2F *thetaverteilung_collerated_taps;
TH2F *coplanarityverteilung_collerated;
TH2F *missingmassverteilung_collerated;

TH1F *coplanarity_collerated;
TH1F *coplanarity_mass_theta_collerated;
TH1F *coplanarity_mass_theta_inv_collerated;
TH1F *coplanarity_mass_collerated;
TH1F* coplanarity_theta_collerated;
TH1F *coplanarity_inv_collerated;

TH1F *thetaproton_collerated;
TH1F *thetaproton_mass_copl_collerated;
TH1F *thetaproton_mass_copl_inv_collerated;
TH1F *thetaproton_mass_copl_inv1_collerated;
TH1F *thetaproton_mass_collerated;
TH1F *thetaproton_copl_collerated;
TH1F *thetaproton_inv_collerated;
TH1F *thetaproton_taps_collerated;
TH1F *thetaproton_mass_copl_taps_collerated;
TH1F *thetaproton_mass_copl_inv_taps_collerated;
TH1F *thetaproton_mass_taps_collerated;
TH1F *thetaproton_copl_taps_collerated;
TH1F *thetaproton_inv_taps_collerated;

TH1F *missingmass_collerated;
TH1F *missingmass_theta_copl_collerated;
TH1F *missingmass_theta_copl_inv_collerated;
TH1F *missingmass_theta_collerated;
TH1F *missingmass_copl_collerated;
TH1F *missingmass_inv_collerated;

TH2F *massesumme_mass_theta_copl_inv_beam_collerated;
TH1F *massesumme_collerated;
TH1F *massesumme_mass_theta_copl_collerated;
TH1F *massesumme_mass_theta_copl_inv1_collerated;
TH1F *massesumme_mass_theta_copl_inv_collerated;
TH1F *massesumme_mass_copl_collerated;
TH1F *massesumme_mass_collerated;
TH1F *massesumme_mass_collerated_proton;
TH1F *massesumme_copl_collerated;
TH1F *massesumme_theta_collerated;
TH2F *massegegenmasse_collerated;
TH2F *massegegenmasse_mass_copl_collerated;
TH2F *massegegenmasse_mass_theta_copl_collerated;
TH2F *massegegenmasse_copl_collerated;
TH2F *massegegenmasse_theta_collerated;
TH2F *massegegenmasse_mass_collerated;
TH2F *massegegenmasse_mass_collerated_proton;


TH2F *cosverteilung_collerated;
TH2F *cosppi0verteilung_collerated; 
TH2F *invverteilung_collerated;
TH2F *invppi0verteilung_collerated;
TH2F *cosverteilung_thetaproton_collerated;
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

TH3F *kristallminus_targetplus_ppi0_collerated;
TH3F *kristallminus_targetminus_ppi0_collerated;
TH3F *kristallplus_targetplus_ppi0_collerated;
TH3F *kristallplus_targetminus_ppi0_collerated;

TH3F *kristallminus_targetplus_ppi0_collerated_pb;
TH3F *kristallminus_targetminus_ppi0_collerated_pb;
TH3F *kristallplus_targetplus_ppi0_collerated_pb;
TH3F *kristallplus_targetminus_ppi0_collerated_pb;

TH3F *kristallminus_targetplus_ppi0_collerated_pt;
TH3F *kristallminus_targetminus_ppi0_collerated_pt;
TH3F *kristallplus_targetplus_ppi0_collerated_pt;
TH3F *kristallplus_targetminus_ppi0_collerated_pt;

TH3F *kristallminus_targetplus_im_collerated;
TH3F *kristallminus_targetminus_im_collerated;
TH3F *kristallplus_targetplus_im_collerated;
TH3F *kristallplus_targetminus_im_collerated;

TH3F *kristallminus_targetplus_im_collerated_pb;
TH3F *kristallminus_targetminus_im_collerated_pb;
TH3F *kristallplus_targetplus_im_collerated_pb;
TH3F *kristallplus_targetminus_im_collerated_pb;

TH3F *kristallminus_targetplus_im_collerated_pt;
TH3F *kristallminus_targetminus_im_collerated_pt;
TH3F *kristallplus_targetplus_im_collerated_pt;
TH3F *kristallplus_targetminus_im_collerated_pt;

TH3F *kristallminus_targetplus_im_ppi0_collerated;
TH3F *kristallminus_targetminus_im_ppi0_collerated;
TH3F *kristallplus_targetplus_im_ppi0_collerated;
TH3F *kristallplus_targetminus_im_ppi0_collerated;

TH3F *kristallminus_targetplus_im_ppi0_collerated_pb;
TH3F *kristallminus_targetminus_im_ppi0_collerated_pb;
TH3F *kristallplus_targetplus_im_ppi0_collerated_pb;
TH3F *kristallplus_targetminus_im_ppi0_collerated_pb;

TH3F *kristallminus_targetplus_im_ppi0_collerated_pt;
TH3F *kristallminus_targetminus_im_ppi0_collerated_pt;
TH3F *kristallplus_targetplus_im_ppi0_collerated_pt;
TH3F *kristallplus_targetminus_im_ppi0_collerated_pt;


//PROTON IDENTIFIED
TH2F *cosverteilung_collerated_proton;
TH2F *cosppi0verteilung_collerated_proton; 
TH2F *invverteilung_collerated_proton;
TH2F *invppi0verteilung_collerated_proton;
TH2F *cosverteilung_thetaproton_collerated_proton;
TH2F *missingmassverteilung_collerated_proton;

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

TH3F *kristallminus_targetplus_ppi0_collerated_proton;
TH3F *kristallminus_targetminus_ppi0_collerated_proton;
TH3F *kristallplus_targetplus_ppi0_collerated_proton;
TH3F *kristallplus_targetminus_ppi0_collerated_proton;

TH3F *kristallminus_targetplus_ppi0_collerated_pb_proton;
TH3F *kristallminus_targetminus_ppi0_collerated_pb_proton;
TH3F *kristallplus_targetplus_ppi0_collerated_pb_proton;
TH3F *kristallplus_targetminus_ppi0_collerated_pb_proton;

TH3F *kristallminus_targetplus_ppi0_collerated_pt_proton;
TH3F *kristallminus_targetminus_ppi0_collerated_pt_proton;
TH3F *kristallplus_targetplus_ppi0_collerated_pt_proton;
TH3F *kristallplus_targetminus_ppi0_collerated_pt_proton;

TH3F *kristallminus_targetplus_im_collerated_proton;
TH3F *kristallminus_targetminus_im_collerated_proton;
TH3F *kristallplus_targetplus_im_collerated_proton;
TH3F *kristallplus_targetminus_im_collerated_proton;

TH3F *kristallminus_targetplus_im_collerated_pb_proton;
TH3F *kristallminus_targetminus_im_collerated_pb_proton;
TH3F *kristallplus_targetplus_im_collerated_pb_proton;
TH3F *kristallplus_targetminus_im_collerated_pb_proton;

TH3F *kristallminus_targetplus_im_collerated_pt_proton;
TH3F *kristallminus_targetminus_im_collerated_pt_proton;
TH3F *kristallplus_targetplus_im_collerated_pt_proton;
TH3F *kristallplus_targetminus_im_collerated_pt_proton;

TH3F *kristallminus_targetplus_im_ppi0_collerated_proton;
TH3F *kristallminus_targetminus_im_ppi0_collerated_proton;
TH3F *kristallplus_targetplus_im_ppi0_collerated_proton;
TH3F *kristallplus_targetminus_im_ppi0_collerated_proton;

TH3F *kristallminus_targetplus_im_ppi0_collerated_pb_proton;
TH3F *kristallminus_targetminus_im_ppi0_collerated_pb_proton;
TH3F *kristallplus_targetplus_im_ppi0_collerated_pb_proton;
TH3F *kristallplus_targetminus_im_ppi0_collerated_pb_proton;

TH3F *kristallminus_targetplus_im_ppi0_collerated_pt_proton;
TH3F *kristallminus_targetminus_im_ppi0_collerated_pt_proton;
TH3F *kristallplus_targetplus_im_ppi0_collerated_pt_proton;
TH3F *kristallplus_targetminus_im_ppi0_collerated_pt_proton;

TH1F *poltable_energy;
TH1F *poltable_energy_weight;
// 	Double_t test1;
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual void fOnEndProcessing();
    virtual void fOnBeforeEventProcessing();
    virtual Bool_t    Write();
//     virtual Double_t targetpol(TFile *f);

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
// 	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Diamond_May14_linpol.txt");
// 	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");


	Double_t runnumber;
	Double_t pol;
	Double_t runnumber_array[1000];
	Double_t pol_array[1000];
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
// 	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Diamond_May14_linpol.txt");



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
// 	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Diamond_May14_linpol.txt");


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


Double_t invariantemasse(TLorentzVector teilchen1, TLorentzVector teilchen2)
{
return (teilchen1+teilchen2).M();
}

Double_t invariantemasseFehler(TLorentzVector teilchen1, TLorentzVector teilchen2)
{
Double_t E1 = teilchen1.E();
Double_t E2 = teilchen2.E();
Double_t E1_error = 0.03*E1;
Double_t E2_error = 0.03*E2;
Double_t phi1 = teilchen1.Phi();
Double_t phi2 = teilchen2.Phi();
Double_t phi1_error = TMath::DegToRad()*2;
Double_t phi2_error = TMath::DegToRad()*2;

Double_t theta1 = teilchen1.Theta();
Double_t theta2 = teilchen2.Theta();
Double_t theta1_error = TMath::DegToRad()*2;
Double_t theta2_error = TMath::DegToRad()*2;

Double_t cos12 = TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi1-phi2) + TMath::Cos(theta1)*TMath::Cos(theta2);


Double_t Dcos12 = TMath::Sqrt(TMath::Power((TMath::Cos(theta1)*TMath::Sin(theta2)*TMath::Cos(phi1-phi2) - TMath::Sin(theta1)*TMath::Cos(theta2))*theta1_error,2) + 		     			  TMath::Power((TMath::Cos(theta2)*TMath::Sin(theta1)*TMath::Cos(phi1-phi2) - TMath::Sin(theta2)*TMath::Cos(theta1))*theta2_error,2) + 						          TMath::Power((TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(phi1-phi2))*phi1_error,2) + 
	 	  TMath::Power((TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(phi1-phi2))*phi2_error,2));

Double_t m12 = TMath::Sqrt(2*E1*E2*(1-cos12));

Double_t Dinv = TMath::Sqrt(TMath::Power(E2*(1-cos12)*E1_error,2)+TMath::Power(E1*(1-cos12)*E2_error,2) + TMath::Power(E2*E1*Dcos12,2))/m12;

return Dinv;
}

Double_t ChiPionPion(Double_t inv11, TLorentzVector teilchen11, TLorentzVector teilchen22, Double_t inv22, TLorentzVector teilchen33, TLorentzVector teilchen44)
{
return TMath::Power((inv11-134.9766)/invariantemasseFehler(teilchen11,teilchen22),2) +  TMath::Power((inv22-134.9766)/invariantemasseFehler(teilchen33,teilchen44),2);
}

Double_t ChiPionEta(Double_t inv111, TLorentzVector teilchen111, TLorentzVector teilchen222, Double_t inv222, TLorentzVector teilchen333, TLorentzVector teilchen444)
{
return TMath::Power((inv111-134.9766)/invariantemasseFehler(teilchen111,teilchen222),2) +  TMath::Power((inv222-547.853)/invariantemasseFehler(teilchen333,teilchen444),2);
}


			
public:
    P2Pi0Analyse();
    virtual ~P2Pi0Analyse();

    //virtual Bool_t	Init(const char* configfile);

};
#endif

