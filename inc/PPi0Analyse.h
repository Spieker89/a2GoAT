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
#include "HistoManu.h"
#include "HistoManu2.h"
#include "HistoManu3.h"
#include "TCutG.h"
#include <sstream>
#include "HistoManu.h"
#include "HistoManu2.h"
#include "HistoManu3.h"
#include "TTree.h"

class	PPi0Analyse  : public PPhysics
{
protected:
    TH1F*	time1;
    TH1F*	time_prompt;
    TH1F*	time_side;

    HistoManu*        IM_all;
    HistoManu*        MM_all;  
    HistoManu*        IM;
    HistoManu*        MM;

    HistoManu*        theta_all_proton;
    HistoManu*        coplanarity_all_proton;
    HistoManu*        IM_all_proton;
    HistoManu*        MM_all_proton;  

    HistoManu*        IM_proton;
    HistoManu*        MM_proton;
     HistoManu*        theta_proton;
     HistoManu*        coplanarity_proton;

    HistoManu2* coplanarity_theta_CMS;
    HistoManu2* coplanarity_theta_CMS_all;

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
HistoManu3 *IM_energy_kplustplus;
HistoManu3 *IM_energy_kplustminus;
HistoManu3 *IM_energy_kminustplus;
HistoManu3 *IM_energy_kminustminus;

HistoManu3 *MM_energy_kplustplus;
HistoManu3 *MM_energy_kplustminus;
HistoManu3 *MM_energy_kminustplus;
HistoManu3 *MM_energy_kminustminus;

HistoManu2 *invmassverteilung_collerated;
HistoManu2 *invmassverteilung_collerated_proton;

// //kinematic variables in dependence of energy
HistoManu3 *IM_energy_kplustplus_proton;
HistoManu3 *IM_energy_kplustminus_proton;
HistoManu3 *IM_energy_kminustplus_proton;
HistoManu3 *IM_energy_kminustminus_proton;

HistoManu3 *MM_energy_kplustplus_proton;
HistoManu3 *MM_energy_kplustminus_proton;
HistoManu3 *MM_energy_kminustplus_proton;
HistoManu3 *MM_energy_kminustminus_proton;

TH3F *events_all;
TH3F *events_witherror;
HistoManu2 *Check_CBdE_E_nocuts_all;
HistoManu2 *Check_CBdE_E_nocuts_CBandPID;
HistoManu2 *Check_CBdE_E_nocuts_CBandMWPC;
HistoManu2 *Check_TAPSdE_E_nocuts;

HistoManu2 *Check_CBdE_E_all;
HistoManu2 *Check_CBdE_E_CBandPID;
HistoManu2 *Check_CBdE_E_CBandMWPC;
HistoManu2 *Check_TAPSdE_E;

HistoManu2 *Check_TAPS_TOF_proton;
HistoManu2 *Check_TAPS_TOF_photon_3ped;
HistoManu2 *Check_TAPS_TOF_photon_2_3ped;


string polsetting;
string planesetting;
TH1F* test;

ofstream calfile;
ofstream calfile1;
HistoManu2 *thetaproton_detector_collerated;

HistoManu2 *cosverteilung_collerated;
HistoManu2 *cosverteilung_collerated_pt;
HistoManu2 *cosverteilung_collerated_pbminus;
HistoManu2 *cosverteilung_collerated_pbplus;
HistoManu2 *cosverteilung_collerated_minus;
HistoManu2 *cosverteilung_collerated_plus;

HistoManu2 *cosverteilung_forE_collerated;
HistoManu2 *coplanarityverteilung_collerated;
HistoManu2 *coplanarityverteilung_forE_collerated;
HistoManu2 *thetaverteilung_collerated;
HistoManu2 *missingmassverteilung_collerated;
HistoManu2 *missingmassverteilung_forE_collerated;
HistoManu2 *missingmassverteilung_collerated_proton;
HistoManu2 *missingmassverteilung_forE_collerated_proton;
HistoManu3 *kristallminus_targetplus_collerated;
HistoManu3 *kristallminus_targetminus_collerated;
HistoManu3 *kristallplus_targetplus_collerated;
HistoManu3 *kristallplus_targetminus_collerated;

HistoManu3 *kohlenstoffasymmetrie_collerated;
HistoManu3 *kohlenstoffasymmetrie_collerated_proton;

HistoManu3 *kristallminus_targetplus_collerated_pb;
HistoManu3 *kristallminus_targetminus_collerated_pb;
HistoManu3 *kristallplus_targetplus_collerated_pb;
HistoManu3 *kristallplus_targetminus_collerated_pb;

HistoManu3 *kristallminus_targetplus_collerated_pt;
HistoManu3 *kristallminus_targetminus_collerated_pt;
HistoManu3 *kristallplus_targetplus_collerated_pt;
HistoManu3 *kristallplus_targetminus_collerated_pt;

HistoManu2 *cosverteilung_collerated_proton;
HistoManu2 *cosverteilung_collerated_proton_pbminus;
HistoManu2 *cosverteilung_collerated_proton_pbplus;
HistoManu2 *cosverteilung_collerated_proton_minus;
HistoManu2 *cosverteilung_collerated_proton_plus;
HistoManu2 *cosverteilung_collerated_proton_pt;
HistoManu2 *cosverteilung_forE_collerated_proton;
HistoManu3 *kristallminus_targetplus_collerated_proton;
HistoManu3 *kristallminus_targetminus_collerated_proton;
HistoManu3 *kristallplus_targetplus_collerated_proton;
HistoManu3 *kristallplus_targetminus_collerated_proton;

HistoManu3 *kristallminus_targetplus_collerated_pb_proton;
HistoManu3 *kristallminus_targetminus_collerated_pb_proton;
HistoManu3 *kristallplus_targetplus_collerated_pb_proton;
HistoManu3 *kristallplus_targetminus_collerated_pb_proton;

HistoManu3 *kristallminus_targetplus_collerated_pt_proton;
HistoManu3 *kristallminus_targetminus_collerated_pt_proton;
HistoManu3 *kristallplus_targetplus_collerated_pt_proton;
HistoManu3 *kristallplus_targetminus_collerated_pt_proton;

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

HistoManu2 *cos_beam_hel1_targetplus_proton;
HistoManu2 *cos_beam_hel0_targetplus_proton;
HistoManu2 *cos_beam_hel1_targetminus_proton;
HistoManu2 *cos_beam_hel0_targetminus_proton;

HistoManu2 *cos_beam_hel1_targetplus_proton_pc;
HistoManu2 *cos_beam_hel0_targetplus_proton_pc;
HistoManu2 *cos_beam_hel1_targetminus_proton_pc;
HistoManu2 *cos_beam_hel0_targetminus_proton_pc;

HistoManu2 *cos_beam_hel1_targetplus_proton_pt;
HistoManu2 *cos_beam_hel0_targetplus_proton_pt;
HistoManu2 *cos_beam_hel1_targetminus_proton_pt;
HistoManu2 *cos_beam_hel0_targetminus_proton_pt;

TH1F *triggertest;
TH1F *poltable_energy;
TH1F *poltable_energy_weight;

HistoManu2 *openingangle_ptopi0_energy_allcuts;
HistoManu2 *openingangle_gammatogamma_energy_allcuts;
HistoManu2 *openingangle_gammatogamma_energy_allcuts_proton;
HistoManu2 *openingangle_ptopi0_energy_nocuts;
HistoManu2 *openingangle_gammatogamma_energy_nocuts;

HistoManu3 *coplanarityverteilung_cospi0_collerated;
HistoManu3 *thetaverteilung_cospi0_collerated;
HistoManu3 *missingmassverteilung_cospi0_collerated;
HistoManu3 *missingmassverteilung_cospi0_collerated_proton;
HistoManu3 *invmassverteilung_cospi0_collerated;
HistoManu3 *invmassverteilung_cospi0_collerated_proton;
HistoManu *h_energy_sum_pi0_3ped;
HistoManu *h_energy_sum_pi0_2_3ped;
HistoManu *beamenergy_gen;
HistoManu2 *cosmeson_beamenergy_monte;
HistoManu2 *cosmeson_beamenergy_rek_3ped;
HistoManu2 *cosmeson_beamenergy_rek_2_3ped;

HistoManu3 *cosmeson_beamenergy_energysum_rek;
HistoManu3 *cosmeson_beamenergy_energysum_rek_proton;
HistoManu3 *cosmeson_beamenergy_energysum_monte;

HistoManu3 *clustersize_cospi0_energy_collerated;
HistoManu3 *clustersize_cospi0_energy_collerated_proton;
HistoManu3 *thetaproton_cospi0_energy_collerated;
HistoManu3 *thetaproton_cospi0_energy_collerated_withoutcuts;
HistoManu2 *thetaproton_phiproton_collerated;
TH1F* h_anzahl_pi0;
HistoManu2* Check_PSA_3ped;


//Copl
TH2F* cut_copldown;
TH2F* cut_coplup;

//Theta
TH2F* cut_thetadown;
TH2F* cut_thetaup;

//Missing Mass
TH2F* cut_mmdown;
TH2F* cut_mmup;

//inv mass
TH2F* cut_invdown;
TH2F* cut_invup;

// outputFile->cd();
// 
// TTree treeanalysis("treeanalysis","a Tree with selected data");
// 
int isPARA, isPERP;
int isSg, isBg;
double Phi; 
double cosT, Eg; 

// 
// treeanalysis.Branch("isPARA", &isPARA);  
// treeanalysis.Branch("isPERP", &isPERP);
// 
// treeanalysis.Branch("isSg", &isSg);  
// treeanalysis.Branch("isBg", &isBg);
// 
// treeanalysis.Branch("cosT",&cosT); 
// treeanalysis.Branch("Eg",&Eg);
// treeanalysis.Branch("Phi", &Phi);

Int_t nparticles=0;

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
			
double dwq(double energy, double cosvalue, TString meson, TString variable) {
	Double_t wq_value=0.;

///////////////////////////////////////////////////////////DWQ from Bonn Gatchina PWA//////////////////////////////////////////////////////////
	Int_t e_pwa[110]= {200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1040,1060,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1340,1360,1380,1400,1420,1440,1460,1480,1500,1520,1540,1560,1580,1600,1620,1640,1660,1680,1700,1720,1740,1760,1780,1800,1820,1840,1860,1880,2000,2020,2040,2060,2080,2100,2120,2140,2160,2180,2200,2220,2240,2260,2280,2300,2320,2340,2360,2380,2400};

	TString se_pwa[110]= {"e200","e220","e240","e260","e280","e300","e320","e340","e360","e380","e400","e420","e440","e460","e480","e500","e520","e540","e560","e580","e600","e620","e640","e660","e680","e700","e720","e740","e760","e780","e800","e820","e840","e860","e880","e900","e920","e940","e960","e980","e1000","e1020","e1040","e1060","e1080","e1100","e1120","e1140","e1160","e1180","e1200","e1220","e1240","e1260","e1280","e1300","e1320","e1340","e1360","e1380","e1400","e1420","e1440","e1460","e1480","e1500","e1520","e1540","e1560","e1580","e1600","e1620","e1640","e1660","e1680","e1700","e1720","e1740","e1760","e1780","e1800","e1820","e1840","e1860","e1880","e2000","e2020","e2040","e2060","e2080","e2100","e2120","e2140","e2160","e2180","e2200","e2220","e2240","e2260","e2280","e2300","e2320","e2340","e2360","e2380","e2400"};

	Int_t e_pwa_eta[84]= {720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1040,1060,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1340,1360,1380,1400,1420,1440,1460,1480,1500,1520,1540,1560,1580,1600,1620,1640,1660,1680,1700,1720,1740,1760,1780,1800,1820,1840,1860,1880,2000,2020,2040,2060,2080,2100,2120,2140,2160,2180,2200,2220,2240,2260,2280,2300,2320,2340,2360,2380,2400};

	TString se_pwa_eta[84]= {"e720","e740","e760","e780","e800","e820","e840","e860","e880","e900","e920","e940","e960","e980","e1000","e1020","e1040","e1060","e1080","e1100","e1120","e1140","e1160","e1180","e1200","e1220","e1240","e1260","e1280","e1300","e1320","e1340","e1360","e1380","e1400","e1420","e1440","e1460","e1480","e1500","e1520","e1540","e1560","e1580","e1600","e1620","e1640","e1660","e1680","e1700","e1720","e1740","e1760","e1780","e1800","e1820","e1840","e1860","e1880","e2000","e2020","e2040","e2060","e2080","e2100","e2120","e2140","e2160","e2180","e2200","e2220","e2240","e2260","e2280","e2300","e2320","e2340","e2360","e2380","e2400"};

	Int_t e_pwa_2pi0[22]= {375,425,475,525,575,625,675,725,775,825,875,925,975,1025,1075,1125,1175,1225,1275,1325,1375,1425};

	TString se_pwa_2pi0[22]= {"e375","e425","e475","e525","e575","e625","e675","e725","e775","e825","e875","e925","e975","e1025","e1075","e1125","e1175","e1225","e1275","e1325","e1375","e1425"};


	Double_t cosbins[43] = {  -1.00,-0.99,-0.95,-0.90,-0.85,-0.80 ,-0.75,-0.70 ,-0.65 ,-0.60 ,-0.55,-0.50 ,-0.45 ,-0.40 ,-0.35,-0.30,-0.25,-0.20,-0.15,-0.10 ,-0.05 , 0.00 ,0.05 ,0.10,0.15,0.20 ,0.25,0.30 ,0.35, 0.40,0.45,0.50, 0.55, 0.60,0.65 ,0.70 , 0.75, 0.80 , 0.85, 0.90, 0.95 , 0.99, 1.00 };

	Double_t cosbins_2pi0[20] = {-0.95,-0.85,-0.75,-0.65 ,-0.55,-0.45 ,-0.35,-0.25,-0.15,0.05 ,0.15,0.25 ,0.35,0.45, 0.55,0.65, 0.75, 0.85, 0.95};
///////////////////////////////////////////////////////////DWQ from Bonn Gatchina PWA//////////////////////////////////////////////////////////


	if(meson=="pi0"){

		for(Int_t i=0;i<109;i++) {
			if((e_pwa[i]<energy)&&(e_pwa[i+1]>energy)){
				Double_t diffOgrenze = TMath::Abs(e_pwa[i+1]-energy);
				Double_t diffUgrenze = TMath::Abs(e_pwa[i]-energy);

				if(diffOgrenze > diffUgrenze){

					TString Result= se_pwa[i];

					FILE*  pwa;
  					pwa=fopen("/hadron/afzal/root/Aktuell/pi0wq/"+Result,"r");//+index //+index num2str(e_pwa[i])//+index num2str(e_pwa[i])num2str(e_pwa[i])
  					Float_t cos_pwa[2][43];

  					for(Int_t i=0;i<43;i++){//+index num2str(e_pwa[i])
						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);
					}

  					for(Int_t i=0;i<42;i++){//+index num2str(e_pwa[i])
// 						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);
						if((cosbins[i]<cosvalue)&&(cosbins[i+1]>cosvalue)){
							Double_t diffOg_cos = TMath::Abs(cosbins[i+1]-cosvalue);
							Double_t diffUg_cos = TMath::Abs(cosbins[i]-cosvalue);
							if(diffOg_cos > diffUg_cos){
								wq_value= cos_pwa[1][i];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
							if(diffOg_cos < diffUg_cos){
								wq_value= cos_pwa[1][i+1];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;

							}
						}
					}
  					fclose(pwa);
				}

				else{
					TString Result= se_pwa[i+1];

					FILE*  pwa;
  					pwa=fopen("/hadron/afzal/root/Aktuell/pi0wq/"+Result,"r");
  					Float_t cos_pwa[2][43];

  					for(Int_t i=0;i<43;i++){//+index num2str(e_pwa[i])
						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);
					}

  					for(Int_t i=0;i<42;i++){
// 						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);

						if((cosbins[i]<cosvalue)&&(cosbins[i+1]>cosvalue)){
							Double_t diffOg_cos = TMath::Abs(cosbins[i+1]-cosvalue);
							Double_t diffUg_cos = TMath::Abs(cosbins[i]-cosvalue);
							if(diffOg_cos > diffUg_cos){
								wq_value= cos_pwa[1][i];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
							if(diffOg_cos < diffUg_cos){
								wq_value= cos_pwa[1][i+1];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
						}
					}
					fclose(pwa);
					
				}
			}
		}
	}
	if(meson=="eta"){
		for(Int_t i=0;i<83;i++) {
			if((e_pwa_eta[i]<energy)&&(e_pwa_eta[i+1]>energy)){
				Double_t diffOgrenze = TMath::Abs(e_pwa_eta[i+1]-energy);
				Double_t diffUgrenze = TMath::Abs(e_pwa_eta[i]-energy);

				if(diffOgrenze > diffUgrenze){

					TString Result= se_pwa_eta[i];

					FILE*  pwa;
  					pwa=fopen("/hadron/afzal/root/Aktuell/etawq/"+Result,"r");//+index //+index num2str(e_pwa[i])//+index num2str(e_pwa[i])num2str(e_pwa[i])
  					Float_t cos_pwa[2][43];

  					for(Int_t i=0;i<43;i++){//+index num2str(e_pwa[i])
						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);
					}

  					for(Int_t i=0;i<42;i++){//+index num2str(e_pwa[i])
// 						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);
						if((cosbins[i]<cosvalue)&&(cosbins[i+1]>cosvalue)){
							Double_t diffOg_cos = TMath::Abs(cosbins[i+1]-cosvalue);
							Double_t diffUg_cos = TMath::Abs(cosbins[i]-cosvalue);
							if(diffOg_cos > diffUg_cos){
								wq_value= cos_pwa[1][i];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
							if(diffOg_cos < diffUg_cos){
								wq_value= cos_pwa[1][i+1];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;

							}
						}
					}
  					fclose(pwa);
				}

				else{
					TString Result= se_pwa_eta[i+1];

					FILE*  pwa;
  					pwa=fopen("/hadron/afzal/root/Aktuell/etawq/"+Result,"r");
  					Float_t cos_pwa[2][43];

  					for(Int_t i=0;i<43;i++){//+index num2str(e_pwa[i])
						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);
					}

  					for(Int_t i=0;i<42;i++){
// 						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);

						if((cosbins[i]<cosvalue)&&(cosbins[i+1]>cosvalue)){
							Double_t diffOg_cos = TMath::Abs(cosbins[i+1]-cosvalue);
							Double_t diffUg_cos = TMath::Abs(cosbins[i]-cosvalue);
							if(diffOg_cos > diffUg_cos){
								wq_value= cos_pwa[1][i];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
							if(diffOg_cos < diffUg_cos){
								wq_value= cos_pwa[1][i+1];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
						}
					}
					fclose(pwa);
					
				}
			}
		}

	}
	if(meson=="2pi0"){

		for(Int_t i=0;i<21;i++) {
			if((e_pwa_2pi0[i]<energy)&&(e_pwa_2pi0[i+1]>energy)){
				Double_t diffOgrenze = TMath::Abs(e_pwa_2pi0[i+1]-energy);
				Double_t diffUgrenze = TMath::Abs(e_pwa_2pi0[i]-energy);

				if(diffOgrenze > diffUgrenze){

					TString Result= se_pwa_2pi0[i];

					FILE*  pwa;
  					if(variable=="cosppi0"){pwa=fopen("/hadron/spieker/Mainz_Analyse/2pi0wq/"+Result+"_cosppi0","r");//+index //+index num2str(e_pwa[i])//+index num2str(e_pwa[i])num2str(e_pwa[i])
					
					}else{
  					pwa=fopen("/hadron/spieker/Mainz_Analyse/2pi0wq/"+Result,"r");//+index //+index num2str(e_pwa[i])//+index num2str(e_pwa[i])num2str(e_pwa[i])
  					}
					Float_t cos_pwa[3][20];

  					for(Int_t i=0;i<20;i++){//+index num2str(e_pwa[i])
						fscanf(pwa,"%f %f %f\n" ,&cos_pwa[0][i],&cos_pwa[1][i],&cos_pwa[2][i]);
					}

  					for(Int_t i=0;i<19;i++){//+index num2str(e_pwa[i])
// 						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);
						if((cosbins_2pi0[i]<cosvalue)&&(cosbins_2pi0[i+1]>cosvalue)){
							Double_t diffOg_cos = TMath::Abs(cosbins_2pi0[i+1]-cosvalue);
							Double_t diffUg_cos = TMath::Abs(cosbins_2pi0[i]-cosvalue);
							if(diffOg_cos > diffUg_cos){
								wq_value= cos_pwa[2][i];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
							if(diffOg_cos < diffUg_cos){
								wq_value= cos_pwa[2][i+1];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;

							}
						}
					}
  					fclose(pwa);
				}

				else{
					TString Result= se_pwa_2pi0[i+1];

					FILE*  pwa;
  					pwa=fopen("/hadron/afzal/root/Aktuell/2pi0wq/"+Result,"r");
  					Float_t cos_pwa[3][20];

  					for(Int_t i=0;i<20;i++){//+index num2str(e_pwa[i])
						fscanf(pwa,"%f %f %f\n" ,&cos_pwa[0][i],&cos_pwa[1][i],&cos_pwa[2][i]);
					}

  					for(Int_t i=0;i<19;i++){
// 						fscanf(pwa,"%f %f \n" ,&cos_pwa[0][i],&cos_pwa[1][i]);

						if((cosbins_2pi0[i]<cosvalue)&&(cosbins_2pi0[i+1]>cosvalue)){
							Double_t diffOg_cos = TMath::Abs(cosbins_2pi0[i+1]-cosvalue);
							Double_t diffUg_cos = TMath::Abs(cosbins_2pi0[i]-cosvalue);
							if(diffOg_cos > diffUg_cos){
								wq_value= cos_pwa[2][i];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
							if(diffOg_cos < diffUg_cos){
								wq_value= cos_pwa[2][i+1];
// 								if(wq_value>200.) std::cerr<< i<< " "<< wq_value<< std::endl;
							}
						}
					}
					fclose(pwa);
					
				}
			}
		}
	}


// 		std::cerr<< wq_value<<std::endl;
	return wq_value;
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
    PPi0Analyse();
    virtual ~PPi0Analyse();

    //virtual Bool_t	Init(const char* configfile);

};
#endif

