#ifndef __P2Pi0Analyse_h__
#define __P2Pi0Analyse_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 
#include <sstream>
#include "GTreeManager.h"
#include "GTreeLinPol.h"
#include "GTreeTrigger.h"
#include "PPhysics.h"
#include "HistoManu.h"
#include "HistoManu2.h"
#include "HistoManu3.h"
#include "GFit.h"
#include "HistoManu.h"
#include "GTreeParticle.h"

class	P2Pi0Analyse  : public PPhysics
{
private:

GFit	fitter;
GFit	fitter1;
GFit	fitter2;

    TH1F*	time1;
    TH1F*	time_prompt;
    TH1F*	time_side;

TH3F *events_all;
TH3F *events_witherror;
HistoManu2 *Check_CBdE_E;
HistoManu2 *Check_TAPSdE_E;
Double_t unten_mass, oben_mass, unten_copl, oben_copl, oben_inv, unten_inv, oben_theta, unten_theta;
Double_t timebackground;
Double_t pionmasse;
Double_t pt1;
Int_t runNumber;
ofstream calfile;
ofstream calfile1;
TString* path11;

string polsetting;
string planesetting;
Double_t pt;
TH1F* test;


    TFile *photonfile = new TFile("/disk/nobackup/001/spieker/Mainz1/a2GoAT/configfiles/uncertainties/photon_uncertainties.root","OPEN");
    TFile *protonfile = new TFile("/disk/nobackup/001/spieker/Mainz1/a2GoAT/configfiles/uncertainties/proton_uncertainties.root","OPEN");
	
    TH2F *photon_CB_energy = (TH2F*)photonfile->Get("CB_E");
    TH2F *photon_CB_theta = (TH2F*)photonfile->Get("CB_theta");
    TH2F *photon_CB_phi = (TH2F*)photonfile->Get("CB_phi");
    TH2F *photon_TAPS_energy = (TH2F*)photonfile->Get("TAPS_E"); 
    TH2F *photon_TAPS_theta = (TH2F*)photonfile->Get("TAPS_theta");
    TH2F *photon_TAPS_phi = (TH2F*)photonfile->Get("TAPS_phi");
	
    TH2F *proton_TAPS_energy = (TH2F*)protonfile->Get("TAPS_E");
    TH2F *proton_TAPS_theta = (TH2F*)protonfile->Get("TAPS_theta");
    TH2F *proton_TAPS_phi = (TH2F*)protonfile->Get("TAPS_phi");


HistoManu3 *openangle_pi01_to_pi02_energy;
HistoManu2 *thetaverteilung_collerated;
HistoManu2 *thetaverteilung_collerated_taps;
HistoManu2 *coplanarityverteilung_collerated;
HistoManu2 *missingmassverteilung_collerated;

HistoManu *coplanarity_collerated;
HistoManu *coplanarity_mass_theta_collerated;
HistoManu *coplanarity_mass_theta_inv_collerated;
HistoManu *coplanarity_mass_collerated;
HistoManu* coplanarity_theta_collerated;
HistoManu *coplanarity_inv_collerated;
HistoManu *coplanarity_fit_vs_rek;

HistoManu *thetaproton_collerated;
HistoManu *thetaproton_mass_copl_collerated;
HistoManu *thetaproton_mass_copl_inv_collerated;
HistoManu *thetaproton_mass_copl_inv1_collerated;
HistoManu *thetaproton_mass_collerated;
HistoManu *thetaproton_copl_collerated;
HistoManu *thetaproton_inv_collerated;
HistoManu *thetaproton_taps_collerated;
HistoManu *thetaproton_mass_copl_taps_collerated;
HistoManu *thetaproton_mass_copl_inv_taps_collerated;
HistoManu *thetaproton_mass_taps_collerated;
HistoManu *thetaproton_copl_taps_collerated;
HistoManu *thetaproton_inv_taps_collerated;
HistoManu *thetaproton_fit_vs_rek_CB;
HistoManu *thetaproton_fit_vs_rek_TAPS;

HistoManu *missingmass_collerated;
HistoManu *missingmass_collerated_proton;
HistoManu *missingmass_theta_copl_collerated;
HistoManu *missingmass_theta_copl_inv_collerated;
HistoManu *missingmass_theta_collerated;
HistoManu *missingmass_copl_collerated;
HistoManu *missingmass_inv_collerated;
HistoManu *missingmass_inv_collerated_proton;


HistoManu2 *massesumme_mass_theta_copl_inv_beam_collerated;
HistoManu *massesumme_collerated;
HistoManu *massesumme_collerated_proton;
HistoManu *massesumme_mass_theta_copl_collerated;
HistoManu *massesumme_mass_theta_copl_inv1_collerated;
HistoManu *massesumme_mass_theta_copl_inv_collerated;
HistoManu *massesumme_mass_copl_collerated;
HistoManu *massesumme_mass_collerated;
HistoManu *massesumme_mass_collerated_proton;
HistoManu *massesumme_copl_collerated;
HistoManu *massesumme_theta_collerated;
HistoManu2 *massegegenmasse_collerated;
HistoManu2 *massegegenmasse_collerated_proton;
HistoManu2 *massegegenmasse_mass_copl_collerated;
HistoManu2 *massegegenmasse_mass_theta_copl_collerated;
HistoManu2 *massegegenmasse_copl_collerated;
HistoManu2 *massegegenmasse_theta_collerated;
HistoManu2 *massegegenmasse_mass_collerated;
HistoManu2 *massegegenmasse_mass_collerated_proton;
HistoManu2 *massegegenmasse_allcutsandinv_collerated;
HistoManu2 *massegegenmasse_allcuts_collerated;


HistoManu2 *cosverteilung_collerated;
HistoManu2 *cosppi0verteilung_collerated; 
HistoManu2 *invverteilung_collerated;
HistoManu2 *invppi0verteilung_collerated;
HistoManu2 *cosverteilung_thetaproton_collerated;
HistoManu3 *kristallminus_targetplus_collerated;
HistoManu3 *kristallminus_targetminus_collerated;
HistoManu3 *kristallplus_targetplus_collerated;
HistoManu3 *kristallplus_targetminus_collerated;

HistoManu3 *kristallminus_targetplus_collerated_pb;
HistoManu3 *kristallminus_targetminus_collerated_pb;
HistoManu3 *kristallplus_targetplus_collerated_pb;
HistoManu3 *kristallplus_targetminus_collerated_pb;

HistoManu3 *kristallminus_targetplus_collerated_pt;
HistoManu3 *kristallminus_targetminus_collerated_pt;
HistoManu3 *kristallplus_targetplus_collerated_pt;
HistoManu3 *kristallplus_targetminus_collerated_pt;

HistoManu3 *kristallminus_targetplus_ppi0_collerated;
HistoManu3 *kristallminus_targetminus_ppi0_collerated;
HistoManu3 *kristallplus_targetplus_ppi0_collerated;
HistoManu3 *kristallplus_targetminus_ppi0_collerated;

HistoManu3 *kristallminus_targetplus_ppi0_collerated_pb;
HistoManu3 *kristallminus_targetminus_ppi0_collerated_pb;
HistoManu3 *kristallplus_targetplus_ppi0_collerated_pb;
HistoManu3 *kristallplus_targetminus_ppi0_collerated_pb;

HistoManu3 *kristallminus_targetplus_ppi0_collerated_pt;
HistoManu3 *kristallminus_targetminus_ppi0_collerated_pt;
HistoManu3 *kristallplus_targetplus_ppi0_collerated_pt;
HistoManu3 *kristallplus_targetminus_ppi0_collerated_pt;

HistoManu3 *kristallminus_targetplus_im_collerated;
HistoManu3 *kristallminus_targetminus_im_collerated;
HistoManu3 *kristallplus_targetplus_im_collerated;
HistoManu3 *kristallplus_targetminus_im_collerated;

HistoManu3 *kristallminus_targetplus_im_collerated_pb;
HistoManu3 *kristallminus_targetminus_im_collerated_pb;
HistoManu3 *kristallplus_targetplus_im_collerated_pb;
HistoManu3 *kristallplus_targetminus_im_collerated_pb;

HistoManu3 *kristallminus_targetplus_im_collerated_pt;
HistoManu3 *kristallminus_targetminus_im_collerated_pt;
HistoManu3 *kristallplus_targetplus_im_collerated_pt;
HistoManu3 *kristallplus_targetminus_im_collerated_pt;

HistoManu3 *kristallminus_targetplus_im_ppi0_collerated;
HistoManu3 *kristallminus_targetminus_im_ppi0_collerated;
HistoManu3 *kristallplus_targetplus_im_ppi0_collerated;
HistoManu3 *kristallplus_targetminus_im_ppi0_collerated;

HistoManu3 *kristallminus_targetplus_im_ppi0_collerated_pb;
HistoManu3 *kristallminus_targetminus_im_ppi0_collerated_pb;
HistoManu3 *kristallplus_targetplus_im_ppi0_collerated_pb;
HistoManu3 *kristallplus_targetminus_im_ppi0_collerated_pb;

HistoManu3 *kristallminus_targetplus_im_ppi0_collerated_pt;
HistoManu3 *kristallminus_targetminus_im_ppi0_collerated_pt;
HistoManu3 *kristallplus_targetplus_im_ppi0_collerated_pt;
HistoManu3 *kristallplus_targetminus_im_ppi0_collerated_pt;


//PROTON IDENTIFIED
HistoManu2 *cosverteilung_collerated_proton;
HistoManu2 *cosppi0verteilung_collerated_proton; 
HistoManu2 *invverteilung_collerated_proton;
HistoManu2 *invppi0verteilung_collerated_proton;
HistoManu2 *cosverteilung_thetaproton_collerated_proton;
HistoManu2 *missingmassverteilung_collerated_proton;

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

HistoManu3 *kristallminus_targetplus_ppi0_collerated_proton;
HistoManu3 *kristallminus_targetminus_ppi0_collerated_proton;
HistoManu3 *kristallplus_targetplus_ppi0_collerated_proton;
HistoManu3 *kristallplus_targetminus_ppi0_collerated_proton;

HistoManu3 *kristallminus_targetplus_ppi0_collerated_pb_proton;
HistoManu3 *kristallminus_targetminus_ppi0_collerated_pb_proton;
HistoManu3 *kristallplus_targetplus_ppi0_collerated_pb_proton;
HistoManu3 *kristallplus_targetminus_ppi0_collerated_pb_proton;

HistoManu3 *kristallminus_targetplus_ppi0_collerated_pt_proton;
HistoManu3 *kristallminus_targetminus_ppi0_collerated_pt_proton;
HistoManu3 *kristallplus_targetplus_ppi0_collerated_pt_proton;
HistoManu3 *kristallplus_targetminus_ppi0_collerated_pt_proton;

HistoManu3 *kristallminus_targetplus_im_collerated_proton;
HistoManu3 *kristallminus_targetminus_im_collerated_proton;
HistoManu3 *kristallplus_targetplus_im_collerated_proton;
HistoManu3 *kristallplus_targetminus_im_collerated_proton;

HistoManu3 *kristallminus_targetplus_im_collerated_pb_proton;
HistoManu3 *kristallminus_targetminus_im_collerated_pb_proton;
HistoManu3 *kristallplus_targetplus_im_collerated_pb_proton;
HistoManu3 *kristallplus_targetminus_im_collerated_pb_proton;

HistoManu3 *kristallminus_targetplus_im_collerated_pt_proton;
HistoManu3 *kristallminus_targetminus_im_collerated_pt_proton;
HistoManu3 *kristallplus_targetplus_im_collerated_pt_proton;
HistoManu3 *kristallplus_targetminus_im_collerated_pt_proton;

HistoManu3 *kristallminus_targetplus_im_ppi0_collerated_proton;
HistoManu3 *kristallminus_targetminus_im_ppi0_collerated_proton;
HistoManu3 *kristallplus_targetplus_im_ppi0_collerated_proton;
HistoManu3 *kristallplus_targetminus_im_ppi0_collerated_proton;

HistoManu3 *kristallminus_targetplus_im_ppi0_collerated_pb_proton;
HistoManu3 *kristallminus_targetminus_im_ppi0_collerated_pb_proton;
HistoManu3 *kristallplus_targetplus_im_ppi0_collerated_pb_proton;
HistoManu3 *kristallplus_targetminus_im_ppi0_collerated_pb_proton;

HistoManu3 *kristallminus_targetplus_im_ppi0_collerated_pt_proton;
HistoManu3 *kristallminus_targetminus_im_ppi0_collerated_pt_proton;
HistoManu3 *kristallplus_targetplus_im_ppi0_collerated_pt_proton;
HistoManu3 *kristallplus_targetminus_im_ppi0_collerated_pt_proton;

TH1F *poltable_energy;
TH1F *poltable_energy_weight;

HistoManu *h_energy_sum_2pi0_5ped;
HistoManu2 *h_energy_taps_2pi0_5ped;
HistoManu *h_energy_sum_2pi0_5ped_weightcos2pi0;
HistoManu *h_energy_sum_2pi0_5ped_weightcosppi0;


//for kinematic fit and so
TH1F *triggertest;
HistoManu *CL;
HistoManu2 *pulls_4g_CB;
HistoManu2 *pulls_4g_TAPS;
HistoManu2 *pulls_beam;
HistoManu2 *pulls_proton_CB;
HistoManu2 *pulls_proton_TAPS;
HistoManu *invmass_diff_rec_kin;

//generated
TH2F *cos2pi0_beamphoton_monte;
TH3F *cosppi0_beamphoton_energysum_monte;

TH3F *m2pi0_beamphoton_energysum_monte;
TH3F *mppi0_beamphoton_energysum_monte;

//reconstructed
HistoManu2 *cos2pi0_beamphoton_rek;;
HistoManu3 *cosppi0_beamphoton_energysum_rek;

HistoManu3 *m2pi0_beamphoton_energysum_rek;
HistoManu3 *mppi0_beamphoton_energysum_rek;

TH1F *h_bph_e_gen;
// 	Double_t test1;
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual void fOnEndProcessing();
    virtual void fOnBeforeEventProcessing();
    virtual Bool_t    Write();
//     virtual Double_t targetpol(TFile *f);

			
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
// 	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Carbon_2013_linpol_neu.txt");
	in2.open("/hadron/spieker/Mainz_Analyse/Butanol_Diamond_May14_linpol.txt");
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

Double_t invariantemasse(const GTreeParticle *photon, const Int_t particle1,const Int_t particle2)
{
	return (photon->Particle(particle1)+photon->Particle(particle2)).M();
}



Double_t invariantemasseFehler(const GTreeParticle *photon, const Int_t particle1,const Int_t particle2)
{
	Double_t E1 = photon->Particle(particle1).E();
	Double_t E2 = photon->Particle(particle2).E();
	
	Double_t phi1 = photon->Particle(particle1).Phi();
	Double_t phi2 = photon->Particle(particle2).Phi();
	
	Double_t theta1 = photon->Particle(particle1).Theta();
	Double_t theta2 = photon->Particle(particle2).Theta();

	Double_t E1_error,E2_error,phi1_error,phi2_error,theta1_error, theta2_error;

	if(photon->HasCB(particle1)){

		E1_error = E1*1.5*photon_CB_energy->GetBinContent(photon_CB_energy->FindBin(E1,theta1*TMath::RadToDeg()));
		phi1_error = 0.83*TMath::DegToRad()*photon_CB_phi->GetBinContent(photon_CB_phi->FindBin(E1,theta1*TMath::RadToDeg()));
		theta1_error = 0.68*TMath::DegToRad()*photon_CB_theta->GetBinContent(photon_CB_theta->FindBin(E1,theta1*TMath::RadToDeg()));

		//if histo is not finding something
		if(E1_error==0) E1_error=(0.02*(E1/1.0e3)*pow((E1/1.0e3),-0.36))*1.0e3*2.0;
		if(theta1_error==0) theta1_error = 5.2*TMath::DegToRad();
 		if(phi1_error==0) phi1_error = 0.7*theta1_error/TMath::Sin(theta1);

	}

	if(photon->HasTAPS(particle1)){

		E1_error = E1*1*photon_TAPS_energy->GetBinContent(photon_TAPS_energy->FindBin(E1,theta1*TMath::RadToDeg()));
		phi1_error = 1*TMath::DegToRad();
		theta1_error = 2.5*TMath::DegToRad();

		//if histo is not finding something
		if(E1_error==0) E1_error = ((0.018 + 0.008*TMath::Sqrt(E1/1.0e3))*(E1/1.0e3))*1.0e3*2.0;
		if(theta1_error==0) theta1_error = 2*TMath::DegToRad();
 		if(phi1_error==0) phi1_error = 2*TMath::DegToRad();

	}

	if(photon->HasCB(particle2)){

		E2_error = E2*1.5*photon_CB_energy->GetBinContent(photon_CB_energy->FindBin(E2,theta2*TMath::RadToDeg()));
		phi2_error = 0.83*TMath::DegToRad()*photon_CB_phi->GetBinContent(photon_CB_phi->FindBin(E2,theta2*TMath::RadToDeg()));
		theta2_error = 0.68*TMath::DegToRad()*photon_CB_theta->GetBinContent(photon_CB_theta->FindBin(E2,theta2*TMath::RadToDeg()));

		//if histo is not finding something
		if(E2_error==0) E2_error=(0.02*(E2/1.0e3)*pow((E2/1.0e3),-0.36))*1.0e3*2.0;
		if(theta2_error==0) theta2_error = 5.2*TMath::DegToRad();
 		if(phi2_error==0) phi2_error = 0.7*theta2_error/TMath::Sin(theta2);

	}

	if(photon->HasTAPS(particle2)){

		E2_error = E2*1*photon_TAPS_energy->GetBinContent(photon_TAPS_energy->FindBin(E2,theta2*TMath::RadToDeg()));
		phi2_error = 1*TMath::DegToRad();
		theta2_error = 2.5*TMath::DegToRad();

		//if histo is not finding something
		if(E2_error==0) E2_error = ((0.018 + 0.008*TMath::Sqrt(E2/1.0e3))*(E2/1.0e3))*1.0e3*2.0;
		if(theta2_error==0) theta2_error = 2*TMath::DegToRad();
 		if(phi2_error==0) phi2_error = 2*TMath::DegToRad();

	}


	Double_t cos12 = TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi1-phi2) + TMath::Cos(theta1)*TMath::Cos(theta2);
	
	
	Double_t Dcos12 = TMath::Sqrt(TMath::Power((TMath::Cos(theta1)*TMath::Sin(theta2)*TMath::Cos(phi1-phi2) - TMath::Sin(theta1)*TMath::Cos(theta2))*theta1_error,2) + 		     			  TMath::Power((TMath::Cos(theta2)*TMath::Sin(theta1)*TMath::Cos(phi1-phi2) - TMath::Sin(theta2)*TMath::Cos(theta1))*theta2_error,2) + 						          TMath::Power((TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(phi1-phi2))*phi1_error,2) + 
			TMath::Power((TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(phi1-phi2))*phi2_error,2));
	
	Double_t m12 = TMath::Sqrt(2*E1*E2*(1-cos12));
	
	Double_t Dinv = TMath::Sqrt(TMath::Power(E2*(1-cos12)*E1_error,2)+TMath::Power(E1*(1-cos12)*E2_error,2) + TMath::Power(E2*E1*Dcos12,2))/m12;

return Dinv;
}

Double_t ChiPionPion(const GTreeParticle *photon,Double_t inv11, Int_t teilchen11, Int_t teilchen22, Double_t inv22, Int_t teilchen33, Int_t teilchen44)
{
return TMath::Power((inv11-134.9766)/invariantemasseFehler(photon,teilchen11,teilchen22),2) +  TMath::Power((inv22-134.9766)/invariantemasseFehler(photon,teilchen33,teilchen44),2);
}

// Double_t ChiPionEta(Double_t inv111, TLorentzVector teilchen111, TLorentzVector teilchen222, Double_t inv222, TLorentzVector teilchen333, TLorentzVector teilchen444)
// {
// return TMath::Power((inv111-134.9766)/invariantemasseFehler(teilchen111,teilchen222),2) +  TMath::Power((inv222-547.853)/invariantemasseFehler(teilchen333,teilchen444),2);
// }


			
public:
    P2Pi0Analyse();
    virtual ~P2Pi0Analyse();

    //virtual Bool_t	Init(const char* configfile);

};
#endif
