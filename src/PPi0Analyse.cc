#include "PPi0Analyse.h"
#include <iostream> 

PPi0Analyse::PPi0Analyse()
{ 

time_prompt = new TH1F("time_prompt", 	"time_prompt", 	1400, -700, 700);
time1= new TH1F("time1", 	"time1", 	1400, -700, 700);
time_side 	= new TH1F("time_side", 	"time_side", 	1400, -700, 700);

 
IM 		= new HistoManu("IM", 	"IM", 		1000,   0, 800);
IM_all          = new HistoManu("IM_all",         "IM_all",           1000,   0, 1000);
IM_proton 		= new HistoManu("IM_proton", 	"IM_proton", 		1000,   0, 1000);
IM_all_proton          = new HistoManu("IM_all_proton",         "IM_all_proton",           1000,   0, 1000);

MM	= new HistoManu("MM", 	"MM", 	 	400,   800, 1200);
MM_all          = new HistoManu("MM_all",         "MM_all",           400,   800, 1200);	
MM_proton		= new HistoManu("MM_proton", 	"MM_proton", 	 	400,   800, 1200);
MM_all_proton          = new HistoManu("MM_all_proton",         "MM_all_proton",           400,   800, 1200);

theta_proton 	= new HistoManu("theta_proton", 	"theta_proton", 		400,   -200, 200);
theta_all_proton       = new HistoManu("theta_all_proton",      "theta_all_proton",                400,   -200, 200);

coplanarity_proton	= new HistoManu("coplanarity_proton","coplanarity_proton",400,0,360);
coplanarity_all_proton = new HistoManu("coplanarity_all_proton","coplanarity_all_proton",400,0,360);

coplanarity_theta_CMS	= new HistoManu2("coplanarity_theta_CMS","coplanarity_theta_CMS",400,0,360,30,220,1420);
coplanarity_theta_CMS_all = new HistoManu2("coplanarity_theta_CMS_all","coplanarity_theta_CMS_all",400,0,360,30,220,1420);
 
poltable_energy          = new HistoManu("poltable_energy",         "poltable_energy",           352,   0, 1448);
poltable_energy_weight          = new HistoManu("poltable_energy_weight",         "poltable_energy_weight",           352,   0, 1448);

h_anzahl_pi0  = new TH1F("h_anzahl_pi0 ", "h_anzahl_pi0 ",  100, -4, 10);
h_anzahl_pi0->Sumw2();

Check_PSA_3ped= new HistoManu2("Check_PSA_3ped ", "Check_PSA_3ped ",  180, 0, 180, 100, 0, 600);




//dE over E spectra
Check_CBdE_E_nocuts_all		= new HistoManu2("Check_CBdE_E_nocuts_all", "dE_E (all detectors compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID}*sin(#theta) [MeV]", 	400, 0, 400, 100, 0, 10);
Check_CBdE_E_nocuts_CBandPID		= new HistoManu2("Check_CBdE_E_nocuts_CBandPID", "dE_E (CB+PID compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID} [MeV]", 	400, 0, 400, 100, 0, 10);
Check_CBdE_E_nocuts_CBandMWPC		= new HistoManu2("Check_CBdE_E_nocuts_CBandMWPC", "dE_E (CB+MWPC compared to MWPC hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{MWPC}*sin(#theta) [MeV]", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E_nocuts		= new HistoManu2("Check_TAPSdE_E_nocuts", "dE_E (all TAPS clusters compared to Veto hits);E_{dep}^{TAPS} [MeV];#Delta E_{dep}^{VETO} [MeV]", 	400, 0, 400, 100, 0, 10);

//dE over E spectra
Check_CBdE_E_all		= new HistoManu2("Check_CBdE_E_all", "dE_E (all detectors compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID}*sin(#theta) [MeV]", 	400, 0, 400, 100, 0, 10);
Check_CBdE_E_CBandPID		= new HistoManu2("Check_CBdE_E_CBandPID", "dE_E (CB+PID compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID} [MeV]", 	400, 0, 400, 100, 0, 10);
Check_CBdE_E_CBandMWPC		= new HistoManu2("Check_CBdE_E_CBandMWPC", "dE_E (CB+MWPC compared to MWPC hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{MWPC}*sin(#theta) [MeV]", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E		= new HistoManu2("Check_TAPSdE_E", "dE_E (all TAPS clusters compared to Veto hits);E_{dep}^{TAPS} [MeV];#Delta E_{dep}^{VETO} [MeV]", 	400, 0, 400, 100, 0, 10);


Check_TAPS_TOF_proton		= new HistoManu2("Check_TAPS_TOF_proton", "TOF analysis (proton);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);
Check_TAPS_TOF_photon_3ped		= new HistoManu2("Check_TAPS_TOF_photon_3ped", "TOF analysis (photon_3ped);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);
Check_TAPS_TOF_photon_2_3ped		= new HistoManu2("Check_TAPS_TOF_photon_2_3ped", "TOF analysis (photon_2_3ped);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);

openingangle_ptopi0_energy_allcuts = new HistoManu2("openingangle_ptopi0_energy_allcuts","openingangle_ptopi0_energy_allcuts;#angle #pi,p;E_{beam} [MeV]",180,0,180,37,200,1421);
openingangle_gammatogamma_energy_allcuts = new HistoManu2("openingangle_gammatogamma_energy_allcuts","openingangle_gammatogamma_energy_allcuts;#angle #gamma,#gamma;E_{beam} [MeV]",180,0,180,37,200,1421);

openingangle_gammatogamma_energy_allcuts_proton = new HistoManu2("openingangle_gammatogamma_energy_allcuts_proton","openingangle_gammatogamma_energy_allcuts;#angle #gamma,#gamma;E_{beam} [MeV]",180,0,180,37,200,1421);

openingangle_ptopi0_energy_nocuts = new HistoManu2("openingangle_ptopi0_energy_nocuts","openingangle_ptopi0_energy_nocuts;#angle #pi,p;E_{beam} [MeV]",180,0,180,37,200,1421);
openingangle_gammatogamma_energy_nocuts = new HistoManu2("openingangle_gammatogamma_energy_nocuts","openingangle_gammatogamma_energy_nocuts;#angle #gamma,#gamma;E_{beam} [MeV]",180,0,180,37,200,1421);


test = new TH1F("test","test",1100,-5.5,5.5);

coplanarityverteilung_collerated = new HistoManu2("coplanarityverteilung_collerated", "Coplanarity-Verteilung;#phi_{#pi}-#phi_{p}[deg]", 400,0,360,37,200,1421);

thetaverteilung_collerated = new HistoManu2("thetaverteilung_collerated", "Polar-Verteilung;#Delt #theta_{p}[deg]", 400,-200, 200,37,200,1421);

missingmassverteilung_collerated = new HistoManu2("missingmassverteilung_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,37,200,1421);
missingmassverteilung_forE_collerated = new HistoManu2("missingmassverteilung_forE_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,30, 220, 1420);

missingmassverteilung_collerated_proton = new HistoManu2("missingmassverteilung_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,37,200,1421);
missingmassverteilung_forE_collerated_proton = new HistoManu2("missingmassverteilung_forE_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,30, 220, 1420);

invmassverteilung_collerated = new HistoManu2("invmassverteilung_collerated", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,37,200,1421);

invmassverteilung_collerated_proton = new HistoManu2("invmassverteilung_collerated_proton", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,37,200,1421);

//3-dimensional cuts
coplanarityverteilung_cospi0_collerated = new HistoManu3("coplanarityverteilung_cospi0_collerated", "Coplanarity-Verteilung;#phi_{#pi}-#phi_{p}[deg]",400,0,360,72,-1,1,37,200,1421);
coplanarityverteilung_forE_collerated = new HistoManu2("coplanarityverteilung_forE_collerated", "Coplanarity-Verteilung;#phi_{#pi}-#phi_{p}[deg]", 400,0,360,30, 220, 1420);
clustersize_cospi0_energy_collerated = new HistoManu3("clustersize_cospi0_energy_collerated", "clustersize_cospi0_energy_collerated;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,37,200,1421);

clustersize_cospi0_energy_collerated_proton = new HistoManu3("clustersize_cospi0_energy_collerated_proton", "clustersize_cospi0_energy_collerated;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,37,200,1421);

thetaproton_cospi0_energy_collerated = new HistoManu3("thetaproton_cospi0_energy_collerated", "#theta_{p} [deg]",360,0,180,72,-1,1,37,200,1421);
thetaproton_phiproton_collerated = new HistoManu2("thetaproton_phiproton_collerated", "#theta_{p} [deg]",180,0,180,360,-180,180);

thetaproton_cospi0_energy_collerated_withoutcuts = new HistoManu3("thetaproton_cospi0_energy_collerated_withoutcuts", "#theta_{p} [deg]",360,0,180,72,-1,1,37,200,1421);
thetaproton_detector_collerated = new HistoManu2("thetaproton_detector_collerated", "#theta_{p} [deg];detector",360,0,180,40,0,45);

thetaverteilung_cospi0_collerated = new HistoManu3("thetaverteilung_cospi0_collerated", "Polar-Verteilung;#Delt #theta_{p}[deg]", 400,-200, 200,72,-1,1,37,200,1421);

missingmassverteilung_cospi0_collerated = new HistoManu3("missingmassverteilung_cospi0_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,72,-1,1,37,200,1421);

missingmassverteilung_cospi0_collerated_proton = new HistoManu3("missingmassverteilung_cospi0_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,72,-1,1,37,200,1421);

invmassverteilung_cospi0_collerated = new HistoManu3("invmassverteilung_cospi0_collerated", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,72,-1,1,37,200,1421);;

invmassverteilung_cospi0_collerated_proton = new HistoManu3("invmassverteilung_cospi0_collerated_proton", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,72,-1,1,37,200,1421);

//kinematic variables in dependence of energy
cosverteilung_collerated = new HistoManu2("cosverteilung_collerated", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_minus = new HistoManu2("cosverteilung_collerated_minus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_plus = new HistoManu2("cosverteilung_collerated_plus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);

cosverteilung_collerated_pbminus = new HistoManu2("cosverteilung_collerated_pbminus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_pbplus = new HistoManu2("cosverteilung_collerated_pbplus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_pt = new HistoManu2("cosverteilung_collerated_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_forE_collerated = new HistoManu2("cosverteilung_forE_collerated", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

// //kinematic variables in dependence of energy
cosverteilung_collerated_proton = new HistoManu2("cosverteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_proton_plus = new HistoManu2("cosverteilung_collerated_proton_plus ", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_proton_minus = new HistoManu2("cosverteilung_collerated_proton_minus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);

cosverteilung_collerated_proton_pt = new HistoManu2("cosverteilung_collerated_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_proton_pbplus = new HistoManu2("cosverteilung_collerated_proton_pbplus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_collerated_proton_pbminus = new HistoManu2("cosverteilung_collerated_proton_pbminus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);


cosverteilung_forE_collerated_proton = new HistoManu2("cosverteilung_forE_collerated_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);
// 
events_all = new TH3F("events_all","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);
events_witherror = new TH3F("events_witherror","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//helicity analysis for proton
cos_beam_hel1_targetplus_proton = new HistoManu2("cos_beam_hel1_targetplus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_proton = new HistoManu2("cos_beam_hel0_targetplus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_proton = new HistoManu2("cos_beam_hel1_targetminus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_proton = new HistoManu2("cos_beam_hel0_targetminus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetplus_proton_pc = new HistoManu2("cos_beam_hel1_targetplus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_proton_pc = new HistoManu2("cos_beam_hel0_targetplus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_proton_pc = new HistoManu2("cos_beam_hel1_targetminus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_proton_pc = new HistoManu2("cos_beam_hel0_targetminus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetplus_proton_pt = new HistoManu2("cos_beam_hel1_targetplus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_proton_pt = new HistoManu2("cos_beam_hel0_targetplus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_proton_pt = new HistoManu2("cos_beam_hel1_targetminus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_proton_pt = new HistoManu2("cos_beam_hel0_targetminus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

//sigma and g analysis for proton
kristallminus_targetplus_collerated_proton = new HistoManu3("kristallminus_targetplus_collerated_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kohlenstoffasymmetrie_collerated_proton = new HistoManu3("kohlenstoffasymmetrie_collerated_proton","kohlenstoffasymmetrie_collerated_proton; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_proton = new HistoManu3("kristallminus_targetminus_collerated_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_proton = new HistoManu3("kristallplus_targetplus_collerated_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_proton = new HistoManu3("kristallplus_targetminus_collerated_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb_proton = new HistoManu3("kristallminus_targetplus_collerated_pb_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pb_proton = new HistoManu3("kristallminus_targetminus_collerated_pb_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pb_proton = new HistoManu3("kristallplus_targetplus_collerated_pb_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pb_proton = new HistoManu3("kristallplus_targetminus_collerated_pb_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt_proton = new HistoManu3("kristallminus_targetplus_collerated_pt_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pt_proton = new HistoManu3("kristallminus_targetminus_collerated_pt_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pt_proton = new HistoManu3("kristallplus_targetplus_collerated_pt_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pt_proton = new HistoManu3("kristallplus_targetminus_collerated_pt_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//helicity analysis without proton 

cos_beam_hel1_targetplus = new HistoManu2("cos_beam_hel1_targetplus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus = new HistoManu2("cos_beam_hel0_targetplus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus = new HistoManu2("cos_beam_hel1_targetminus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus = new HistoManu2("cos_beam_hel0_targetminus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetplus_pc = new HistoManu2("cos_beam_hel1_targetplus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_pc = new HistoManu2("cos_beam_hel0_targetplus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_pc = new HistoManu2("cos_beam_hel1_targetminus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_pc = new HistoManu2("cos_beam_hel0_targetminus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetplus_pt = new HistoManu2("cos_beam_hel1_targetplus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_pt = new HistoManu2("cos_beam_hel0_targetplus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_pt = new HistoManu2("cos_beam_hel1_targetminus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_pt = new HistoManu2("cos_beam_hel0_targetminus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

//analysis for sigma and g without proton

kohlenstoffasymmetrie_collerated = new HistoManu3("kohlenstoffasymmetrie_collerated","kohlenstoffasymmetrie_collerated; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetplus_collerated = new HistoManu3("kristallminus_targetplus_collerated","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated = new HistoManu3("kristallminus_targetminus_collerated","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated = new HistoManu3("kristallplus_targetplus_collerated","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated = new HistoManu3("kristallplus_targetminus_collerated","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb = new HistoManu3("kristallminus_targetplus_collerated_pb","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pb = new HistoManu3("kristallminus_targetminus_collerated_pb","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pb = new HistoManu3("kristallplus_targetplus_collerated_pb","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pb = new HistoManu3("kristallplus_targetminus_collerated_pb","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt = new HistoManu3("kristallminus_targetplus_collerated_pt","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pt = new HistoManu3("kristallminus_targetminus_collerated_pt","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pt = new HistoManu3("kristallplus_targetplus_collerated_pt","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);
kristallplus_targetminus_collerated_pt = new HistoManu3("kristallplus_targetminus_collerated_pt","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);


beamenergy_gen  = new HistoManu("beamenergy_gen ", "beamenergy_gen ",  300, 0, 1557);
beamenergy_gen->Sumw2();

//for monte carlo studies

cosmeson_beamenergy_energysum_rek = new HistoManu3("cosmeson_beamenergy_energysum_rek", "cosmeson_beamenergy_energysum_rek;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1, 30, 220, 1420,200, 0, 2000);
cosmeson_beamenergy_energysum_rek->Sumw2();

cosmeson_beamenergy_energysum_rek_proton = new HistoManu3("cosmeson_beamenergy_energysum_rek_proton", "cosmeson_beamenergy_energysum_rek_proton;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",18, -1, 1,  30, 220, 1420,200, 0, 2000);
cosmeson_beamenergy_energysum_rek_proton->Sumw2();

cosmeson_beamenergy_energysum_monte = new HistoManu3("cosmeson_beamenergy_energysum_monte", "cosmeson_beamenergy_energysum_monte;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);
cosmeson_beamenergy_energysum_monte->Sumw2();

cosmeson_beamenergy_monte = new HistoManu2("cosmeson_beamenergy_monte", "cosmeson_beamenergy_monte;E_{#gamma} [MeV];cos #theta_{#pi}", 18, -1, 1, 30, 220, 1420);
cosmeson_beamenergy_monte->Sumw2();

cosmeson_beamenergy_rek_3ped = new HistoManu2("cosmeson_beamenergy_rek_3ped", "cosmeson_beamenergy_rek_3ped;E_{#gamma} [MeV];cos #theta_{#pi}", 18, -1, 1, 30, 220, 1420);
cosmeson_beamenergy_rek_3ped->Sumw2();

cosmeson_beamenergy_rek_2_3ped = new HistoManu2("cosmeson_beamenergy_rek_2_3ped", "cosmeson_beamenergy_rek_2_3ped;E_{#gamma} [MeV];cos #theta_{#pi}",  18, -1, 1,30, 220, 1420);
cosmeson_beamenergy_rek_2_3ped->Sumw2();

h_energy_sum_pi0_3ped = new HistoManu("h_energy_sum_pi0_3ped ", "h_energy_sum_pi0_3ped;E_{sum} [MeV]",  300, 0, 1557);

h_energy_sum_pi0_2_3ped = new HistoManu("h_energy_sum_pi0_2_3ped ", "h_energy_sum_pi0_3ped;E_{sum} [MeV]",  300, 0, 1557);

triggertest = new TH1F("triggertest","triggerpattern",34,0,34);

pionmasse = 134.9766;
unten_copl = -25.;
oben_copl = 25.;
unten_theta=-10.;
oben_theta=10.;

//DEFINE TCUTG for dE_E!

   cuttaps->SetVarX("dE_E (all TAPS clusters compared to Veto hits)");
   cuttaps->SetVarY("");
   cuttaps->SetTitle("Graph");
   cuttaps->SetTitle("CutProton");
   cuttaps->SetFillColor(1);
   cuttaps->SetPoint(0,1.84744,9.4837-0.5);
   cuttaps->SetPoint(1,1.84744,7.3947-0.5);
   cuttaps->SetPoint(2,69.1895,4.09986-1);
   cuttaps->SetPoint(3,197.318,2.40149-0.5);
   cuttaps->SetPoint(4,299.225,2.43546-0.5);
   cuttaps->SetPoint(5,399.344,2.43546-0.5);
   cuttaps->SetPoint(6,399.344,4.01495-0.5);
   cuttaps->SetPoint(7,227.712,4.84715-0.5);
   cuttaps->SetPoint(8,121.633,6.52853-0.5);
   cuttaps->SetPoint(9,85.876,8.10802-0.5);
   cuttaps->SetPoint(10,71.5733,9.4837-0.5);
   cuttaps->SetPoint(11,1.84744,9.4837-0.5);
   cuttaps->Draw("");

   cutcb->SetVarX("dE_E (all CB clusters compared to PID hits)");
   cutcb->SetVarY("");
   cutcb->SetTitle("Graph");
   cutcb->SetFillColor(1);
   cutcb->SetPoint(0,13.2223,4.30711);
   cutcb->SetPoint(1,113.154,2.30543);
   cutcb->SetPoint(2,397.995,1.95549);
   cutcb->SetPoint(3,400.714,3.15929);
   cutcb->SetPoint(4,278.688,3.43925);
   cutcb->SetPoint(5,162.101,5.66489);
   cutcb->SetPoint(6,116.893,7.41461);
   cutcb->SetPoint(7,78.8239,9.9902);
   cutcb->SetPoint(8,14.242,9.9902);
   cutcb->SetPoint(9,14.242,4.27912);
   cutcb->SetPoint(10,12.5425,4.27912);
   cutcb->SetPoint(11,13.2223,4.30711);
HistoManu::InitCuts(-8, 8,100, 300);




}


PPi0Analyse::~PPi0Analyse()
{
}



Bool_t	PPi0Analyse::Start()
{
//targetpol
pt=PPi0Analyse::targetpol(inputFile);

//filenumber
TString* filename1 = new TString(inputFile->GetPath());
path11= new TString(filename1->Tokenize("_")->At(filename1->Tokenize("_")->GetEntries()-1)->GetName());
path11->Resize(path11->Length()-7);

//runnumber
stringstream ss(path11->Data());
ss >> runNumber;
if(runNumber>=1363 && runNumber<=3374)planesetting=PPi0Analyse::polplane(inputFile);
else{
if(GetLinpol()->GetPolarizationPlane()==1){planesetting="PERP";}
if(GetLinpol()->GetPolarizationPlane()==0){planesetting="PARA";}

calfile.open(Form("/disk/user/spieker/GoAT/txtfiles/1/liste_PID_May14_old_prompt_%i.txt",runNumber));
calfile1.open(Form("/disk/user/spieker/GoAT/txtfiles/1/liste_PID_May14_old_side_%i.txt",runNumber));


}

TFile *g = new TFile("/disk/user/spieker/Mainz/cutranges_3sigma_mainz_withMWPC_nov13may14_test2.root", "OPEN");

//Copl
cut_copldown = (TH2F*)g->Get("Selection/copl_down1");
cut_coplup = (TH2F*)g->Get("Selection/copl_up1");

//Theta
cut_thetadown = (TH2F*)g->Get(Form("Selection/theta_down1"));
cut_thetaup = (TH2F*)g->Get(Form("Selection/theta_up1"));

//Missing Mass
cut_mmdown = (TH2F*)g->Get(Form("Selection/mm_down1"));
cut_mmup = (TH2F*)g->Get(Form("Selection/mm_up1"));

//inv mass
cut_invdown = (TH2F*)g->Get(Form("Selection/inv_down1"));
cut_invup = (TH2F*)g->Get(Form("Selection/inv_up1"));

    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();


return kTRUE;
}


void	PPi0Analyse::ProcessEvent()	
{
// cout << GetTracks()->GetTAPSDistance() << endl;

	for(Int_t i=0;i<GetTrigger()->GetNTriggerPattern();i++){	
	triggertest->Fill(GetTrigger()->GetTriggerPattern(i));
	}
	
	Bool_t cbtrigger=kFALSE, cbtriggercond=kFALSE, tapstrigger=kFALSE;
	
		for(Int_t i=0;i<GetTrigger()->GetNTriggerPattern();i++){
			if(GetTrigger()->GetTriggerPattern(i)==16){
				cbtrigger=kTRUE;
				if(GetTrigger()->GetEnergySum()>100.){
				cbtriggercond=kTRUE;
				}
			}
		}

// 	if(GetScalers()->GetNEntries()==0){//MC has no GetScalers()
// 
// 		TLorentzVector beam_mc = GetGeant()->GetBeam();
// 		Double_t e_beam_mc = beam_mc.E()*1000.;
// 		Double_t time1=1.0;	
// 		beamenergy_gen->Fill(e_beam_mc,time1);
// 	
// 		TLorentzVector proton_mc = GetGeant()->GetTrueVector(0);
// 		proton_mc.SetPxPyPzE(proton_mc.Px()*1000,proton_mc.Py()*1000,proton_mc.Pz()*1000,proton_mc.E()*1000);
// 		
// 		TLorentzVector g1_mc = GetGeant()->GetTrueVector(1);
// 		TLorentzVector g2_mc = GetGeant()->GetTrueVector(2);
// 		TLorentzVector meson_mc = g1_mc+g2_mc;
// 		meson_mc.SetPxPyPzE(meson_mc.Px()*1000,meson_mc.Py()*1000,meson_mc.Pz()*1000,meson_mc.E()*1000);
// 
// 		//Boost-System
// 		TLorentzVector meson_mc_cms = CMVector(meson_mc, beam_mc, proton_mc);
// 		Double_t cosT_mc_cms= meson_mc_cms.CosTheta();
// 		if(e_beam_mc>=200&& e_beam_mc<=1420){
// 			cosmeson_beamenergy_energysum_monte->Fill(cosT_mc_cms,e_beam_mc,1000*GetGeant()->GetCBESum(),time1);
// 			cosmeson_beamenergy_monte->Fill(cosT_mc_cms,e_beam_mc,time1);
// 		}
// 	}

	Double_t anzahlpi0_prompt=0;
	Double_t anzahlpi0_sideband=0;
	Bool_t cutsfullfilled;
	TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
	TLorentzVector  beam_4vect;
	TLorentzVector  missingp_4vect;
	TLorentzVector pi0_4vect;

	if(GetPhotons()->GetNParticles()==2){
		pi0_4vect = GetPhotons()->Particle(0)+GetPhotons()->Particle(1);
		//invMass
		Double_t inv=pi0_4vect.M();
	
		//TEST the target polarization (no entry with pt=0)
		test->Fill(pt);

		for(Int_t j=0; j < GetTagger()->GetNTagged();j++)
		{

// 			if(GetPhotons()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;
// 			if(GetPhotons()->GetDetectors(1)==GTreeTrack::DETECTOR_PbWO4) continue;
// 			if(GetRootinos()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;
// 			cout << runNumber << endl;
 			if(runNumber>=1363 && runNumber<=3374 && pt!=5){continue;}//only carbon runs
 			if((runNumber<1363 || runNumber>3374) && pt==5){continue;}//only moeller/diamond runs

			Double_t time= GetTagger()->GetTaggedTime(j) - 0.5*(GetPhotons()->GetTime(1)+GetPhotons()->GetTime(0));
			time1->Fill(time);

			if(HistoManu::IsPrompt(time)){
				time_prompt->Fill(time);
			}

			if(HistoManu::IsRandom(time)){	
				time_side->Fill(time);
			}

			Double_t timediff_photons=  GetPhotons()->GetTime(0)- GetPhotons()->GetTime(1);


 			if(GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j))==-1){continue;}//skip events with wrong polarisation!

			if(!(HistoManu::IsRandom(time) || HistoManu::IsPrompt(time))){continue;} //skip events which are not prompt or random

			poltable_energy->Fill(GetTagger()->GetTaggedEnergy(j));
			poltable_energy_weight->Fill(GetTagger()->GetTaggedEnergy(j),GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j)));

			//get target, beam and missng particle information
			beam_4vect = TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(j),GetTagger()->GetTaggedEnergy(j));
			missingp_4vect = beam_4vect + protonvektor_target - pi0_4vect;
			Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
			Double_t beamphoton1E=GetTagger()->GetTaggedEnergy(j);
			Double_t pb=GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j));
			Double_t circpol = GetCircPol(beamphoton1E,runNumber);
			Bool_t helicity =GetTrigger()->GetHelicity();

			//due to kinematic not possible -> kick them out
			if(anglethetaproton_rek > 70){continue;}
			
			if(beamphoton1E<200 || beamphoton1E>1430){continue;}

			//Missing Mass
			Double_t missingmass=missingp_4vect.M();

			//Boost-System
			TLorentzVector pi0_4vect_boost = CMVector(pi0_4vect, beam_4vect, protonvektor_target);
			TVector3 pi0_3vektor_boost = pi0_4vect_boost.Vect();
			Double_t cospi0 = pi0_3vektor_boost.CosTheta();
			Double_t phimeson = TMath::RadToDeg()*(pi0_4vect.Vect().Phi());
			Double_t phimeson_CMS = TMath::RadToDeg()*(pi0_4vect_boost.Vect().Phi());
			Double_t thetameson_CMS = TMath::RadToDeg()*(pi0_4vect_boost.Vect().Theta());

			//energy dependent cuts
			/*Double_t oben_copl=233.412+(-0.13196)*beamphoton1E+(0.000137884)*beamphoton1E*beamphoton1E+(-4.70835e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_copl=122.821+(0.136561)*beamphoton1E+(-0.000140207)*beamphoton1E*beamphoton1E+(4.74016e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_theta=-7.49325+(0.0624572)*beamphoton1E+(-5.35677e-05)*beamphoton1E*beamphoton1E+(9.25958e-09)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_theta=-8.26291+(-0.000512115)*beamphoton1E+(-1.12875e-05)*beamphoton1E*beamphoton1E+(8.75687e-09)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_mass=923.325+(0.264511)*beamphoton1E+(-0.000372629)*beamphoton1E*beamphoton1E+(2.10017e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_mass=947.805+(-0.182183)*beamphoton1E+(7.12484e-05)*beamphoton1E*beamphoton1E+(-4.50161e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_mass_proton=962.88+(0.0458514)*beamphoton1E+(3.80859e-05)*beamphoton1E*beamphoton1E+(-9.62237e-09)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_mass_proton=918.615+(-0.00536361)*beamphoton1E+(-0.000256642)*beamphoton1E*beamphoton1E+(1.20651e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_inv=151.025+(0.0490305)*beamphoton1E+(-5.08058e-05)*beamphoton1E*beamphoton1E+(2.07695e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_inv=109.144+(0.00456935)*beamphoton1E+(-9.11259e-06)*beamphoton1E*beamphoton1E+(-3.0289e-09)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_inv_proton=159.356+(-0.0111137)*beamphoton1E+(5.53727e-05)*beamphoton1E*beamphoton1E+(-3.00086e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_inv_proton=102.522+(0.0582579)*beamphoton1E+(-0.000110087)*beamphoton1E*beamphoton1E+(4.95622e-08)*beamphoton1E*beamphoton1E*beamphoton1E;*/

		/*	//energy independent cuts

			Double_t oben_copl=190+20;
			Double_t unten_copl=190-20;
			Double_t oben_theta=20;
			Double_t unten_theta=-20;
			Double_t oben_mass=938+77;
			Double_t unten_mass=938-77;
			Double_t oben_inv=135+25;
			Double_t unten_inv=135-25;
			Double_t oben_mass_proton=938+77;
			Double_t unten_mass_proton=938-77;
			Double_t oben_inv_proton=135+25;
			Double_t unten_inv_proton=135-25;
*/
			//three dimensional cut
 			Double_t unten_copl=cut_copldown->GetBinContent(cut_copldown->GetXaxis()->FindBin(cospi0),cut_copldown->GetYaxis()->FindBin(beamphoton1E));
 			Double_t oben_copl=cut_coplup->GetBinContent(cut_coplup->GetXaxis()->FindBin(cospi0),cut_coplup->GetYaxis()->FindBin(beamphoton1E));
 
 			Double_t unten_inv=cut_invdown->GetBinContent(cut_invdown->GetXaxis()->FindBin(cospi0),cut_invdown->GetYaxis()->FindBin(beamphoton1E));
 			Double_t oben_inv=cut_invup->GetBinContent(cut_invup->GetXaxis()->FindBin(cospi0),cut_invup->GetYaxis()->FindBin(beamphoton1E));
 			Double_t unten_inv_proton=unten_inv;
 			Double_t oben_inv_proton=oben_inv;
 
 			Double_t unten_mass=cut_mmdown->GetBinContent(cut_mmdown->GetXaxis()->FindBin(cospi0),cut_mmdown->GetYaxis()->FindBin(beamphoton1E));
 			Double_t oben_mass=cut_mmup->GetBinContent(cut_mmup->GetXaxis()->FindBin(cospi0),cut_mmup->GetYaxis()->FindBin(beamphoton1E));
 
 			Double_t unten_mass_proton=unten_mass;
 			Double_t oben_mass_proton=oben_mass;
 
 			Double_t unten_theta=cut_thetadown->GetBinContent(cut_thetadown->GetXaxis()->FindBin(cospi0),cut_thetadown->GetYaxis()->FindBin(beamphoton1E));
 			Double_t oben_theta=cut_thetaup->GetBinContent(cut_thetaup->GetXaxis()->FindBin(cospi0),cut_thetaup->GetYaxis()->FindBin(beamphoton1E));

//_________________________________only one charged particle__________________________________________________________________________

			Double_t anglephimeson;
			Double_t anglephiproton;
			Double_t phi=1000000000000;
			Double_t anglethetaproton_meas,theta_proton_CMS,phi_pmeasured,theta_pmeasured,Ekin_pmissing,phi_proton_CMS;
			Double_t thetadiff=1000000000000;
			Double_t thetadiff_CMS=1000000000000;

			Double_t opening_angle_gammatogamma=TMath::ACos(TMath::Cos(GetPhotons()->Particle(0).Phi()-(GetPhotons()->Particle(1)).Phi())*TMath::Sin(GetPhotons()->Particle(0).Theta())*TMath::Sin((GetPhotons()->Particle(1)).Theta())+ TMath::Cos(GetPhotons()->Particle(0).Theta())*TMath::Cos((GetPhotons()->Particle(1)).Theta()) )*TMath::RadToDeg();

			if(GetRootinos()->GetNParticles()==1){

				if(GetRootinos()->GetClusterSize(0)>4){continue;}//skip clusters with smaller than 4

				//get charged particle information__________________________________
				TLorentzVector proton_4vect_meas;
				if(GetRootinos()->GetNParticles()==1){proton_4vect_meas=GetRootinos()->Particle(0);}

				//new Lorentzvector for missing proton
				Ekin_pmissing=TMath::Sqrt(TMath::Power(missingp_4vect.E(),2)-TMath::Power(938.272046,2));
				phi_pmeasured=GetRootinos()->GetPhi(0);
				theta_pmeasured=GetRootinos()->GetTheta(0);
				TLorentzVector proton_4vect_constructed=TLorentzVector(Ekin_pmissing*TMath::Sin(GetRootinos()->GetThetaRad(0))*TMath::Cos(GetRootinos()->GetPhiRad(0)), Ekin_pmissing*TMath::Sin(GetRootinos()->GetThetaRad(0))*TMath::Sin(GetRootinos()->GetPhiRad(0)), Ekin_pmissing*TMath::Cos(GetRootinos()->GetThetaRad(0)), missingp_4vect.E());

				//Boost-System of Proton
				TLorentzVector proton_4vect_constructed_boost = CMVector(proton_4vect_constructed, beam_4vect, protonvektor_target);
				theta_proton_CMS=TMath::RadToDeg()*proton_4vect_constructed_boost.Vect().Theta();
				phi_proton_CMS=TMath::RadToDeg()*proton_4vect_constructed_boost.Vect().Phi();


				if(phi_proton_CMS<0. && phi_proton_CMS> -180.) theta_proton_CMS=  theta_proton_CMS*-1.;
				if(phimeson_CMS<0. && phimeson_CMS> -180.) thetameson_CMS=  thetameson_CMS*-1.;

				//Theta-Difference in CMS system
				thetadiff_CMS=theta_proton_CMS-thetameson_CMS;
				if(thetadiff_CMS < 0){thetadiff_CMS = thetadiff_CMS + 360;}


				//apply TCUTGs
				if(GetRootinos()->HasCB(0) && cutcb->IsInside(proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0))==1){kCB=10;}else{kCB=10;}
				if(GetRootinos()->HasTAPS(0) && cuttaps->IsInside(proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0))==1){kTAPS=10;}else{kTAPS=10;}

				Double_t opening_angle_ptopi0= TMath::ACos(TMath::Cos(proton_4vect_meas.Phi()-(pi0_4vect).Phi())*TMath::Sin(proton_4vect_meas.Theta())*TMath::Sin((pi0_4vect).Theta())+ TMath::Cos(proton_4vect_meas.Theta())*TMath::Cos((pi0_4vect).Theta()) )*TMath::RadToDeg();

				Double_t timediff_pmeson=GetRootinos()->GetTime(0) - 0.5*(GetPhotons()->GetTime(0)+GetPhotons()->GetTime(1));

				//Phi-Difference
				anglephimeson = TMath::RadToDeg()*pi0_4vect.Vect().Phi();
				anglephiproton = phi_pmeasured;
				phi = anglephimeson - anglephiproton;
				if(phi < 0){phi = phi + 360;}
			
				//theta differenz
				anglethetaproton_meas = theta_pmeasured;
				thetadiff=anglethetaproton_rek-anglethetaproton_meas;
// 				cout << anglethetaproton_meas << "\t" << anglephiproton << endl;

				thetaproton_cospi0_energy_collerated_withoutcuts->Fill(anglethetaproton_meas,cospi0,beamphoton1E,time);

				//due to kinematic not possible -> kick them out
				if(anglethetaproton_meas > 70){continue;}

				//without any cuts!
				openingangle_ptopi0_energy_nocuts->Fill(opening_angle_ptopi0,beamphoton1E,time);
				openingangle_gammatogamma_energy_nocuts->Fill(opening_angle_gammatogamma,beamphoton1E,time);
			
				//visualisation of TCutG
				if(GetRootinos()->GetDetectors(0)==7)Check_CBdE_E_nocuts_all->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0)*TMath::Sin(GetRootinos()->GetTheta(0)),time);
		
				if(GetRootinos()->GetDetectors(0)==3)Check_CBdE_E_nocuts_CBandPID->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0),time);
		
				if(GetRootinos()->GetDetectors(0)==5)Check_CBdE_E_nocuts_CBandMWPC->Fill(GetRootinos()->Particle(0).E(),0.5*(GetRootinos()->GetMWPC1Energy(0)+GetRootinos()->GetMWPC0Energy(0))*TMath::Sin(GetRootinos()->GetTheta(0)),time);
		
				if(GetRootinos()->HasTAPS(0)){
					Check_TAPSdE_E_nocuts->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0),time);
				}

				MM_proton->Fill(missingmass,time);
				theta_proton->Fill(thetadiff,time);
				coplanarity_proton->Fill(phi,time);
				IM_proton->Fill(inv,time);

    				coplanarity_theta_CMS->Fill(thetadiff_CMS,beamphoton1E,time);

				//visualisation of several cuts!

				if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > 	(unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

					theta_all_proton->Fill(thetadiff,time);
					thetaverteilung_collerated->Fill(thetadiff,beamphoton1E,time);
					thetaverteilung_cospi0_collerated->Fill(thetadiff,cospi0,beamphoton1E,time);
				}
	
				if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (thetadiff > (unten_theta) && thetadiff < (oben_theta)) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

					coplanarity_all_proton->Fill(phi,time);
					coplanarityverteilung_collerated->Fill(phi,beamphoton1E,time);
					if(helicity==0 || helicity==1){coplanarityverteilung_forE_collerated->Fill(phi,beamphoton1E,time);}
					coplanarityverteilung_cospi0_collerated->Fill(phi,cospi0,beamphoton1E,time);
				}

				if((thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

					MM_all_proton->Fill(missingmass,time);
					missingmassverteilung_collerated_proton->Fill(missingmass,beamphoton1E,time);
					if(helicity==1 || helicity==0){missingmassverteilung_forE_collerated_proton->Fill(missingmass,beamphoton1E,time);}
					missingmassverteilung_cospi0_collerated_proton->Fill(missingmass,cospi0,beamphoton1E,time);
					}

				if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) &&  (phi > ((unten_copl)) && phi < ((oben_copl)))&&(kTAPS==10 || kCB==10)){
	
					IM_all_proton->Fill(inv,time);
					invmassverteilung_collerated_proton->Fill(inv,beamphoton1E,time);
					invmassverteilung_cospi0_collerated_proton->Fill(inv,cospi0,beamphoton1E,time);			
				}


				//stuff after all cuts and creation of needed analysis histograms_______________________________
	
				if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)&&thetadiff_CMS>160 && thetadiff_CMS<200){

    					coplanarity_theta_CMS_all->Fill(thetadiff_CMS,beamphoton1E,time);


					//PULSE SHAPE ANALYSIS
					for(Int_t i=0; i<GetTracks()->GetNTracks() ;i++){//loop over all tracks

						Double_t clusterenergy = GetTracks()->GetClusterEnergy(i);
						
						Double_t PSA_angle=GetTracks()->GetPSAAngle(i);
						Double_t PSA_radius=GetTracks()->GetPSARadius(i);			

						if(GetPhotons()->HasTAPS(0) && clusterenergy== (GetPhotons()->Particle(0)).E() ){
							Check_PSA_3ped->Fill(PSA_angle,PSA_radius,time);
						}

						if(GetPhotons()->HasTAPS(1) && clusterenergy== (GetPhotons()->Particle(1)).E()){	
							Check_PSA_3ped->Fill(PSA_angle,PSA_radius,time);
						}

						if(GetRootinos()->HasTAPS(0) && clusterenergy== (GetRootinos()->Particle(0)).E()){
							Check_PSA_3ped->Fill(PSA_angle,PSA_radius,time);
						}

					}
					//NUMBER OF PIOs
					if(HistoManu::IsPrompt(time)){anzahlpi0_prompt++; cutsfullfilled=kTRUE;}
					if(HistoManu::IsRandom(time)){anzahlpi0_sideband++;}


					//visualisation of TCutG
					if(GetRootinos()->GetDetectors(0)==7)  Check_CBdE_E_all->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0)*TMath::Sin(GetRootinos()->GetTheta(0)),time);
			
					if(GetRootinos()->GetDetectors(0)==3)  Check_CBdE_E_CBandPID->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0),time);
			
					if(GetRootinos()->GetDetectors(0)==5)  Check_CBdE_E_CBandMWPC->Fill(GetRootinos()->Particle(0).E(),0.5*(GetRootinos()->GetMWPC1Energy(0)+GetRootinos()->GetMWPC0Energy(0))*TMath::Sin(GetRootinos()->GetTheta(0)),time);
			
					if(GetRootinos()->HasTAPS(0)){
						Check_TAPSdE_E->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0),time);
					}


					if(GetRootinos()->HasTAPS(0)){
						Check_TAPS_TOF_proton->Fill(((GetRootinos()->GetTime(0)-GetTagger()->GetTaggedTime(j))/(1.75700/TMath::Cos(GetRootinos()->GetThetaRad(0))))+1/0.299792458,proton_4vect_meas.E(),time);
					}
					
					//TAPS TOF 
					if(GetPhotons()->HasTAPS(0)){
						Check_TAPS_TOF_photon_3ped->Fill(((GetPhotons()->GetTime(0)-GetTagger()->GetTaggedTime(j))/(1.75700/TMath::Cos(GetPhotons()->GetThetaRad(0))))+1/0.299792458,GetPhotons()->Particle(0).E(),time);
					}

					if(GetPhotons()->HasTAPS(1)){
						Check_TAPS_TOF_photon_3ped->Fill(((GetPhotons()->GetTime(1)-GetTagger()->GetTaggedTime(j))/(1.75700/TMath::Cos(GetPhotons()->GetThetaRad(1))))+1/0.299792458,GetPhotons()->Particle(1).E(),time);
					}
					
					//some check histograms concerning proton angle and clustersizes
					thetaproton_cospi0_energy_collerated->Fill(anglethetaproton_meas,cospi0,beamphoton1E,time);
// 					if(GetRootinos()->GetDetectors(0)==40 && HistoManu::IsPrompt(time))cout << GetRootinos()->GetDetectors(0) << "\t" << time << "\t"<< anglethetaproton_meas <<endl;


// 					if(GetRootinos()->HasTAPS(0))anglethetaproton_meas= (TMath::RadToDeg()*TMath::ATan2((145.700*TMath::Tan(anglethetaproton_meas*TMath::DegToRad())),(175.700)));

					thetaproton_detector_collerated->Fill(anglethetaproton_meas,GetRootinos()->GetDetectors(0),time);
					thetaproton_phiproton_collerated->Fill(anglethetaproton_meas,anglephiproton*TMath::RadToDeg(),time);
					clustersize_cospi0_energy_collerated->Fill(GetPhotons()->GetClusterSize(0),cospi0,beamphoton1E,time);
					clustersize_cospi0_energy_collerated->Fill(GetPhotons()->GetClusterSize(1),cospi0,beamphoton1E,time);
					clustersize_cospi0_energy_collerated_proton->Fill(GetRootinos()->GetClusterSize(0),cospi0,beamphoton1E,time);

     					if(HistoManu::IsPrompt(time)) calfile << anglethetaproton_meas << "\t" << GetEventParameters()->GetEventNumber() << "\t" << runNumber << "\t" <<GetRootinos()->GetDetectors(0)<< "\n";
     					if(HistoManu::IsRandom(time)) calfile1 << anglethetaproton_meas<< "\t" << GetEventParameters()->GetEventNumber() << "\t" << runNumber << "\t" << GetRootinos()->GetDetectors(0) <<"\n";

					cosverteilung_collerated_proton->Fill(cospi0,beamphoton1E,time);
					if(planesetting=="PERP"){

						cosverteilung_collerated_proton_pbplus->Fill(cospi0,beamphoton1E,time,pb);
						cosverteilung_collerated_proton_plus->Fill(cospi0,beamphoton1E,time);


					}
					if(planesetting=="PARA"){

						cosverteilung_collerated_proton_pbminus->Fill(cospi0,beamphoton1E,time,pb);
						cosverteilung_collerated_proton_minus->Fill(cospi0,beamphoton1E,time);

					}

					if(helicity==1 || helicity==0){cosverteilung_forE_collerated_proton->Fill(cospi0,beamphoton1E,time);}

// 					if(GetScalers()->GetNEntries()==0){
// 
// 						h_energy_sum_pi0_3ped->Fill(GetTrigger()->GetEnergySum(),time);
// 						cosmeson_beamenergy_energysum_rek_proton->Fill(cospi0,beamphoton1E,GetTrigger()->GetEnergySum(),time);
// 						cosmeson_beamenergy_rek_3ped->Fill(cospi0,beamphoton1E,time);
// 					}
						
					openingangle_ptopi0_energy_allcuts->Fill(opening_angle_ptopi0,beamphoton1E,time);
					openingangle_gammatogamma_energy_allcuts_proton->Fill(opening_angle_gammatogamma,beamphoton1E,time);


					//starting histogram for extraction of G and E
					if(pt>0 ){

						cosverteilung_collerated_proton_pt->Fill(cospi0,beamphoton1E,time,pt);

						if(helicity==1){
							cos_beam_hel1_targetplus_proton->Fill(cospi0,beamphoton1E,time);
							cos_beam_hel1_targetplus_proton_pc->Fill(cospi0,beamphoton1E,time,circpol);
							cos_beam_hel1_targetplus_proton_pt->Fill(cospi0,beamphoton1E,time,pt);
						}

						if(helicity==0){
							cos_beam_hel0_targetplus_proton->Fill(cospi0,beamphoton1E,time);
							cos_beam_hel0_targetplus_proton_pc->Fill(cospi0,beamphoton1E,time,circpol);
							cos_beam_hel0_targetplus_proton_pt->Fill(cospi0,beamphoton1E,time,pt);
						}
		
						if(planesetting=="PARA"){
							kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,time);
							kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,time,pb);
							kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,time,pt);
						}
	
						if(planesetting=="PERP"){
							kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,time);
							kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,time,pb);
							kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,time,pt);	
						}
					}

					if(pt<0 ){

						cosverteilung_collerated_proton_pt->Fill(cospi0,beamphoton1E,time,-1*pt);

						cosT=cospi0;
						Eg=beamphoton1E;
						Phi=phimeson;
	
						if(HistoManu::IsPrompt(time)){isSg=1;isBg=0;}else{isBg=1;isSg=0;}
						if(planesetting=="PARA"){isPARA=1;isPERP=0;}
						if(planesetting=="PERP"){isPARA=0;isPERP=1;}
	
// 						calfile << cosT << "\t" << Eg<< "\t" <<Phi<< "\t" << isSg << "\t" <<  isBg << "\t" << isPARA << "\t" << isPERP <<"\n";


						if(helicity==1){
							cos_beam_hel1_targetminus_proton->Fill(cospi0,beamphoton1E,time);
							cos_beam_hel1_targetminus_proton_pc->Fill(cospi0,beamphoton1E,time,circpol);
							cos_beam_hel1_targetminus_proton_pt->Fill(cospi0,beamphoton1E,time,pt);
						}

						if(helicity==0){
							cos_beam_hel0_targetminus_proton->Fill(cospi0,beamphoton1E,time);
							cos_beam_hel0_targetminus_proton_pc->Fill(cospi0,beamphoton1E,time,circpol);
							cos_beam_hel0_targetminus_proton_pt->Fill(cospi0,beamphoton1E,time,pt);
						}

						if(planesetting=="PARA"){
							kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,time);
							kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,time,pb);
							kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,time,pt);
						}
	
						if(planesetting=="PERP"){

							kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,time);
							kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,time,pb);
							kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,time,pt);
				
						}
	
					}
				}						
// 				
			}
// 				//all events with proton identified



			//without any cuts!

			MM->Fill(missingmass,time);
			IM->Fill(inv,time);


			if(GetRootinos()->GetNParticles()==1){

				if(!(thetadiff > (unten_theta) && thetadiff < (oben_theta)) ||  !(phi > (unten_copl) && phi < (oben_copl))){continue;} 

			}


			//with cuts
		
			if(((inv > (unten_inv) && inv < (oben_inv)))){

				MM_all->Fill(missingmass,time);
				missingmassverteilung_collerated->Fill(missingmass,beamphoton1E,time);
				if(helicity==1 || helicity==0){missingmassverteilung_forE_collerated->Fill(missingmass,beamphoton1E,time);}
				missingmassverteilung_cospi0_collerated->Fill(missingmass,cospi0,beamphoton1E,time);
			}
	
			if(((missingmass > (unten_mass) && missingmass < (oben_mass)))){

				IM_all->Fill(inv,time);
				invmassverteilung_collerated->Fill(inv,beamphoton1E,time);
				invmassverteilung_cospi0_collerated->Fill(inv,cospi0,beamphoton1E,time);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))&&GetTrigger()->GetNErrors()!=0){
				events_witherror->Fill(phimeson,cospi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){
				events_all->Fill(phimeson,cospi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){

	
				if(GetPhotons()->HasTAPS(0)){
					Check_TAPS_TOF_photon_2_3ped->Fill(((GetPhotons()->GetTime(0)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(0).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(0).E(),time);
				}

				if(GetPhotons()->HasTAPS(1)){
					Check_TAPS_TOF_photon_2_3ped->Fill(((GetPhotons()->GetTime(1)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(1).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(1).E(),time);
				}
		
				cosverteilung_collerated->Fill(cospi0,beamphoton1E,time);

				if(planesetting=="PERP"){

					cosverteilung_collerated_pbplus->Fill(cospi0,beamphoton1E,time,pb);
					cosverteilung_collerated_plus->Fill(cospi0,beamphoton1E,time);


				}
				if(planesetting=="PARA"){

					cosverteilung_collerated_pbminus->Fill(cospi0,beamphoton1E,time,pb);
					cosverteilung_collerated_minus->Fill(cospi0,beamphoton1E,time);

				}


				if(helicity==1 || helicity==0){cosverteilung_forE_collerated->Fill(cospi0,beamphoton1E,time);}

				if(cbtrigger==kTRUE || GetScalers()->GetNEntries()==0){

					h_energy_sum_pi0_2_3ped->Fill(GetTrigger()->GetEnergySum(),time);
					//cosmeson_beamenergy_energysum_rek->Fill(cospi0,beamphoton1E,GetTrigger()->GetEnergySum(),time);
					//cosmeson_beamenergy_rek_2_3ped->Fill(beamphoton1E,cospi0,time);
				}
						
				openingangle_gammatogamma_energy_allcuts->Fill(opening_angle_gammatogamma,beamphoton1E,time);


				if(pt>0 ){

					cosverteilung_collerated_pt->Fill(cospi0,beamphoton1E,time,pt);

					if(helicity==1){
						cos_beam_hel1_targetplus->Fill(cospi0,beamphoton1E,time);
						cos_beam_hel1_targetplus_pc->Fill(cospi0,beamphoton1E,time,circpol);
						cos_beam_hel1_targetplus_pt->Fill(cospi0,beamphoton1E,time,pt);
					}

					if(helicity==0){
						cos_beam_hel0_targetplus->Fill(cospi0,beamphoton1E,time);
						cos_beam_hel0_targetplus_pc->Fill(cospi0,beamphoton1E,time,circpol);
						cos_beam_hel0_targetplus_pt->Fill(cospi0,beamphoton1E,time,pt);
					}
		
					if(planesetting=="PARA"){
						kristallminus_targetplus_collerated->Fill(phimeson,cospi0,beamphoton1E,time);
						kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,time,pb);
						kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,time,pt);
					}
	
					if(planesetting=="PERP"){
						kristallplus_targetplus_collerated->Fill(phimeson,cospi0,beamphoton1E,time);
						kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,time,pb);
						kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,time,pt);	
					}
				}

				if(pt<0 ){

					cosverteilung_collerated_pt->Fill(cospi0,beamphoton1E,time,-1.*pt);

					if(helicity==1){
						cos_beam_hel1_targetminus->Fill(cospi0,beamphoton1E,time);
						cos_beam_hel1_targetminus_pc->Fill(cospi0,beamphoton1E,time,circpol);
						cos_beam_hel1_targetminus_pt->Fill(cospi0,beamphoton1E,time,pt);
					}

					if(helicity==0){
						cos_beam_hel0_targetminus->Fill(cospi0,beamphoton1E,time);
						cos_beam_hel0_targetminus_pc->Fill(cospi0,beamphoton1E,time,circpol);
						cos_beam_hel0_targetminus_pt->Fill(cospi0,beamphoton1E,time,pt);
					}

					if(planesetting=="PARA"){
						kristallminus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E,time);
						kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,time,pb);
						kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,time,pt);
					}
	
					if(planesetting=="PERP"){
						kristallplus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E,time);
						kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,time,pb);
						kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,time,pt);
				
					}
	
				}

			}
		
		}//GetTagger()
	}//2 GetPhotons()

	anzahlpi0_prompt=anzahlpi0_prompt-HistoManu::backgroundSubstractionFactor*anzahlpi0_sideband;
	if(cutsfullfilled==kTRUE){h_anzahl_pi0->Fill(anzahlpi0_prompt);}
}


void PPi0Analyse::fOnEndProcessing() {

}


void	PPi0Analyse::ProcessScalerRead()
{
    //time.ScalerReadCorrection(5);
}

Bool_t PPi0Analyse::Write()
{

TDirectory* curDir8  = outputFile->mkdir("Observable_E");
curDir8->cd();
cos_beam_hel1_targetminus_proton->Write();
cos_beam_hel1_targetplus_proton->Write();
cos_beam_hel0_targetminus_proton->Write();
cos_beam_hel0_targetplus_proton->Write();

cos_beam_hel1_targetminus_proton_pc->Write();
cos_beam_hel1_targetplus_proton_pc->Write();
cos_beam_hel0_targetminus_proton_pc->Write();
cos_beam_hel0_targetplus_proton_pc->Write();

cos_beam_hel1_targetminus_proton_pt->Write();
cos_beam_hel1_targetplus_proton_pt->Write();
cos_beam_hel0_targetminus_proton_pt->Write();
cos_beam_hel0_targetplus_proton_pt->Write();

cos_beam_hel1_targetminus->Write();
cos_beam_hel1_targetplus->Write();
cos_beam_hel0_targetminus->Write();
cos_beam_hel0_targetplus->Write();

cos_beam_hel1_targetminus_pc->Write();
cos_beam_hel1_targetplus_pc->Write();
cos_beam_hel0_targetminus_pc->Write();
cos_beam_hel0_targetplus_pc->Write();

cos_beam_hel1_targetminus_pt->Write();
cos_beam_hel1_targetplus_pt->Write();
cos_beam_hel0_targetminus_pt->Write();
cos_beam_hel0_targetplus_pt->Write();

cosverteilung_forE_collerated->Write();
cosverteilung_forE_collerated_proton->Write();

missingmassverteilung_forE_collerated_proton->Write();
missingmassverteilung_forE_collerated->Write();

coplanarityverteilung_forE_collerated->Write();
TDirectory* curDir1;
if(runNumber>=1363 && runNumber<=3374){//for carbon
	polsetting=PPi0Analyse::poledge(inputFile);
	TString *data=new TString(polsetting);
	curDir1  = outputFile->mkdir(Form("%s",data->Data()));
}
else{curDir1  = outputFile->mkdir(Form("Butanol%i",(Int_t)GetLinpol()->GetEdgeSetting()));}


TNamed filenumber=TNamed(path11->Data(), "Filenumber");


curDir1->cd();
if(GetScalers()->GetNScalers() == 0){


	TNamed filenumber1=TNamed(Form("%s_0scaler",path11->Data()), "Filenumber1");
	filenumber1.Write();
}
poltable_energy->Write();
poltable_energy_weight->Write();
filenumber.Write();
triggertest->Write();
h_energy_sum_pi0_2_3ped->Write();
h_energy_sum_pi0_3ped->Write();
test->Write();


beamenergy_gen->Write();
cosmeson_beamenergy_energysum_rek->Write();
cosmeson_beamenergy_energysum_monte->Write();
cosmeson_beamenergy_energysum_rek_proton->Write();
cosmeson_beamenergy_monte->Write();
cosmeson_beamenergy_rek_2_3ped->Write();
cosmeson_beamenergy_rek_3ped->Write();


curDir1->cd();
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

Check_TAPS_TOF_photon_2_3ped->Write();

openingangle_gammatogamma_energy_allcuts->Write();


openingangle_gammatogamma_energy_nocuts->Write();

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
coplanarity_theta_CMS->Write();
coplanarity_theta_CMS_all->Write();

h_anzahl_pi0->Write();

Check_CBdE_E_all->Write();
Check_CBdE_E_CBandPID->Write();
Check_CBdE_E_CBandMWPC->Write();
Check_TAPSdE_E->Write();

Check_CBdE_E_nocuts_all->Write();
Check_CBdE_E_nocuts_CBandPID->Write();
Check_CBdE_E_nocuts_CBandMWPC->Write();
Check_TAPSdE_E_nocuts->Write();

Check_TAPS_TOF_photon_3ped->Write();
Check_TAPS_TOF_proton->Write();

openingangle_gammatogamma_energy_allcuts_proton->Write();
openingangle_ptopi0_energy_nocuts->Write();
openingangle_ptopi0_energy_allcuts->Write();

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
invmassverteilung_collerated->Write();
missingmassverteilung_cospi0_collerated->Write();
invmassverteilung_cospi0_collerated->Write();

cosverteilung_collerated->Write();
cosverteilung_collerated_pt->Write();
cosverteilung_collerated_pbminus->Write();
cosverteilung_collerated_pbplus->Write();
cosverteilung_collerated_minus->Write();
cosverteilung_collerated_plus->Write();

events_witherror->Write();
events_all->Write();

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

thetaproton_cospi0_energy_collerated_withoutcuts->Write();
thetaproton_phiproton_collerated->Write();
thetaproton_detector_collerated->Write();

cosverteilung_collerated_proton->Write();
cosverteilung_collerated_proton_pbplus->Write();
cosverteilung_collerated_proton_pbminus->Write();
cosverteilung_collerated_proton_plus->Write();
cosverteilung_collerated_proton_minus->Write();

cosverteilung_collerated_proton_pt->Write();

thetaverteilung_collerated->Write();
coplanarityverteilung_collerated->Write();
missingmassverteilung_collerated_proton->Write();
invmassverteilung_collerated_proton->Write();
coplanarityverteilung_cospi0_collerated->Write();
thetaverteilung_cospi0_collerated->Write();
missingmassverteilung_cospi0_collerated_proton->Write();
invmassverteilung_cospi0_collerated_proton->Write();

clustersize_cospi0_energy_collerated->Write();
clustersize_cospi0_energy_collerated_proton->Write();
thetaproton_cospi0_energy_collerated->Write();
outputFile->cd();
outputFile->Close();
return 0;
}
