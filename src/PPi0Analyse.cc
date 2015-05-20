#include "PPi0Analyse.h"

PPi0Analyse::PPi0Analyse()
{ 

time_prompt = new TH1F("time_prompt", 	"time_prompt", 	1400, -700, 700);
time1= new TH1F("time1", 	"time1", 	1400, -700, 700);
time_side 	= new TH1F("time_side", 	"time_side", 	1400, -700, 700);

 
IM 		= new TH1F("IM", 	"IM", 		1000,   0, 1000);
IM_all          = new TH1F("IM_all",         "IM_all",           1000,   0, 1000);
IM_proton 		= new TH1F("IM_proton", 	"IM_proton", 		1000,   0, 1000);
IM_all_proton          = new TH1F("IM_all_proton",         "IM_all_proton",           1000,   0, 1000);

MM	= new TH1F("MM", 	"MM", 	 	400,   800, 1200);
MM_all          = new TH1F("MM_all",         "MM_all",           400,   800, 1200);	
MM_proton		= new TH1F("MM_proton", 	"MM_proton", 	 	400,   800, 1200);
MM_all_proton          = new TH1F("MM_all_proton",         "MM_all_proton",           400,   800, 1200);

theta_proton 	= new TH1F("theta_proton", 	"theta_proton", 		400,   -200, 200);
theta_all_proton       = new TH1F("theta_all_proton",      "theta_all_proton",                400,   -200, 200);

coplanarity_proton	= new TH1F("coplanarity_proton","coplanarity_proton",400,0,360);
coplanarity_all_proton = new TH1F("coplanarity_all_proton","coplanarity_all_proton",400,0,360);
 
poltable_energy          = new TH1F("poltable_energy",         "poltable_energy",           352,   0, 1448);
poltable_energy_weight          = new TH1F("poltable_energy_weight",         "poltable_energy_weight",           352,   0, 1448);


Check_CBdE_E_nocuts		= new TH2F("Check_CBdE_E_nocuts", "dE_E (all CB clusters compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID} [MeV]", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E_nocuts		= new TH2F("Check_TAPSdE_E_nocuts", "dE_E (all TAPS clusters compared to Veto hits);E_{dep}^{TAPS} [MeV];#Delta E_{dep}^{VETO} [MeV]", 	400, 0, 400, 100, 0, 10);
Check_CBdE_E		= new TH2F("Check_CBdE_E", "dE_E (all CB clusters compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID} [MeV]", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E		= new TH2F("Check_TAPSdE_E", "dE_E (all TAPS clusters compared to Veto hits);E_{dep}^{TAPS} [MeV];#Delta E_{dep}^{VETO} [MeV]", 	400, 0, 400, 100, 0, 10);
Check_TAPS_TOF_proton		= new TH2F("Check_TAPS_TOF_proton", "TOF analysis (proton);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);
Check_TAPS_TOF_photon_3ped		= new TH2F("Check_TAPS_TOF_photon_3ped", "TOF analysis (photon_3ped);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);
Check_TAPS_TOF_photon_2_3ped		= new TH2F("Check_TAPS_TOF_photon_2_3ped", "TOF analysis (photon_2_3ped);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);

openingangle_ptopi0_energy_allcuts = new TH2F("openingangle_ptopi0_energy_allcuts","openingangle_ptopi0_energy_allcuts;#angle #pi,p;E_{beam} [MeV]",180,0,180,410,200,1430);
openingangle_gammatogamma_energy_allcuts = new TH2F("openingangle_gammatogamma_energy_allcuts","openingangle_gammatogamma_energy_allcuts;#angle #gamma,#gamma;E_{beam} [MeV]",180,0,180,410,200,1430);

openingangle_gammatogamma_energy_allcuts_proton = new TH2F("openingangle_gammatogamma_energy_allcuts_proton","openingangle_gammatogamma_energy_allcuts;#angle #gamma,#gamma;E_{beam} [MeV]",180,0,180,410,200,1430);

openingangle_ptopi0_energy_nocuts = new TH2F("openingangle_ptopi0_energy_nocuts","openingangle_ptopi0_energy_nocuts;#angle #pi,p;E_{beam} [MeV]",180,0,180,410,200,1430);
openingangle_gammatogamma_energy_nocuts = new TH2F("openingangle_gammatogamma_energy_nocuts","openingangle_gammatogamma_energy_nocuts;#angle #gamma,#gamma;E_{beam} [MeV]",180,0,180,410,200,1430);


test = new TH1F("test","test",1100,-5.5,5.5);

//kinematic variables in dependence of energy
IM_energy_kplustplus = new TH3F("IM_energy_kplustplus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
IM_energy_kplustminus = new TH3F("IM_energy_kplustminus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
IM_energy_kminustplus = new TH3F("IM_energy_kminustplus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
IM_energy_kminustminus = new TH3F("IM_energy_kminustminus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);

MM_energy_kplustplus = new TH3F("MM_energy_kplustplus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);
MM_energy_kplustminus = new TH3F("MM_energy_kplustminus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);
MM_energy_kminustplus = new TH3F("MM_energy_kminustplus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);
MM_energy_kminustminus = new TH3F("MM_energy_kminustminus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);

// //kinematic variables in dependence of energy
IM_energy_kplustplus_proton = new TH3F("IM_energy_kplustplus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
IM_energy_kplustminus_proton = new TH3F("IM_energy_kplustminus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
IM_energy_kminustplus_proton = new TH3F("IM_energy_kminustplus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
IM_energy_kminustminus_proton = new TH3F("IM_energy_kminustminus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);

MM_energy_kplustplus_proton = new TH3F("MM_energy_kplustplus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);
MM_energy_kplustminus_proton = new TH3F("MM_energy_kplustminus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);
MM_energy_kminustplus_proton = new TH3F("MM_energy_kminustplus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);
MM_energy_kminustminus_proton = new TH3F("MM_energy_kminustminus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,410,200,1430);

coplanarityverteilung_collerated = new TH2F("coplanarityverteilung_collerated", "Coplanarity-Verteilung;#phi_{#pi}-#phi_{p}[deg]", 400,0,360,410,200,1430);
coplanarityverteilung_collerated->Sumw2();

thetaverteilung_collerated = new TH2F("thetaverteilung_collerated", "Polar-Verteilung;#Delt #theta_{p}[deg]", 400,-200, 200,410,200,1430);
thetaverteilung_collerated->Sumw2();

missingmassverteilung_collerated = new TH2F("missingmassverteilung_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,410,200,1430);
missingmassverteilung_collerated->Sumw2();

missingmassverteilung_collerated_proton = new TH2F("missingmassverteilung_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,410,200,1430);
missingmassverteilung_collerated_proton->Sumw2();

invmassverteilung_collerated = new TH2F("invmassverteilung_collerated", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,410,200,1430);
invmassverteilung_collerated->Sumw2();

invmassverteilung_collerated_proton = new TH2F("invmassverteilung_collerated_proton", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,410,200,1430);
invmassverteilung_collerated_proton->Sumw2();

//3-dimensional cuts
coplanarityverteilung_cospi0_collerated = new TH3F("coplanarityverteilung_cospi0_collerated", "Coplanarity-Verteilung;#phi_{#pi}-#phi_{p}[deg]",400,0,360,72,-1,1,410,200,1430);
coplanarityverteilung_cospi0_collerated->Sumw2();

clustersize_cospi0_energy_collerated = new TH3F("clustersize_cospi0_energy_collerated", "clustersize_cospi0_energy_collerated;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,410,200,1430);
clustersize_cospi0_energy_collerated->Sumw2();

clustersize_cospi0_energy_collerated_proton = new TH3F("clustersize_cospi0_energy_collerated_proton", "clustersize_cospi0_energy_collerated;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,410,200,1430);
clustersize_cospi0_energy_collerated_proton->Sumw2();

thetaproton_cospi0_energy_collerated = new TH3F("thetaproton_cospi0_energy_collerated", "#theta_{p} [deg]",360,0,180,72,-1,1,410,200,1430);
thetaproton_cospi0_energy_collerated->Sumw2();

thetaverteilung_cospi0_collerated = new TH3F("thetaverteilung_cospi0_collerated", "Polar-Verteilung;#Delt #theta_{p}[deg]", 400,-200, 200,72,-1,1,410,200,1430);
thetaverteilung_cospi0_collerated->Sumw2();

missingmassverteilung_cospi0_collerated = new TH3F("missingmassverteilung_cospi0_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,72,-1,1,410,200,1430);
missingmassverteilung_cospi0_collerated->Sumw2();

missingmassverteilung_cospi0_collerated_proton = new TH3F("missingmassverteilung_cospi0_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,72,-1,1,410,200,1430);
missingmassverteilung_cospi0_collerated_proton->Sumw2();

invmassverteilung_cospi0_collerated = new TH3F("invmassverteilung_cospi0_collerated", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
invmassverteilung_cospi0_collerated->Sumw2();

invmassverteilung_cospi0_collerated_proton = new TH3F("invmassverteilung_cospi0_collerated_proton", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,72,-1,1,410,200,1430);
invmassverteilung_cospi0_collerated_proton->Sumw2();

//kinematic variables in dependence of energy
cosverteilung_collerated = new TH2F("cosverteilung_collerated", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cosverteilung_collerated->Sumw2();

// //kinematic variables in dependence of energy
cosverteilung_collerated_proton = new TH2F("cosverteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cosverteilung_collerated_proton->Sumw2();
// 

events_all = new TH3F("events_all","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
events_witherror = new TH3F("events_witherror","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);

//helicity analysis for proton
cos_beam_hel1_targetplus_proton = new TH2F("cos_beam_hel1_targetplus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetplus_proton->Sumw2();

cos_beam_hel0_targetplus_proton = new TH2F("cos_beam_hel0_targetplus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetplus_proton->Sumw2();

cos_beam_hel1_targetminus_proton = new TH2F("cos_beam_hel1_targetminus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetminus_proton->Sumw2();

cos_beam_hel0_targetminus_proton = new TH2F("cos_beam_hel0_targetminus_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetminus_proton->Sumw2();

cos_beam_hel1_targetplus_proton_pc = new TH2F("cos_beam_hel1_targetplus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetplus_proton_pc->Sumw2();

cos_beam_hel0_targetplus_proton_pc = new TH2F("cos_beam_hel0_targetplus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetplus_proton_pc->Sumw2();

cos_beam_hel1_targetminus_proton_pc = new TH2F("cos_beam_hel1_targetminus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetminus_proton_pc->Sumw2();

cos_beam_hel0_targetminus_proton_pc = new TH2F("cos_beam_hel0_targetminus_proton_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetminus_proton_pc->Sumw2();

cos_beam_hel1_targetplus_proton_pt = new TH2F("cos_beam_hel1_targetplus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetplus_proton_pt->Sumw2();

cos_beam_hel0_targetplus_proton_pt = new TH2F("cos_beam_hel0_targetplus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetplus_proton_pt->Sumw2();

cos_beam_hel1_targetminus_proton_pt = new TH2F("cos_beam_hel1_targetminus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetminus_proton_pt->Sumw2();

cos_beam_hel0_targetminus_proton_pt = new TH2F("cos_beam_hel0_targetminus_proton_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetminus_proton_pt->Sumw2();

//sigma and g analysis for proton
kristallminus_targetplus_collerated_proton = new TH3F("kristallminus_targetplus_collerated_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetplus_collerated_proton->Sumw2();

kristallminus_targetminus_collerated_proton = new TH3F("kristallminus_targetminus_collerated_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetminus_collerated_proton->Sumw2();

kristallplus_targetplus_collerated_proton = new TH3F("kristallplus_targetplus_collerated_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetplus_collerated_proton->Sumw2();

kristallplus_targetminus_collerated_proton = new TH3F("kristallplus_targetminus_collerated_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetminus_collerated_proton->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb_proton = new TH3F("kristallminus_targetplus_collerated_pb_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetplus_collerated_pb_proton->Sumw2();

kristallminus_targetminus_collerated_pb_proton = new TH3F("kristallminus_targetminus_collerated_pb_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetminus_collerated_pb_proton->Sumw2();

kristallplus_targetplus_collerated_pb_proton = new TH3F("kristallplus_targetplus_collerated_pb_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetplus_collerated_pb_proton->Sumw2();

kristallplus_targetminus_collerated_pb_proton = new TH3F("kristallplus_targetminus_collerated_pb_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetminus_collerated_pb_proton->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt_proton = new TH3F("kristallminus_targetplus_collerated_pt_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetplus_collerated_pt_proton->Sumw2();

kristallminus_targetminus_collerated_pt_proton = new TH3F("kristallminus_targetminus_collerated_pt_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetminus_collerated_pt_proton->Sumw2();

kristallplus_targetplus_collerated_pt_proton = new TH3F("kristallplus_targetplus_collerated_pt_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetplus_collerated_pt_proton->Sumw2();

kristallplus_targetminus_collerated_pt_proton = new TH3F("kristallplus_targetminus_collerated_pt_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetminus_collerated_pt_proton->Sumw2();

//helicity analysis without proton 

cos_beam_hel1_targetplus = new TH2F("cos_beam_hel1_targetplus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetplus->Sumw2();

cos_beam_hel0_targetplus = new TH2F("cos_beam_hel0_targetplus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetplus->Sumw2();

cos_beam_hel1_targetminus = new TH2F("cos_beam_hel1_targetminus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetminus->Sumw2();

cos_beam_hel0_targetminus = new TH2F("cos_beam_hel0_targetminus", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetminus->Sumw2();

cos_beam_hel1_targetplus_pc = new TH2F("cos_beam_hel1_targetplus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetplus_pc->Sumw2();

cos_beam_hel0_targetplus_pc = new TH2F("cos_beam_hel0_targetplus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetplus_pc->Sumw2();

cos_beam_hel1_targetminus_pc = new TH2F("cos_beam_hel1_targetminus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetminus_pc->Sumw2();

cos_beam_hel0_targetminus_pc = new TH2F("cos_beam_hel0_targetminus_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetminus_pc->Sumw2();

cos_beam_hel1_targetplus_pt = new TH2F("cos_beam_hel1_targetplus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetplus_pt->Sumw2();

cos_beam_hel0_targetplus_pt = new TH2F("cos_beam_hel0_targetplus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetplus_pt->Sumw2();

cos_beam_hel1_targetminus_pt = new TH2F("cos_beam_hel1_targetminus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel1_targetminus_pt->Sumw2();

cos_beam_hel0_targetminus_pt = new TH2F("cos_beam_hel0_targetminus_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,410,200,1430);
cos_beam_hel0_targetminus_pt->Sumw2();

//analysis for sigma and g without proton
kristallminus_targetplus_collerated = new TH3F("kristallminus_targetplus_collerated","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetplus_collerated->Sumw2();

kristallminus_targetminus_collerated = new TH3F("kristallminus_targetminus_collerated","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetminus_collerated->Sumw2();

kristallplus_targetplus_collerated = new TH3F("kristallplus_targetplus_collerated","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetplus_collerated->Sumw2();

kristallplus_targetminus_collerated = new TH3F("kristallplus_targetminus_collerated","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetminus_collerated->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb = new TH3F("kristallminus_targetplus_collerated_pb","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetplus_collerated_pb->Sumw2();

kristallminus_targetminus_collerated_pb = new TH3F("kristallminus_targetminus_collerated_pb","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetminus_collerated_pb->Sumw2();

kristallplus_targetplus_collerated_pb = new TH3F("kristallplus_targetplus_collerated_pb","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetplus_collerated_pb->Sumw2();

kristallplus_targetminus_collerated_pb = new TH3F("kristallplus_targetminus_collerated_pb","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetminus_collerated_pb->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt = new TH3F("kristallminus_targetplus_collerated_pt","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetplus_collerated_pt->Sumw2();

kristallminus_targetminus_collerated_pt = new TH3F("kristallminus_targetminus_collerated_pt","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallminus_targetminus_collerated_pt->Sumw2();

kristallplus_targetplus_collerated_pt = new TH3F("kristallplus_targetplus_collerated_pt","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetplus_collerated_pt->Sumw2();

kristallplus_targetminus_collerated_pt = new TH3F("kristallplus_targetminus_collerated_pt","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,410,200,1430);
kristallplus_targetminus_collerated_pt->Sumw2();


beamenergy_gen  = new TH1F("beamenergy_gen ", "beamenergy_gen ",  300, 0, 1557);
beamenergy_gen->Sumw2();

cosmeson_beamenergy_energysum_rek = new TH3F("cosmeson_beamenergy_energysum_rek", "cosmeson_beamenergy_energysum_rek;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1, 30, 220, 1420,200, 0, 2000);
cosmeson_beamenergy_energysum_rek->Sumw2();

cosmeson_beamenergy_energysum_rek_proton = new TH3F("cosmeson_beamenergy_energysum_rek_proton", "cosmeson_beamenergy_energysum_rek_proton;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",18, -1, 1,  30, 220, 1420,200, 0, 2000);
cosmeson_beamenergy_energysum_rek_proton->Sumw2();

cosmeson_beamenergy_energysum_monte = new TH3F("cosmeson_beamenergy_energysum_monte", "cosmeson_beamenergy_energysum_monte;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);
cosmeson_beamenergy_energysum_monte->Sumw2();

cosmeson_beamenergy_monte = new TH2F("cosmeson_beamenergy_monte", "cosmeson_beamenergy_monte;E_{#gamma} [MeV];cos #theta_{#pi}",  30, 220, 1420,18, -1, 1);
cosmeson_beamenergy_monte->Sumw2();

cosmeson_beamenergy_rek_3ped = new TH2F("cosmeson_beamenergy_rek_3ped", "cosmeson_beamenergy_rek_3ped;E_{#gamma} [MeV];cos #theta_{#pi}",  30, 220, 1420,18, -1, 1);
cosmeson_beamenergy_rek_3ped->Sumw2();

cosmeson_beamenergy_rek_2_3ped = new TH2F("cosmeson_beamenergy_rek_2_3ped", "cosmeson_beamenergy_rek_2_3ped;E_{#gamma} [MeV];cos #theta_{#pi}",  30, 220, 1420,18, -1, 1);
cosmeson_beamenergy_rek_2_3ped->Sumw2();

h_energy_sum_pi0_3ped = new TH1F("h_energy_sum_pi0_3ped ", "h_energy_sum_pi0_3ped;E_{sum} [MeV]",  300, 0, 1557);
h_energy_sum_pi0_3ped->Sumw2();

h_energy_sum_pi0_2_3ped = new TH1F("h_energy_sum_pi0_2_3ped ", "h_energy_sum_pi0_3ped;E_{sum} [MeV]",  300, 0, 1557);
h_energy_sum_pi0_2_3ped->Sumw2();

triggertest = new TH1F("triggertest","triggerpattern",34,0,34);

pionmasse = 134.9766;
unten_copl = -25.;
oben_copl = 25.;
unten_theta=-10.;
oben_theta=10.;

//DEFINE TCUTG for dE_E!

//    cuttaps->SetVarX("dE_E (all TAPS clusters compared to Veto hits)");
//    cuttaps->SetVarY("");
//    cuttaps->SetTitle("Graph");
//    cuttaps->SetFillColor(1);
//    cuttaps->SetPoint(0,20.1219,8.68497);
//    cuttaps->SetPoint(1,20.1219,5.10838);
//    cuttaps->SetPoint(2,104.81,1.45954);
//    cuttaps->SetPoint(3,400.542,1.35116);
//    cuttaps->SetPoint(4,399.865,3.98844);
//    cuttaps->SetPoint(5,184.756,4.18714);
//    cuttaps->SetPoint(6,19.1057,8.95592);
//    cuttaps->SetPoint(7,20.1219,8.54046);
//    cuttaps->SetPoint(8,20.4607,8.08887);
//    cuttaps->SetPoint(9,20.1219,8.68497);

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


//    cutcb->SetVarX("dE_E (all CB clusters compared to PID hits)");
//    cutcb->SetVarY("");
//    cutcb->SetTitle("Graph");
//    cutcb->SetTitle("ProtonCut");
//    cutcb->SetFillColor(1);
//    cutcb->SetPoint(0,1.84744,4.57541+0.5);
//    cutcb->SetPoint(1,60.8462,2.26562+0.5);
//    cutcb->SetPoint(2,261.085,1.0428+0.5);
//    cutcb->SetPoint(3,396.961,0.872962+0.5);
//    cutcb->SetPoint(4,396.961,1.50136+0.5);
//    cutcb->SetPoint(5,246.186,2.91101+0.5);
//    cutcb->SetPoint(6,131.764,4.76223+0.5);
//    cutcb->SetPoint(7,31.6448,9.12704+0.5);
//    cutcb->SetPoint(8,1.84744,9.12704+0.5);
//    cutcb->SetPoint(9,1.84744,4.57541+0.5);
//    cutcb->Draw("");


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
Int_t runNumber;
ss >> runNumber;
// planesetting=PPi0Analyse::polplane(inputFile);
if(GetLinpol()->GetPolarizationPlane()==1){planesetting="PERP";}
if(GetLinpol()->GetPolarizationPlane()==0){planesetting="PARA";}
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();


	return kTRUE;
}

void PPi0Analyse::fOnBeforeEventProcessing() {

}

Double_t cutPromptMin=-8.;
Double_t cutPromptMax=8.;
Double_t cutSideMin=100.;
Double_t cutSideMax=300.;
Double_t backgroundSubstractionFactor = (cutPromptMax - cutPromptMin)/(2*(cutSideMax - cutSideMin));

bool IsPromptt(Double_t value)
{
   if ((value >= cutPromptMin) && (value <= cutPromptMax)){
       return 1;}else{
	   return 0;}
	
}
	
bool IsRandomm(Double_t value)
{

	if((value > (-cutSideMax) && value < (-cutSideMin)) ||(value > cutSideMin && value < cutSideMax)){	
            return 1;}else{
		    return 0;}
}

void FillTH1_timeweighted(TH1F *h1, Double_t fillvalue,Double_t timing)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue);
	}


}

void FillTH1_timeandvalueweighted(TH1F *h1, Double_t fillvalue,Double_t timing, Double_t weight)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,weight*backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,weight);
	}


}

void FillTH2_timeweighted(TH2F *h1, Double_t fillvalue,Double_t fillvalue2,Double_t timing)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2);
	}


}

void FillTH2_timeandvalueweighted(TH2F *h1, Double_t fillvalue,Double_t fillvalue2, Double_t timing, Double_t weight)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,weight*backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2,weight);
	}


}

void FillTH3_timeweighted(TH3F *h1, Double_t fillvalue,Double_t fillvalue2, Double_t fillvalue3,Double_t timing)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3,backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3);
	}


}

void FillTH3_timeandvalueweighted(TH3F *h1, Double_t fillvalue,Double_t fillvalue2,Double_t fillvalue3,Double_t timing, Double_t weight)
{
	if(IsRandomm(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3,weight*backgroundSubstractionFactor);
	}

	if(IsPromptt(timing)){
		h1->Fill(fillvalue,fillvalue2,fillvalue3,weight);
	}


}



void	PPi0Analyse::ProcessEvent()	
{

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

/*	if(GetScalers()->GetNEntries()==0){//MC has no GetScalers()
		TLorentzVector beam_mc = GetGeant()->GetBeam();
		Double_t e_beam_mc = beam_mc.E()*1000.;
		//    cout <<  beam_mc.E()*1000<< endl;        
		beamenergy_gen->Fill(e_beam_mc);
		
		//    cout << " Number of generated particles: "<<GetGeant()->GetNTrueParticles() << endl;
		
		TLorentzVector proton_mc = GetGeant()->GetTrueVector(0);
		proton_mc.SetPxPyPzE(proton_mc.Px()*1000,proton_mc.Py()*1000,proton_mc.Pz()*1000,proton_mc.E()*1000);
// 		cout << proton_mc.E() << endl;
		TLorentzVector g1_mc = GetGeant()->GetTrueVector(1);
		TLorentzVector g2_mc = GetGeant()->GetTrueVector(2);
		TLorentzVector meson_mc = g1_mc+g2_mc;
		meson_mc.SetPxPyPzE(meson_mc.Px()*1000,meson_mc.Py()*1000,meson_mc.Pz()*1000,meson_mc.E()*1000);

		//Boost-System
		TLorentzVector meson_mc_cms = CMVector(meson_mc, beam_mc, proton_mc);
		Double_t cosT_mc_cms= meson_mc_cms.CosTheta();
		
		//    cout << " invm_12: "<< (proton_mc+g1_mc).M() << " invm23: "<< (g1_mc+g2_mc).M() << endl;
		if(e_beam_mc>=200&& e_beam_mc<=800){
			cosmeson_beamenergy_energysum_monte->Fill(cosT_mc_cms,e_beam_mc,1000*GetGeant()->GetCBESum());
			cosmeson_beamenergy_monte->Fill(e_beam_mc,cosT_mc_cms);
		}
	}*/

	if(GetPhotons()->GetNParticles()==2){
		TLorentzVector pi0_4vect = GetPhotons()->Particle(0)+GetPhotons()->Particle(1);
		//invMass
		Double_t inv=pi0_4vect.M();
	
		//TEST the target polarization (no entry with pt=0)
		test->Fill(pt);

		for(Int_t j=0; j < GetTagger()->GetNTagged();j++)
		{

			if(GetPhotons()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;
			if(GetPhotons()->GetDetectors(1)==GTreeTrack::DETECTOR_PbWO4) continue;
			if(GetRootinos()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;

	 		if(pt==5){continue;}
// 	  		if(pt!=5){continue;}//only carbon runs

			Double_t time= GetTagger()->GetTaggedTime(j) - 0.5*(GetPhotons()->GetTime(1)+GetPhotons()->GetTime(0));
			time1->Fill(time);

			if(time > -8 && time < 8){
				time_prompt->Fill(time);
			}

				if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
					time_side->Fill(time);
				}
			if(GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j))==-1){continue;}//skip events with wrong polarisation!

			if(!(IsRandomm(time) || IsPromptt(time))){continue;} //skip events which are not prompt or random

			poltable_energy->Fill(GetTagger()->GetTaggedEnergy(j));
			poltable_energy_weight->Fill(GetTagger()->GetTaggedEnergy(j),GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j)));

			//get target, beam and missng particle information
			TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
			TLorentzVector  beam_4vect = TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(j),GetTagger()->GetTaggedEnergy(j));
			TLorentzVector  missingp_4vect = beam_4vect + protonvektor_target - pi0_4vect;
			Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
			Double_t beamphoton1E=GetTagger()->GetTaggedEnergy(j);
			Double_t pb=GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j));
			Double_t circpol = GetCircPol(beamphoton1E,runNumber);
			Bool_t helicity =GetTrigger()->GetHelicity();

			//due to kinematic not possible -> kick them out
			if(anglethetaproton_rek > 90){continue;}
			
		//	if(beamphoton1E<200 || beamphoton1E>800){continue;}

			//energy dependent cuts

		/*	Double_t oben_copl=286.412+(-0.481279)*beamphoton1E+(0.000852105)*beamphoton1E*beamphoton1E+(-5.07571e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_copl=76.5693+(0.440704)*beamphoton1E+(-0.00076097)*beamphoton1E*beamphoton1E+(4.4698e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_theta=-53.3809+(0.359972)*beamphoton1E+(-0.000655192)*beamphoton1E*beamphoton1E+(3.93683e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_theta=41.081+(-0.324925)*beamphoton1E+(0.000653989)*beamphoton1E*beamphoton1E+(-4.22857e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_mass=905.506+(0.394287)*beamphoton1E+(-0.000674192)*beamphoton1E*beamphoton1E+(4.45793e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_mass=1004.28+(-0.57981)*beamphoton1E+(0.000946365)*beamphoton1E*beamphoton1E+(-6.617e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_mass_proton=940.549+(0.222691)*beamphoton1E+(-0.000363125)*beamphoton1E*beamphoton1E+(2.34718e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_mass_proton=911.29+(0.0271286)*beamphoton1E+(-0.000312506)*beamphoton1E*beamphoton1E+(1.84e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_inv=166.049+(-0.0587534)*beamphoton1E+(0.000150841)*beamphoton1E*beamphoton1E+(-9.2621e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_inv=102.426+(0.0656183)*beamphoton1E+(-0.000141924)*beamphoton1E*beamphoton1E+(8.53659e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_inv_proton=150.767+(0.0311453)*beamphoton1E+(-2.27935e-05)*beamphoton1E*beamphoton1E+(1.10307e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_inv_proton=111.269+(0.0168011)*beamphoton1E+(-4.45002e-05)*beamphoton1E*beamphoton1E+(2.29099e-08)*beamphoton1E*beamphoton1E*beamphoton1E;*/

			Double_t oben_copl=190+15;
			Double_t unten_copl=190-15;
			Double_t oben_theta=10;
			Double_t unten_theta=-10;
			Double_t oben_mass=938+67;
			Double_t unten_mass=938-67;
			Double_t oben_inv=135+20;
			Double_t unten_inv=135-20;
			Double_t oben_mass_proton=938+67;
			Double_t unten_mass_proton=938-67;
			Double_t oben_inv_proton=135+20;
			Double_t unten_inv_proton=135-20;


			//Missing Mass
			Double_t missingmass=missingp_4vect.M();

			//Boost-System
			TLorentzVector pi0_4vect_boost = CMVector(pi0_4vect, beam_4vect, protonvektor_target);
			TVector3 pi0_3vektor_boost = pi0_4vect_boost.Vect();
			Double_t cospi0 = pi0_3vektor_boost.CosTheta();
			Double_t phimeson = TMath::RadToDeg()*(pi0_4vect.Vect().Phi());

//_________________________________only one charged particle__________________________________________________________________________

			Double_t anglephimeson;
			Double_t anglephiproton;
			Double_t phi;
			Double_t anglethetaproton_meas;
			Double_t thetadiff;

			Double_t opening_angle_gammatogamma=TMath::ACos(TMath::Cos(GetPhotons()->Particle(0).Phi()-(GetPhotons()->Particle(1)).Phi())*TMath::Sin(GetPhotons()->Particle(0).Theta())*TMath::Sin((GetPhotons()->Particle(1)).Theta())+ TMath::Cos(GetPhotons()->Particle(0).Theta())*TMath::Cos((GetPhotons()->Particle(1)).Theta()) )*TMath::RadToDeg();

			if(GetRootinos()->GetNParticles()==1){

				//get charged particle information
				TLorentzVector proton_4vect_meas;
				if(GetRootinos()->GetNParticles()==1){proton_4vect_meas=GetRootinos()->Particle(0);}

				//apply TCUTGs
				if(GetRootinos()->HasCB(0) && cutcb->IsInside(proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0))==1){kCB=10;}else{kCB=10;}
				if(GetRootinos()->HasTAPS(0) && cuttaps->IsInside(proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0))==1){kTAPS=10;}else{kTAPS=10;}

				Double_t opening_angle_ptopi0= TMath::ACos(TMath::Cos(proton_4vect_meas.Phi()-(pi0_4vect).Phi())*TMath::Sin(proton_4vect_meas.Theta())*TMath::Sin((pi0_4vect).Theta())+ TMath::Cos(proton_4vect_meas.Theta())*TMath::Cos((pi0_4vect).Theta()) )*TMath::RadToDeg();

				//Phi-Difference
				anglephimeson = pi0_4vect.Vect().Phi();
				anglephiproton = proton_4vect_meas.Vect().Phi();
				phi = 360*(anglephimeson - anglephiproton)/(2*TMath::Pi());
				if(phi < 0){phi = phi + 360;}
			
				//theta differenz
				anglethetaproton_meas = TMath::RadToDeg()*proton_4vect_meas.Vect().Theta();
				thetadiff=anglethetaproton_rek-anglethetaproton_meas;

				//due to kinematic not possible -> kick them out
				if(anglethetaproton_meas > 90){continue;}

				//without any cuts!
				FillTH2_timeweighted(openingangle_ptopi0_energy_nocuts,opening_angle_ptopi0,beamphoton1E,time);
				FillTH2_timeweighted(openingangle_gammatogamma_energy_nocuts,opening_angle_gammatogamma,beamphoton1E,time);

				if(GetRootinos()->HasCB(0)){
					FillTH2_timeweighted(Check_CBdE_E_nocuts,proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
				}
				if(GetRootinos()->HasTAPS(0)){
					FillTH2_timeweighted(Check_TAPSdE_E_nocuts,proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
				}

				FillTH1_timeweighted(MM_proton,missingmass,time);
				FillTH1_timeweighted(theta_proton,thetadiff,time);
				FillTH1_timeweighted(coplanarity_proton,phi,time);
				FillTH1_timeweighted(IM_proton,inv,time);

				//visualisation of several cuts!

				if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > 	(unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

					FillTH1_timeweighted(theta_all_proton,thetadiff,time);
					FillTH2_timeweighted(thetaverteilung_collerated,thetadiff,beamphoton1E,time);
					FillTH3_timeweighted(thetaverteilung_cospi0_collerated,thetadiff,cospi0,beamphoton1E,time);
				}
	
				if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (thetadiff > (unten_theta) && thetadiff < (oben_theta)) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

					FillTH1_timeweighted(coplanarity_all_proton,phi,time);
					FillTH2_timeweighted(coplanarityverteilung_collerated,phi,beamphoton1E,time);
					FillTH3_timeweighted(coplanarityverteilung_cospi0_collerated,phi,cospi0,beamphoton1E,time);
				}

				if((thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

					if(planesetting=="PARA"){
						if(pt>0 ){
							FillTH3_timeweighted(MM_energy_kminustplus_proton,missingmass,cospi0,beamphoton1E,time);
						}
	
						if(pt<0 ){
							FillTH3_timeweighted(MM_energy_kminustminus_proton,missingmass,cospi0,beamphoton1E,time);
						}
					}	

					if(planesetting=="PERP"){
						if(pt>0 ){
							FillTH3_timeweighted(MM_energy_kplustplus_proton,missingmass,cospi0,beamphoton1E,time);
						}
					
						if(pt<0 ){
							FillTH3_timeweighted(MM_energy_kplustminus_proton,missingmass,cospi0,beamphoton1E,time);
						}
					}	

					FillTH1_timeweighted(MM_all_proton,missingmass,time);
					FillTH2_timeweighted(missingmassverteilung_collerated_proton,missingmass,beamphoton1E,time);
					FillTH3_timeweighted(missingmassverteilung_cospi0_collerated_proton,missingmass,cospi0,beamphoton1E,time);
					}
	
				if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) &&  (phi > ((unten_copl)) && phi < ((oben_copl)))&&(kTAPS==10 || kCB==10)){

					if(planesetting=="PARA"){
						if(pt>0 ){
							FillTH3_timeweighted(IM_energy_kminustplus_proton,inv,cospi0,beamphoton1E,time);
						}
						if(pt<0 ){
							FillTH3_timeweighted(IM_energy_kminustminus_proton,inv,cospi0,beamphoton1E,time);
						}
					}

					if(planesetting=="PERP"){
						if(pt>0 ){
							FillTH3_timeweighted(IM_energy_kplustplus_proton,inv,cospi0,beamphoton1E,time);
						}
						if(pt<0 ){
							FillTH3_timeweighted(IM_energy_kplustminus_proton,inv,cospi0,beamphoton1E,time);
						}
					}
		
					FillTH1_timeweighted(IM_all_proton,inv,time);
					FillTH2_timeweighted(invmassverteilung_collerated_proton,inv,beamphoton1E,time);
					FillTH3_timeweighted(invmassverteilung_cospi0_collerated_proton,inv,cospi0,beamphoton1E,time);			
				}


				//stuff after all cuts and creation of needed analysis histograms_______________________________
	
				if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

					if(GetRootinos()->HasCB(0)){
						FillTH2_timeweighted(Check_CBdE_E,proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
					}

					if(GetRootinos()->HasTAPS(0)){
						FillTH2_timeweighted(Check_TAPSdE_E,proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
						FillTH2_timeweighted(Check_TAPS_TOF_proton,((GetRootinos()->GetTime(0)-GetTagger()->GetTaggedTime(j))/proton_4vect_meas.Vect().Mag())+1/0.299792458,proton_4vect_meas.E(),time);
					}

					if(GetPhotons()->HasTAPS(0)){
						FillTH2_timeweighted(Check_TAPS_TOF_photon_3ped,((GetPhotons()->GetTime(0)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(0).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(0).E(),time);
					}

					if(GetPhotons()->HasTAPS(1)){
						FillTH2_timeweighted(Check_TAPS_TOF_photon_3ped,((GetPhotons()->GetTime(1)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(1).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(1).E(),time);
					}

					FillTH3_timeweighted(thetaproton_cospi0_energy_collerated,anglethetaproton_meas,cospi0,beamphoton1E,time);
					FillTH3_timeweighted(clustersize_cospi0_energy_collerated,GetPhotons()->GetClusterSize(0),cospi0,beamphoton1E,time);
					FillTH3_timeweighted(clustersize_cospi0_energy_collerated,GetPhotons()->GetClusterSize(1),cospi0,beamphoton1E,time);
					FillTH3_timeweighted(clustersize_cospi0_energy_collerated_proton,GetRootinos()->GetClusterSize(0),cospi0,beamphoton1E,time);
		
					FillTH2_timeweighted(cosverteilung_collerated_proton,cospi0,beamphoton1E,time);

					if(cbtrigger==kTRUE || GetScalers()->GetNEntries()==0){

						FillTH1_timeweighted(h_energy_sum_pi0_3ped,GetTrigger()->GetEnergySum(),time);
						FillTH3_timeweighted(cosmeson_beamenergy_energysum_rek_proton,cospi0,beamphoton1E,GetTrigger()->GetEnergySum(),time);
						FillTH2_timeweighted(cosmeson_beamenergy_rek_3ped,beamphoton1E,cospi0,time);
					}
						
					FillTH2_timeweighted(openingangle_ptopi0_energy_allcuts,opening_angle_ptopi0,beamphoton1E,time);
					FillTH2_timeweighted(openingangle_gammatogamma_energy_allcuts_proton,opening_angle_gammatogamma,beamphoton1E,time);

					if(pt>0 ){

						if(helicity==1){
							FillTH2_timeweighted(cos_beam_hel1_targetplus_proton,cospi0,beamphoton1E,time);
							FillTH2_timeandvalueweighted(cos_beam_hel1_targetplus_proton_pc,cospi0,beamphoton1E,time,pb);
							FillTH2_timeandvalueweighted(cos_beam_hel1_targetplus_proton_pt,cospi0,beamphoton1E,time,pt);
						}

						if(helicity==0){
							FillTH2_timeweighted(cos_beam_hel0_targetplus_proton,cospi0,beamphoton1E,time);
							FillTH2_timeandvalueweighted(cos_beam_hel0_targetplus_proton_pc,cospi0,beamphoton1E,time,pb);
							FillTH2_timeandvalueweighted(cos_beam_hel0_targetplus_proton_pt,cospi0,beamphoton1E,time,pt);
						}
		
						if(planesetting=="PARA"){
							FillTH3_timeweighted(kristallminus_targetplus_collerated_proton,phimeson,cospi0,beamphoton1E,time);
							FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pb_proton,phimeson,cospi0,beamphoton1E,time,pb);
							FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pt_proton,phimeson,cospi0,beamphoton1E,time,pt);
						}
	
						if(planesetting=="PERP"){
							FillTH3_timeweighted(kristallplus_targetplus_collerated_proton,phimeson,cospi0,beamphoton1E,time);
							FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pb_proton,phimeson,cospi0,beamphoton1E,time,pb);
							FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pt_proton,phimeson,cospi0,beamphoton1E,time,pt);	
						}
					}

					if(pt<0 ){

						if(helicity==1){
							FillTH2_timeweighted(cos_beam_hel1_targetminus_proton,cospi0,beamphoton1E,time);
							FillTH2_timeandvalueweighted(cos_beam_hel1_targetminus_proton_pc,cospi0,beamphoton1E,time,pb);
							FillTH2_timeandvalueweighted(cos_beam_hel1_targetminus_proton_pt,cospi0,beamphoton1E,time,pt);
						}

						if(helicity==0){
							FillTH2_timeweighted(cos_beam_hel0_targetminus_proton,cospi0,beamphoton1E,time);
							FillTH2_timeandvalueweighted(cos_beam_hel0_targetminus_proton_pc,cospi0,beamphoton1E,time,pb);
							FillTH2_timeandvalueweighted(cos_beam_hel0_targetminus_proton_pt,cospi0,beamphoton1E,time,pt);
						}

						if(planesetting=="PARA"){
							FillTH3_timeweighted(kristallminus_targetminus_collerated_proton,phimeson,cospi0,beamphoton1E,time);
							FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pb_proton,phimeson,cospi0,beamphoton1E,time,pb);
							FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pt_proton,phimeson,cospi0,beamphoton1E,time,pt);
						}
	
						if(planesetting=="PERP"){

							FillTH3_timeweighted(kristallplus_targetminus_collerated_proton,phimeson,cospi0,beamphoton1E,time);
							FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pb_proton,phimeson,cospi0,beamphoton1E,time,pb);
							FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pt_proton,phimeson,cospi0,beamphoton1E,time,pt);
				
						}
	
					}
				}						
// 				
			}
// 				//all events with proton identified



			//without any cuts!

			FillTH1_timeweighted(MM,missingmass,time);
			FillTH1_timeweighted(IM,inv,time);


			if(GetRootinos()->GetNParticles()==1){

					if(!(thetadiff > (unten_theta) && thetadiff < (oben_theta)) ||  !(phi > (unten_copl) && phi < (oben_copl))){continue;} 

			}

			//with cuts
		
			if(((inv > (unten_inv) && inv < (oben_inv)))){

				if(planesetting=="PARA"){

					if(pt>0){
						FillTH3_timeweighted(MM_energy_kminustplus,missingmass,cospi0,beamphoton1E,time);
					}
	
					if(pt<0){
						FillTH3_timeweighted(MM_energy_kminustminus,missingmass,cospi0,beamphoton1E,time);
					}
				}

				if(planesetting=="PERP"){

					if(pt>0 ){
						FillTH3_timeweighted(MM_energy_kplustplus,missingmass,cospi0,beamphoton1E,time);
					}
				
					if(pt<0 ){
						FillTH3_timeweighted(MM_energy_kplustminus,missingmass,cospi0,beamphoton1E,time);
					}
				}

					FillTH1_timeweighted(MM_all,missingmass,time);
					FillTH2_timeweighted(missingmassverteilung_collerated,missingmass,beamphoton1E,time);
					FillTH3_timeweighted(missingmassverteilung_cospi0_collerated,missingmass,cospi0,beamphoton1E,time);
			}
	
			if(((missingmass > (unten_mass) && missingmass < (oben_mass)))){

				if(planesetting=="PARA"){

					if(pt>0){
						FillTH3_timeweighted(IM_energy_kminustplus,inv,cospi0,beamphoton1E,time);
					}
	
					if(pt<0){
						FillTH3_timeweighted(IM_energy_kminustminus,inv,cospi0,beamphoton1E,time);
					}
				}

				if(planesetting=="PERP"){

					if(pt>0 ){
						FillTH3_timeweighted(IM_energy_kplustplus,inv,cospi0,beamphoton1E,time);
					}
				
					if(pt<0 ){
						FillTH3_timeweighted(IM_energy_kplustminus,inv,cospi0,beamphoton1E,time);
					}
				}
				FillTH1_timeweighted(IM_all,inv,time);
				FillTH2_timeweighted(invmassverteilung_collerated,inv,beamphoton1E,time);
				FillTH3_timeweighted(invmassverteilung_cospi0_collerated,inv,cospi0,beamphoton1E,time);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))&&GetTrigger()->GetNErrors()!=0){
				events_witherror->Fill(phimeson,cospi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){
				events_all->Fill(phimeson,cospi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){

	
				if(GetPhotons()->HasTAPS(0)){
					FillTH2_timeweighted(Check_TAPS_TOF_photon_2_3ped,((GetPhotons()->GetTime(0)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(0).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(0).E(),time);
				}

				if(GetPhotons()->HasTAPS(1)){
					FillTH2_timeweighted(Check_TAPS_TOF_photon_2_3ped,((GetPhotons()->GetTime(1)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(1).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(1).E(),time);
				}
		
				FillTH2_timeweighted(cosverteilung_collerated,cospi0,beamphoton1E,time);

				if(cbtrigger==kTRUE || GetScalers()->GetNEntries()==0){

					FillTH1_timeweighted(h_energy_sum_pi0_2_3ped,GetTrigger()->GetEnergySum(),time);
					FillTH3_timeweighted(cosmeson_beamenergy_energysum_rek,cospi0,beamphoton1E,GetTrigger()->GetEnergySum(),time);
					FillTH2_timeweighted(cosmeson_beamenergy_rek_2_3ped,beamphoton1E,cospi0,time);
				}
						
				FillTH2_timeweighted(openingangle_gammatogamma_energy_allcuts,opening_angle_gammatogamma,beamphoton1E,time);


				if(pt>0 ){

					if(helicity==1){
						FillTH2_timeweighted(cos_beam_hel1_targetplus_proton,cospi0,beamphoton1E,time);
						FillTH2_timeandvalueweighted(cos_beam_hel1_targetplus_pc,cospi0,beamphoton1E,time,pb);
						FillTH2_timeandvalueweighted(cos_beam_hel1_targetplus_pt,cospi0,beamphoton1E,time,pt);
					}

					if(helicity==0){
						FillTH2_timeweighted(cos_beam_hel0_targetplus,cospi0,beamphoton1E,time);
						FillTH2_timeandvalueweighted(cos_beam_hel0_targetplus_pc,cospi0,beamphoton1E,time,pb);
						FillTH2_timeandvalueweighted(cos_beam_hel0_targetplus_pt,cospi0,beamphoton1E,time,pt);
					}
		
					if(planesetting=="PARA"){
						FillTH3_timeweighted(kristallminus_targetplus_collerated,phimeson,cospi0,beamphoton1E,time);
						FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pb,phimeson,cospi0,beamphoton1E,time,pb);
						FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pt,phimeson,cospi0,beamphoton1E,time,pt);
					}
	
					if(planesetting=="PERP"){
						FillTH3_timeweighted(kristallplus_targetplus_collerated,phimeson,cospi0,beamphoton1E,time);
						FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pb,phimeson,cospi0,beamphoton1E,time,pb);
						FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pt,phimeson,cospi0,beamphoton1E,time,pt);	
					}
				}

				if(pt<0 ){

					if(helicity==1){
						FillTH2_timeweighted(cos_beam_hel1_targetminus,cospi0,beamphoton1E,time);
						FillTH2_timeandvalueweighted(cos_beam_hel1_targetminus_pc,cospi0,beamphoton1E,time,pb);
						FillTH2_timeandvalueweighted(cos_beam_hel1_targetminus_pt,cospi0,beamphoton1E,time,pt);
					}

					if(helicity==0){
						FillTH2_timeweighted(cos_beam_hel0_targetminus,cospi0,beamphoton1E,time);
						FillTH2_timeandvalueweighted(cos_beam_hel0_targetminus_pc,cospi0,beamphoton1E,time,pb);
						FillTH2_timeandvalueweighted(cos_beam_hel0_targetminus_pt,cospi0,beamphoton1E,time,pt);
					}

					if(planesetting=="PARA"){
						FillTH3_timeweighted(kristallminus_targetminus_collerated,phimeson,cospi0,beamphoton1E,time);
						FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pb,phimeson,cospi0,beamphoton1E,time,pb);
						FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pt,phimeson,cospi0,beamphoton1E,time,pt);
					}
	
					if(planesetting=="PERP"){
						FillTH3_timeweighted(kristallplus_targetminus_collerated,phimeson,cospi0,beamphoton1E,time);
						FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pb,phimeson,cospi0,beamphoton1E,time,pb);
						FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pt,phimeson,cospi0,beamphoton1E,time,pt);
				
					}
	
				}

			}

		}//GetTagger()
	}//2 GetPhotons()
}


void PPi0Analyse::fOnEndProcessing() {

}


void	PPi0Analyse::ProcessScalerRead()
{
    //time.ScalerReadCorrection(5);
}

Bool_t PPi0Analyse::Write()
{

//FOR CARBON
// polsetting=PPi0Analyse::poledge(inputFile);
// TString *data=new TString(polsetting);
// TDirectory* curDir1  = outputFile->mkdir(Form("%s",data->Data()));

TDirectory* curDir1  = outputFile->mkdir(Form("Butanol%i",(Int_t)GetLinpol()->GetEdgeSetting()));
TNamed filenumber=TNamed(path11->Data(), "Filenumber");


curDir1->cd();
poltable_energy->Write();
poltable_energy_weight->Write();
filenumber.Write();
triggertest->Write();
beamenergy_gen->Write();
cosmeson_beamenergy_energysum_rek->Write();
cosmeson_beamenergy_energysum_monte->Write();
cosmeson_beamenergy_energysum_rek_proton->Write();
cosmeson_beamenergy_monte->Write();
cosmeson_beamenergy_rek_2_3ped->Write();
cosmeson_beamenergy_rek_3ped->Write();
h_energy_sum_pi0_2_3ped->Write();
h_energy_sum_pi0_3ped->Write();
test->Write();
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
IM_energy_kplustplus_proton->Write();
IM_energy_kplustminus_proton->Write();
IM_energy_kminustplus_proton->Write();
IM_energy_kminustminus_proton->Write();
Check_CBdE_E->Write();
Check_CBdE_E_nocuts->Write();
Check_TAPSdE_E->Write();
Check_TAPSdE_E_nocuts->Write();
Check_TAPS_TOF_photon_3ped->Write();
Check_TAPS_TOF_proton->Write();

MM_energy_kplustplus_proton->Write();
MM_energy_kplustminus_proton->Write();
MM_energy_kminustplus_proton->Write();
MM_energy_kminustminus_proton->Write();

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

cos_beam_hel1_targetminus->Write();
cos_beam_hel1_targetplus->Write();
cos_beam_hel0_targetminus->Write();
cos_beam_hel0_targetplus->Write();


missingmassverteilung_collerated->Write();
invmassverteilung_collerated->Write();
missingmassverteilung_cospi0_collerated->Write();
invmassverteilung_cospi0_collerated->Write();

cosverteilung_collerated->Write();
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

cos_beam_hel1_targetminus_proton->Write();
cos_beam_hel1_targetplus_proton->Write();
cos_beam_hel0_targetminus_proton->Write();
cos_beam_hel0_targetplus_proton->Write();

cosverteilung_collerated_proton->Write();

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


outputFile->Close();
return 0;
}
