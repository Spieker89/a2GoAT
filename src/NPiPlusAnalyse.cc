#include "NPiPlusAnalyse.h"

NPiPlusAnalyse::NPiPlusAnalyse()
{ 

time_prompt = new TH1F("time_prompt", 	"time_prompt", 	1400, -700, 700);
time1= new TH1F("time1", 	"time1", 	1400, -700, 700);
time_side 	= new TH1F("time_side", 	"time_side", 	1400, -700, 700);


theta_neutron 	= new HistoManu("theta_neutron", 	"theta_neutron", 		400,   -200, 200);
theta_all_neutron       = new HistoManu("theta_all_neutron",      "theta_all_neutron",                400,   -200, 200);
 
poltable_energy          = new HistoManu("poltable_energy",         "poltable_energy",           352,   0, 1448);
poltable_energy_weight          = new HistoManu("poltable_energy_weight",         "poltable_energy_weight",           352,   0, 1448);

h_anzahl_piplus  = new TH1F("h_anzahl_piplus ", "h_anzahl_piplus ",  100, -4, 10);
h_anzahl_piplus->Sumw2();

coplanarity_withoutcuts	= new HistoManu("coplanarity_withoutcuts","coplanarity_withoutcuts",400,0,360);



//PSA analysis.
Check_PSA_piplus= new HistoManu2("Check_PSA_piplus ", "Check_PSA_piplus ",  360, 0, 90, 1200, 0, 600);
Check_PSA_neutron= new HistoManu2("Check_PSA_neutron ", "Check_PSA_neutron ",  360, 0, 90, 1200, 0, 600);
Check_PSA_together= new HistoManu2("Check_PSA_together ", "Check_PSA_together ",  360, 0, 90, 1200, 0, 600);

//TOF spectra
Check_TAPS_TOF_neutron		= new HistoManu2("Check_TAPS_TOF_neutron", "TOF analysis (neutron);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 1200, 0, 600);
Check_TAPS_TOF_piplus		= new HistoManu2("Check_TAPS_TOF_piplus", "TOF analysis (piplus);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 1200, 0, 600);

Check_TAPS_energy_time_neutron		= new HistoManu2("Check_TAPS_time_energy_neutron", "time vs E analysis (neutron);E_{dep} [Mev];t [ns]", 	1200, 0, 600,700,0,35);
Check_TAPS_energy_time_piplus		= new HistoManu2("Check_TAPS_time_energy_piplus", "time vs E (piplus);E_{dep} [Mev];t [ns]",1200, 0, 600,700,0,35);


//dE over E spectra
Check_CBdE_E_nocuts_all		= new HistoManu2("Check_CBdE_E_nocuts_all", "dE_E (all detectors compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID}*sin(#theta) [MeV]", 	400, 0, 400, 100, 0, 10);
Check_CBdE_E_nocuts_CBandPID		= new HistoManu2("Check_CBdE_E_nocuts_CBandPID", "dE_E (CB+PID compared to PID hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{PID} [MeV]", 	400, 0, 400, 100, 0, 10);
Check_CBdE_E_nocuts_CBandMWPC		= new HistoManu2("Check_CBdE_E_nocuts_CBandMWPC", "dE_E (CB+MWPC compared to MWPC hits);E_{dep}^{CB} [MeV];#Delta E_{dep}^{MWPC}*sin(#theta) [MeV]", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E_nocuts		= new HistoManu2("Check_TAPSdE_E_nocuts", "dE_E (all TAPS clusters compared to Veto hits);E_{dep}^{TAPS} [MeV];#Delta E_{dep}^{VETO} [MeV]", 	400, 0, 400, 100, 0, 10);


test = new TH1F("test","test",1100,-5.5,5.5);


clustersize_cospiplus_energy_collerated_neutron_CB = new HistoManu3("clustersize_cospiplus_energy_collerated_neutron_CB", "clustersize_cospiplus_energy_collerated_neutron_CB;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,37,200,1421);

clustersize_cospiplus_energy_collerated_neutron_TAPS = new HistoManu3("clustersize_cospiplus_energy_collerated_neutron_TAPS", "clustersize_cospiplus_energy_collerated_neutron_TAPS;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,37,200,1421);

clustersize_cospiplus_energy_collerated_piplus_CB = new HistoManu3("clustersize_cospiplus_energy_collerated_piplus_CB", "clustersize_cospiplus_energy_collerated_piplus_CB;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,37,200,1421);

clustersize_cospiplus_energy_collerated_piplus_TAPS = new HistoManu3("clustersize_cospiplus_energy_collerated_piplus_TAPS", "clustersize_cospiplus_energy_collerated_piplus_TAPS;clustersize;cos #theta_{#pi};energy",100,0,100,72,-1,1,37,200,1421);

//kinematic variables in dependence of energy
cosverteilung_collerated = new HistoManu2("cosverteilung_collerated", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_forE_collerated = new HistoManu2("cosverteilung_forE_collerated", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

// //kinematic variables in dependence of energy
cosverteilung_collerated_neutron = new HistoManu2("cosverteilung_collerated_neutron", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosverteilung_forE_collerated_neutron = new HistoManu2("cosverteilung_forE_collerated_neutron", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);
// 
events_all = new TH3F("events_all","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);
events_witherror = new TH3F("events_witherror","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//helicity analysis for neutron
cos_beam_hel1_targetplus_neutron = new HistoManu2("cos_beam_hel1_targetplus_neutron", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_neutron = new HistoManu2("cos_beam_hel0_targetplus_neutron", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_neutron = new HistoManu2("cos_beam_hel1_targetminus_neutron", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_neutron = new HistoManu2("cos_beam_hel0_targetminus_neutron", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetplus_neutron_pc = new HistoManu2("cos_beam_hel1_targetplus_neutron_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_neutron_pc = new HistoManu2("cos_beam_hel0_targetplus_neutron_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_neutron_pc = new HistoManu2("cos_beam_hel1_targetminus_neutron_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_neutron_pc = new HistoManu2("cos_beam_hel0_targetminus_neutron_pc", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetplus_neutron_pt = new HistoManu2("cos_beam_hel1_targetplus_neutron_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetplus_neutron_pt = new HistoManu2("cos_beam_hel0_targetplus_neutron_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel1_targetminus_neutron_pt = new HistoManu2("cos_beam_hel1_targetminus_neutron_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

cos_beam_hel0_targetminus_neutron_pt = new HistoManu2("cos_beam_hel0_targetminus_neutron_pt", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,30, 220, 1420);

//sigma and g analysis for neutron
kristallminus_targetplus_collerated_neutron = new HistoManu3("kristallminus_targetplus_collerated_neutron","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kohlenstoffasymmetrie_collerated_neutron = new HistoManu3("kohlenstoffasymmetrie_collerated_neutron","kohlenstoffasymmetrie_collerated_neutron; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_neutron = new HistoManu3("kristallminus_targetminus_collerated_neutron","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_neutron = new HistoManu3("kristallplus_targetplus_collerated_neutron","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_neutron = new HistoManu3("kristallplus_targetminus_collerated_neutron","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb_neutron = new HistoManu3("kristallminus_targetplus_collerated_pb_neutron","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pb_neutron = new HistoManu3("kristallminus_targetminus_collerated_pb_neutron","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pb_neutron = new HistoManu3("kristallplus_targetplus_collerated_pb_neutron","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pb_neutron = new HistoManu3("kristallplus_targetminus_collerated_pb_neutron","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt_neutron = new HistoManu3("kristallminus_targetplus_collerated_pt_neutron","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pt_neutron = new HistoManu3("kristallminus_targetminus_collerated_pt_neutron","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pt_neutron = new HistoManu3("kristallplus_targetplus_collerated_pt_neutron","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pt_neutron = new HistoManu3("kristallplus_targetminus_collerated_pt_neutron","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//helicity analysis without neutron 

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

//analysis for sigma and g without neutron

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
/*
cosmeson_beamenergy_energysum_rek = new HistoManu3("cosmeson_beamenergy_energysum_rek", "cosmeson_beamenergy_energysum_rek;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1, 30, 220, 1420,200, 0, 2000);
cosmeson_beamenergy_energysum_rek->Sumw2();

cosmeson_beamenergy_energysum_rek_neutron = new HistoManu3("cosmeson_beamenergy_energysum_rek_neutron", "cosmeson_beamenergy_energysum_rek_neutron;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",18, -1, 1,  30, 220, 1420,200, 0, 2000);
cosmeson_beamenergy_energysum_rek_neutron->Sumw2();

cosmeson_beamenergy_energysum_monte = new HistoManu3("cosmeson_beamenergy_energysum_monte", "cosmeson_beamenergy_energysum_monte;cos #theta_{#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);
cosmeson_beamenergy_energysum_monte->Sumw2();

cosmeson_beamenergy_monte = new HistoManu2("cosmeson_beamenergy_monte", "cosmeson_beamenergy_monte;E_{#gamma} [MeV];cos #theta_{#pi}",  30, 220, 1420,18, -1, 1);
cosmeson_beamenergy_monte->Sumw2();

cosmeson_beamenergy_rek_3ped = new HistoManu2("cosmeson_beamenergy_rek_3ped", "cosmeson_beamenergy_rek_3ped;E_{#gamma} [MeV];cos #theta_{#pi}",  30, 220, 1420,18, -1, 1);
cosmeson_beamenergy_rek_3ped->Sumw2();

cosmeson_beamenergy_rek_2_3ped = new HistoManu2("cosmeson_beamenergy_rek_2_3ped", "cosmeson_beamenergy_rek_2_3ped;E_{#gamma} [MeV];cos #theta_{#pi}",  30, 220, 1420,18, -1, 1);
cosmeson_beamenergy_rek_2_3ped->Sumw2();*/

h_energy_sum_piplus = new HistoManu("h_energy_sum_piplus ", "h_energy_sum_piplus;E_{sum} [MeV]",  300, 0, 1557);

h_energy_sum_piplus_2_3ped = new HistoManu("h_energy_sum_piplus_2_3ped ", "h_energy_sum_piplus_3ped;E_{sum} [MeV]",  300, 0, 1557);

triggertest = new TH1F("triggertest","triggerpattern",34,0,34);

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

NPiPlusAnalyse::~NPiPlusAnalyse()
{
}



Bool_t	NPiPlusAnalyse::Start()
{
//targetpol
pt=NPiPlusAnalyse::targetpol(inputFile);

//filenumber
TString* filename1 = new TString(inputFile->GetPath());
path11= new TString(filename1->Tokenize("_")->At(filename1->Tokenize("_")->GetEntries()-1)->GetName());
path11->Resize(path11->Length()-7);

//runnumber
stringstream ss(path11->Data());
ss >> runNumber;
if(runNumber>=1363 && runNumber<=3374)planesetting=NPiPlusAnalyse::polplane(inputFile);
else{
if(GetLinpol()->GetPolarizationPlane()==1){planesetting="PERP";}
if(GetLinpol()->GetPolarizationPlane()==0){planesetting="PARA";}
}

    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();


	return kTRUE;
}


void	NPiPlusAnalyse::ProcessEvent()	
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
		
		TLorentzVector neutron_mc = GetGeant()->GetTrueVector(0);
		neutron_mc.SetPxPyPzE(neutron_mc.Px()*1000,neutron_mc.Py()*1000,neutron_mc.Pz()*1000,neutron_mc.E()*1000);
// 		cout << neutron_mc.E() << endl;
		TLorentzVector g1_mc = GetGeant()->GetTrueVector(1);
		TLorentzVector g2_mc = GetGeant()->GetTrueVector(2);
		TLorentzVector meson_mc = g1_mc+g2_mc;
		meson_mc.SetPxPyPzE(meson_mc.Px()*1000,meson_mc.Py()*1000,meson_mc.Pz()*1000,meson_mc.E()*1000);

		//Boost-System
		TLorentzVector meson_mc_cms = CMVector(meson_mc, beam_mc, neutron_mc);
		Double_t cosT_mc_cms= meson_mc_cms.CosTheta();
		
		//    cout << " invm_12: "<< (neutron_mc+g1_mc).M() << " invm23: "<< (g1_mc+g2_mc).M() << endl;
		if(e_beam_mc>=200&& e_beam_mc<=800){
			cosmeson_beamenergy_energysum_monte->Fill(cosT_mc_cms,e_beam_mc,1000*GetGeant()->GetCBESum());
			cosmeson_beamenergy_monte->Fill(e_beam_mc,cosT_mc_cms);
		}
	}*/

	Double_t anzahlpiplus_prompt=0;
	Double_t anzahlpiplus_sideband=0;
	Bool_t cutsfullfilled;
	TLorentzVector neutronvektor_target(0.,0.,0.,938.272046);
	TLorentzVector  beam_4vect;
	TLorentzVector  missingp_4vect;
	TLorentzVector piplus_4vect;

	if(GetPhotons()->GetNParticles()==1 && GetRootinos()->GetNParticles()==1){ //select one neutral and one charged, neutral treated as neutron, charged as piplus
;
	
		//TEST the target polarization (no entry with pt=0)
		test->Fill(pt);

		for(Int_t j=0; j < GetTagger()->GetNTagged();j++)
		{

			if(GetPhotons()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;
			if(GetRootinos()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue; //reject them if they are detected in pbwo4
// 			cout << runNumber << endl;
			if(runNumber>=1363 && runNumber<=3374 && pt!=5 && GetScalers()->GetNEntries()!=0){continue;}//only carbon runs
			if((runNumber<1363 || runNumber>3374) && pt==5 && GetScalers()->GetNEntries()!=0){continue;}//only moeller/diamond runs

			Double_t time= GetTagger()->GetTaggedTime(j) - (GetRootinos()->GetTime(0));//take the piplus as reference time?
			time1->Fill(time);

			if(HistoManu::IsPrompt(time)){
				time_prompt->Fill(time);
			}

			if(HistoManu::IsRandom(time)){	
				time_side->Fill(time);
			}

			Double_t timediff_particles=  GetPhotons()->GetTime(0)- GetRootinos()->GetTime(0); //time diff between the particles


			if(GetScalers()->GetNEntries()!=0 && GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j))==-1){continue;}//skip events with wrong polarisation!

			if(!(HistoManu::IsRandom(time) || HistoManu::IsPrompt(time))){continue;} //skip events which are not prompt or random

			poltable_energy->Fill(GetTagger()->GetTaggedEnergy(j));
			poltable_energy_weight->Fill(GetTagger()->GetTaggedEnergy(j),GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j)));

			//get beam  particle information
			beam_4vect = TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(j),GetTagger()->GetTaggedEnergy(j));
			Double_t beamphoton1E=GetTagger()->GetTaggedEnergy(j);
			Double_t pb=GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j));
			Double_t circpol = GetCircPol(beamphoton1E,runNumber);
			Bool_t helicity =GetTrigger()->GetHelicity();

			
			if(beamphoton1E<200 || beamphoton1E>1430){continue;}


			//Boost-System
			piplus_4vect=GetRootinos()->Particle(0);
			TLorentzVector piplus_4vect_boost = CMVector(piplus_4vect, beam_4vect, neutronvektor_target);
			TVector3 piplus_3vektor_boost = piplus_4vect_boost.Vect();
			Double_t cospiplus = piplus_3vektor_boost.CosTheta();
			Double_t phipiplus = TMath::RadToDeg()*(piplus_4vect.Vect().Phi());
			Double_t phipiplus_CMS = TMath::RadToDeg()*(piplus_4vect_boost.Vect().Phi());
			Double_t thetapiplus_CMS = TMath::RadToDeg()*(piplus_4vect_boost.Vect().Theta());

			//phi cuts
			Double_t phi = phipiplus - GetPhotons()->Particle(0).Vect().Phi()*TMath::RadToDeg();
			if(phi < 0){phi = phi + 360;}

			coplanarity_withoutcuts->Fill(phi,time);


			//visualisation of TCutG
			if(GetRootinos()->GetDetectors(0)==7)Check_CBdE_E_nocuts_all->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0)*TMath::Sin(GetRootinos()->GetTheta(0)),time);

			if(GetRootinos()->GetDetectors(0)==3)Check_CBdE_E_nocuts_CBandPID->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0),time);

			if(GetRootinos()->GetDetectors(0)==5)Check_CBdE_E_nocuts_CBandMWPC->Fill(GetRootinos()->Particle(0).E(),0.5*(GetRootinos()->GetMWPC1Energy(0)+GetRootinos()->GetMWPC0Energy(0))*TMath::Sin(GetRootinos()->GetTheta(0)),time);

			if(GetRootinos()->HasTAPS(0)){
				Check_TAPSdE_E_nocuts->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0),time);
			}

			//PULSE SHAPE ANALYSIS
			for(Int_t i=0; i<GetTracks()->GetNTracks() ;i++){//loop over all tracks

				Double_t clusterenergy = GetTracks()->GetClusterEnergy(i);
					
				Double_t PSA_angle=GetTracks()->GetPSAAngle(i);
				Double_t PSA_radius=GetTracks()->GetPSARadius(i);

				if(GetPhotons()->HasTAPS(0) && clusterenergy== (GetPhotons()->Particle(0)).E() ){
					Check_PSA_together->Fill(PSA_angle,PSA_radius,time);
					Check_PSA_neutron->Fill(PSA_angle,PSA_radius,time);

				}

				if(GetRootinos()->HasTAPS(0) && clusterenergy== (GetRootinos()->Particle(0)).E()){
					Check_PSA_together->Fill(PSA_angle,PSA_radius,time);
					Check_PSA_piplus->Fill(PSA_angle,PSA_radius,time);

				}

			}

			//TAPS TOF 
			if(GetPhotons()->HasTAPS(0)){
				Double_t length_neutron=1.75700/TMath::Cos(GetPhotons()->GetThetaRad(0));
				Double_t TOF_neutron=-GetPhotons()->GetTime(0)+GetTagger()->GetTaggedTime(j);
		
				Check_TAPS_TOF_neutron->Fill((TOF_neutron/length_neutron)+1/0.299792458,GetPhotons()->Particle(0).E(),time);
			}

			if(GetPhotons()->HasTAPS(0)){
				Check_TAPS_energy_time_neutron->Fill(GetPhotons()->Particle(0).E(),GetPhotons()->GetTime(0),time);
			}

			if(GetRootinos()->HasTAPS(0)){
				Double_t length_piplus=1.75700/TMath::Cos(GetRootinos()->GetThetaRad(0));
				Double_t TOF_piplus=-GetRootinos()->GetTime(0)+GetTagger()->GetTaggedTime(j);

				Check_TAPS_TOF_piplus->Fill((TOF_piplus/length_piplus)+1/0.299792458,GetRootinos()->Particle(0).E(),time);

			}

			if(GetRootinos()->HasTAPS(0)){
				Check_TAPS_energy_time_piplus->Fill(GetRootinos()->Particle(0).E(),GetRootinos()->GetTime(0),time);
			}

			//apply TCUTGs

			//all possible particles
			if(GetRootinos()->GetDetectors(0)==7 && cutcb->IsInside(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0)*TMath::Sin(GetRootinos()->GetTheta(0))) ==1){kCB_all=1;}else{kCB_all=0;}

			//only PID AND CB 
			if(GetRootinos()->GetDetectors(0)==3 && cutcb->IsInside(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0)) ==1){kCB_CBandPID=1;}else{kCB_CBandPID=0;}

			//only MWPC and CB
			if(GetRootinos()->GetDetectors(0)==5 && cutcb->IsInside(GetRootinos()->Particle(0).E(),0.5*(GetRootinos()->GetMWPC1Energy(0)+GetRootinos()->GetMWPC0Energy(0))*TMath::Sin(GetRootinos()->GetTheta(0))) ==1){kCB_CBandMWPC=1;}else{kCB_CBandMWPC=0;}

			//TAPS+Vetos
			if(GetRootinos()->HasTAPS(0) && cuttaps->IsInside(GetRootinos()->Particle(0).E(),GetRootinos()->GetVetoEnergy(0))==1){kTAPS=1;}else{kTAPS=0;}


			//Cut on Clustersize
			Int_t cluster_size_neutron = GetPhotons()->GetClusterSize(0);
			Int_t cluster_size_piplus = GetRootinos()->GetClusterSize(0);


			if(GetPhotons()->HasCB(0))clustersize_cospiplus_energy_collerated_neutron_CB->Fill(cluster_size_neutron,cospiplus,beamphoton1E,time);
			if(GetPhotons()->HasTAPS(0))clustersize_cospiplus_energy_collerated_neutron_TAPS->Fill(cluster_size_neutron,cospiplus,beamphoton1E,time);

			if(GetRootinos()->HasCB(0))clustersize_cospiplus_energy_collerated_piplus_CB->Fill(cluster_size_piplus,cospiplus,beamphoton1E,time);
			if(GetRootinos()->HasTAPS(0))clustersize_cospiplus_energy_collerated_piplus_TAPS->Fill(cluster_size_piplus,cospiplus,beamphoton1E,time);


			if(cluster_size_neutron>4) continue; //no neutron candidate


/*
			if((kCB_CBandMWPC==1 || kCB_all ==1 || kCB_CBandPID==1) && kPSA_neutron==1 && kTOF_neutron==1){
				//NUMBER OF PIOs
				if(HistoManu::IsPrompt(time)){anzahlpiplus_prompt++; cutsfullfilled=kTRUE;}
				if(HistoManu::IsRandom(time)){anzahlpiplus_sideband++;}

				//some check histograms concerning neutron angle and clustersizes
				thetaneutron_cospiplus_energy_collerated->Fill(GetRootinos()->GetTheta(0)*TMath::RadToDeg(),cospiplus,beamphoton1E,time);
				clustersize_cospiplus_energy_collerated_piplus->Fill(GetRootinos()->GetClusterSize(0),cospiplus,beamphoton1E,time);
				clustersize_cospiplus_energy_collerated_neutron->Fill(GetPhotons()->GetClusterSize(0),cospiplus,beamphoton1E,time);
			
				cosverteilung_collerated_neutron->Fill(cospiplus,beamphoton1E,time);
				if(helicity==1 || helicity==0){cosverteilung_forE_collerated_neutron->Fill(cospiplus,beamphoton1E,time);}
	
				if(cbtrigger==kTRUE || GetScalers()->GetNEntries()==0){
	
					h_energy_sum_piplus->Fill(GetTrigger()->GetEnergySum(),time);
					//cosmeson_beamenergy_energysum_rek_neutron->Fill(cospiplus,beamphoton1E,GetTrigger()->GetEnergySum(),time);
					//cosmeson_beamenergy_rek_3ped->Fill(beamphoton1E,cospiplus,time);
				}
						
		

				//starting histogram for extraction of G and E
				if(pt>0 ){

					if(helicity==1){
						cos_beam_hel1_targetplus_neutron->Fill(cospiplus,beamphoton1E,time);
						cos_beam_hel1_targetplus_neutron_pc->Fill(cospiplus,beamphoton1E,time,circpol);
						cos_beam_hel1_targetplus_neutron_pt->Fill(cospiplus,beamphoton1E,time,pt);
					}

					if(helicity==0){
						cos_beam_hel0_targetplus_neutron->Fill(cospiplus,beamphoton1E,time);
						cos_beam_hel0_targetplus_neutron_pc->Fill(cospiplus,beamphoton1E,time,circpol);
						cos_beam_hel0_targetplus_neutron_pt->Fill(cospiplus,beamphoton1E,time,pt);
					}
		
					if(planesetting=="PARA"){
						kristallminus_targetplus_collerated_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time);
						kristallminus_targetplus_collerated_pb_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
						kristallminus_targetplus_collerated_pt_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);
					}
	
					if(planesetting=="PERP"){
						kristallplus_targetplus_collerated_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time);
						kristallplus_targetplus_collerated_pb_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
						kristallplus_targetplus_collerated_pt_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);	
					}
				}

				if(pt<0 ){

					if(helicity==1){
						cos_beam_hel1_targetminus_neutron->Fill(cospiplus,beamphoton1E,time);
						cos_beam_hel1_targetminus_neutron_pc->Fill(cospiplus,beamphoton1E,time,circpol);
						cos_beam_hel1_targetminus_neutron_pt->Fill(cospiplus,beamphoton1E,time,pt);
					}

					if(helicity==0){
						cos_beam_hel0_targetminus_neutron->Fill(cospiplus,beamphoton1E,time);
						cos_beam_hel0_targetminus_neutron_pc->Fill(cospiplus,beamphoton1E,time,circpol);
						cos_beam_hel0_targetminus_neutron_pt->Fill(cospiplus,beamphoton1E,time,pt);
					}

					if(planesetting=="PARA"){
						kristallminus_targetminus_collerated_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time);
						kristallminus_targetminus_collerated_pb_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
						kristallminus_targetminus_collerated_pt_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);
					}
	
					if(planesetting=="PERP"){
						kristallplus_targetminus_collerated_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time);
						kristallplus_targetminus_collerated_pb_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
						kristallplus_targetminus_collerated_pt_neutron->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);
				
					}
	
				}
			}		*/				
// 				
		}//tagger loop

	}//identified neutron

// 			//without any cuts!
// 
// 			MM->Fill(missingmass,time);
// 			IM->Fill(inv,time);
// 
// 
// 			if(GetRootinos()->GetNParticles()==1){
// 
// 				if(!(thetadiff > (unten_theta) && thetadiff < (oben_theta)) ||  !(phi > (unten_copl) && phi < (oben_copl))){continue;} 
// 
// 			}
// 
// 
// 			//with cuts
// 		
// 			if(((inv > (unten_inv) && inv < (oben_inv)))){
// 
// 				MM_all->Fill(missingmass,time);
// 				missingmassverteilung_collerated->Fill(missingmass,beamphoton1E,time);
// 				if(helicity==1 || helicity==0){missingmassverteilung_forE_collerated->Fill(missingmass,beamphoton1E,time);}
// 				missingmassverteilung_cospiplus_collerated->Fill(missingmass,cospiplus,beamphoton1E,time);
// 			}
// 	
// 			if(((missingmass > (unten_mass) && missingmass < (oben_mass)))){
// 
// 				IM_all->Fill(inv,time);
// 				invmassverteilung_collerated->Fill(inv,beamphoton1E,time);
// 				invmassverteilung_cospiplus_collerated->Fill(inv,cospiplus,beamphoton1E,time);
// 			}
// 
// 			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))&&GetTrigger()->GetNErrors()!=0){
// 				events_witherror->Fill(phipiplus,cospiplus, beamphoton1E);
// 			}
// 
// 			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){
// 				events_all->Fill(phipiplus,cospiplus, beamphoton1E);
// 			}
// 
// 			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){
// 
// 	
// 				if(GetPhotons()->HasTAPS(0)){
// 					Check_TAPS_TOF_photon_2_3ped->Fill(((GetPhotons()->GetTime(0)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(0).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(0).E(),time);
// 				}
// 
// 				if(GetPhotons()->HasTAPS(1)){
// 					Check_TAPS_TOF_photon_2_3ped->Fill(((GetPhotons()->GetTime(1)-GetTagger()->GetTaggedTime(j))/GetPhotons()->Particle(1).Vect().Mag())+1/0.299792458,GetPhotons()->Particle(1).E(),time);
// 				}
// 		
// 				cosverteilung_collerated->Fill(cospiplus,beamphoton1E,time);
// 
// 				if(helicity==1 || helicity==0){cosverteilung_forE_collerated->Fill(cospiplus,beamphoton1E,time);}
// 
// 				if(cbtrigger==kTRUE || GetScalers()->GetNEntries()==0){
// 
// 					h_energy_sum_piplus_2_3ped->Fill(GetTrigger()->GetEnergySum(),time);
// 					//cosmeson_beamenergy_energysum_rek->Fill(cospiplus,beamphoton1E,GetTrigger()->GetEnergySum(),time);
// 					//cosmeson_beamenergy_rek_2_3ped->Fill(beamphoton1E,cospiplus,time);
// 				}
// 						
// 				openingangle_gammatogamma_energy_allcuts->Fill(opening_angle_gammatogamma,beamphoton1E,time);
// 
// 
// 				if(pt>0 ){
// 
// 					if(helicity==1){
// 						cos_beam_hel1_targetplus->Fill(cospiplus,beamphoton1E,time);
// 						cos_beam_hel1_targetplus_pc->Fill(cospiplus,beamphoton1E,time,circpol);
// 						cos_beam_hel1_targetplus_pt->Fill(cospiplus,beamphoton1E,time,pt);
// 					}
// 
// 					if(helicity==0){
// 						cos_beam_hel0_targetplus->Fill(cospiplus,beamphoton1E,time);
// 						cos_beam_hel0_targetplus_pc->Fill(cospiplus,beamphoton1E,time,circpol);
// 						cos_beam_hel0_targetplus_pt->Fill(cospiplus,beamphoton1E,time,pt);
// 					}
// 		
// 					if(planesetting=="PARA"){
// 						kristallminus_targetplus_collerated->Fill(phipiplus,cospiplus,beamphoton1E,time);
// 						kristallminus_targetplus_collerated_pb->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
// 						kristallminus_targetplus_collerated_pt->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);
// 					}
// 	
// 					if(planesetting=="PERP"){
// 						kristallplus_targetplus_collerated->Fill(phipiplus,cospiplus,beamphoton1E,time);
// 						kristallplus_targetplus_collerated_pb->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
// 						kristallplus_targetplus_collerated_pt->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);	
// 					}
// 				}
// 
// 				if(pt<0 ){
// 
// 					if(helicity==1){
// 						cos_beam_hel1_targetminus->Fill(cospiplus,beamphoton1E,time);
// 						cos_beam_hel1_targetminus_pc->Fill(cospiplus,beamphoton1E,time,circpol);
// 						cos_beam_hel1_targetminus_pt->Fill(cospiplus,beamphoton1E,time,pt);
// 					}
// 
// 					if(helicity==0){
// 						cos_beam_hel0_targetminus->Fill(cospiplus,beamphoton1E,time);
// 						cos_beam_hel0_targetminus_pc->Fill(cospiplus,beamphoton1E,time,circpol);
// 						cos_beam_hel0_targetminus_pt->Fill(cospiplus,beamphoton1E,time,pt);
// 					}
// 
// 					if(planesetting=="PARA"){
// 						kristallminus_targetminus_collerated->Fill(phipiplus,cospiplus,beamphoton1E,time);
// 						kristallminus_targetminus_collerated_pb->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
// 						kristallminus_targetminus_collerated_pt->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);
// 					}
// 	
// 					if(planesetting=="PERP"){
// 						kristallplus_targetminus_collerated->Fill(phipiplus,cospiplus,beamphoton1E,time);
// 						kristallplus_targetminus_collerated_pb->Fill(phipiplus,cospiplus,beamphoton1E,time,pb);
// 						kristallplus_targetminus_collerated_pt->Fill(phipiplus,cospiplus,beamphoton1E,time,pt);
// 				
// 					}
// 	
// 				}
// 
// 			}
// 
// 		}//GetTagger()
// 	}//2 GetPhotons()

	anzahlpiplus_prompt=anzahlpiplus_prompt-HistoManu::backgroundSubstractionFactor*anzahlpiplus_sideband;
	if(cutsfullfilled==kTRUE){h_anzahl_piplus->Fill(anzahlpiplus_prompt);}
}


void NPiPlusAnalyse::fOnEndProcessing() {

}


void	NPiPlusAnalyse::ProcessScalerRead()
{
    //time.ScalerReadCorrection(5);
}

Bool_t NPiPlusAnalyse::Write()
{

TDirectory* curDir8  = outputFile->mkdir("Observable_E");
curDir8->cd();
cos_beam_hel1_targetminus_neutron->Write();
cos_beam_hel1_targetplus_neutron->Write();
cos_beam_hel0_targetminus_neutron->Write();
cos_beam_hel0_targetplus_neutron->Write();

cos_beam_hel1_targetminus_neutron_pc->Write();
cos_beam_hel1_targetplus_neutron_pc->Write();
cos_beam_hel0_targetminus_neutron_pc->Write();
cos_beam_hel0_targetplus_neutron_pc->Write();

cos_beam_hel1_targetminus_neutron_pt->Write();
cos_beam_hel1_targetplus_neutron_pt->Write();
cos_beam_hel0_targetminus_neutron_pt->Write();
cos_beam_hel0_targetplus_neutron_pt->Write();

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
cosverteilung_forE_collerated_neutron->Write();

TDirectory* curDir1;
if(runNumber>=1363 && runNumber<=3374 && GetScalers()->GetNScalers()!=0){//for carbon
	polsetting=NPiPlusAnalyse::poledge(inputFile);
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
h_energy_sum_piplus_2_3ped->Write();
h_energy_sum_piplus->Write();
test->Write();


/*
beamenergy_gen->Write();
cosmeson_beamenergy_energysum_rek->Write();
cosmeson_beamenergy_energysum_monte->Write();
cosmeson_beamenergy_energysum_rek_neutron->Write();
cosmeson_beamenergy_monte->Write();
cosmeson_beamenergy_rek_2_3ped->Write();
cosmeson_beamenergy_rek_3ped->Write();
*/
curDir1->cd();
TDirectory* curDir3  = curDir1->mkdir("Selektion");
curDir3->cd();
time_prompt->Write();
time1->Write();
time_side->Write();

test->Write();


curDir1->cd();
TDirectory* curDir5  = curDir1->mkdir("Selektion_withNeutron");
curDir5->cd();


h_anzahl_piplus->Write();
//CB over E
Check_CBdE_E_nocuts_all->Write();
Check_CBdE_E_nocuts_CBandPID->Write();
Check_CBdE_E_nocuts_CBandMWPC->Write();

Check_TAPSdE_E_nocuts->Write();

//TOF
Check_TAPS_TOF_neutron->Write();
Check_TAPS_TOF_piplus->Write();

Check_TAPS_energy_time_neutron->Write();
Check_TAPS_energy_time_piplus->Write();


//PSA
Check_PSA_piplus->Write();
Check_PSA_neutron->Write();
Check_PSA_together->Write();

coplanarity_withoutcuts->Write();

clustersize_cospiplus_energy_collerated_piplus_CB->Write();
clustersize_cospiplus_energy_collerated_piplus_TAPS->Write();
clustersize_cospiplus_energy_collerated_neutron_CB->Write();
clustersize_cospiplus_energy_collerated_neutron_TAPS->Write();
/*
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



cosverteilung_collerated->Write();
events_witherror->Write();
events_all->Write();

curDir1->cd();
TDirectory* curDir4  = curDir1->mkdir("Allgemein_withProton");
curDir4->cd();
kristallminus_targetplus_collerated_neutron->Write();
kristallminus_targetminus_collerated_neutron->Write();
kristallplus_targetplus_collerated_neutron->Write();
kristallplus_targetminus_collerated_neutron->Write();

kristallminus_targetplus_collerated_pb_neutron->Write();
kristallminus_targetminus_collerated_pb_neutron->Write();
kristallplus_targetplus_collerated_pb_neutron->Write();
kristallplus_targetminus_collerated_pb_neutron->Write();

kristallminus_targetplus_collerated_pt_neutron->Write();
kristallminus_targetminus_collerated_pt_neutron->Write();
kristallplus_targetplus_collerated_pt_neutron->Write();
kristallplus_targetminus_collerated_pt_neutron->Write();
thetaneutron_cospiplus_energy_collerated_withoutcuts->Write();

cosverteilung_collerated_neutron->Write();

clustersize_cospiplus_energy_collerated_piplus->Write();
clustersize_cospiplus_energy_collerated_neutron->Write();
thetaneutron_cospiplus_energy_collerated->Write();*/


outputFile->Close();
return 0;
}
