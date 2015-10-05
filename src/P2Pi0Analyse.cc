#include "P2Pi0Analyse.h"

P2Pi0Analyse::P2Pi0Analyse()
{ 
time_prompt = new HistoManu("time_prompt", 	"time_prompt", 	1400, -700, 700);
time1= new HistoManu("time1", 	"time1", 	1400, -700, 700);
time_side 	= new HistoManu("time_side", 	"time_side", 	1400, -700, 700);
   
 
poltable_energy          = new TH1F("poltable_energy",         "poltable_energy",           352,   0, 1448);
poltable_energy_weight          = new TH1F("poltable_energy_weight",         "poltable_energy_weight",           352,   0, 1448);

test = new TH1F("test","test",1100,-5.5,5.5);

//TEst for kinematic fit
pulls_4g_CB = new HistoManu2("pulls_4g_CB","pulls_4g_CB", 200, -10, 10, 3, 0, 3);
pulls_proton_CB = new HistoManu2("pulls_proton_CB","pulls_proton_CB", 200, -10, 10, 3, 0, 3);

pulls_4g_TAPS = new HistoManu2("pulls_4g_TAPS","pulls_4g_TAPS", 200, -10, 10, 3, 0, 3);
pulls_proton_TAPS = new HistoManu2("pulls_proton_TAPS","pulls_proton_TAPS", 200, -10, 10, 3, 0, 3);

pulls_beam = new HistoManu2("pulls_beam","pulls_beam", 200, -10, 10, 3, 0, 3);
CL  = new HistoManu("CL ", "CL;CL ",  100, 0, 1);
invmass_diff_rec_kin  = new HistoManu("invmass_diff_rec_kin ", "invmass_diff_rec_kin;m_{2#gamma,rec}-m_{2#gamma,fit} ",  200, -100, 100);
thetaproton_fit_vs_rek_TAPS= new HistoManu("thetaproton_fit_vs_rek_TAPS", "Polar angle (Signal); #theta_{rek}-#theta_{meas} [deg]", 400,-200,200);
thetaproton_fit_vs_rek_CB= new HistoManu("thetaproton_fit_vs_rek_CB", "Polar angle (Signal); #theta_{rek}-#theta_{meas} [deg]", 400,-200,200);

coplanarity_fit_vs_rek = new HistoManu("coplanarity_fit_vs_rek", "Azimutwinkel (Signal); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);

openangle_pi01_to_pi02_energy = new HistoManu3("openangle_pi01_to_pi02_energy","openangle_pi01_to_pi02_energy;#angle #gamma1,#gamma2;#angle #gamma3,#gamma4;E_{beam} [MeV]",180,0,180,180,0,180,37,200,1421);

events_all = new TH3F("events_all","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);
events_witherror = new TH3F("events_witherror","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//.......................................histograms for  Coplanarity.............................

coplanarity_collerated = new HistoManu("coplanarity_collerated", "Azimutwinkel (Signal); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_mass_theta_collerated = new HistoManu("coplanarity_mass_theta_collerated", "Azimutwinkel (Cut: PolarWinkel+MissingMass(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_mass_theta_inv_collerated = new HistoManu("coplanarity_mass_theta_inv_collerated", "Azimutwinkel (Cut: All(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_mass_collerated = new HistoManu("coplanarity_mass_collerated", "Azimutwinkel (Cut: MissingMass(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_theta_collerated = new HistoManu("coplanarity_theta_collerated", "Azimutwinkel (Cut: PolarWinkel+MissingMass(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_inv_collerated = new HistoManu("coplanarity_inv_collerated", "Azimutwinkel (Cut: invariante Masse(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);


//.........................histograms for  ThetaProton Polarwinkel...........................

thetaproton_collerated = new HistoManu("thetaproton_collerated", "Polarwinkel (Signal); #theta_{rek}-#theta_{meas} [deg]", 400,-200,200);
thetaproton_mass_copl_collerated = new HistoManu("thetaproton_mass_copl_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_copl_inv_collerated = new HistoManu("thetaproton_mass_copl_inv_collerated", "Polarwinkel (Cut: All(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_collerated = new HistoManu("thetaproton_mass_collerated", "Polarwinkel (Cut: Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_copl_collerated = new HistoManu("thetaproton_copl_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_inv_collerated = new HistoManu("thetaproton_inv_collerated", "Polarwinkel (Cut: invariante Masse(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);



//detektiert in TAPS
thetaproton_taps_collerated = new HistoManu("thetaproton_taps_collerated", "Polarwinkel (Signal); #theta_{rek}-#theta_{meas} [deg]", 400,-200,200);
thetaproton_mass_copl_taps_collerated = new HistoManu("thetaproton_mass_copl_taps_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_copl_inv_taps_collerated = new HistoManu("thetaproton_mass_copl_inv_taps_collerated", "Polarwinkel (Cut: All(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_taps_collerated = new HistoManu("thetaproton_mass_taps_collerated", "Polarwinkel (Cut: Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_copl_taps_collerated = new HistoManu("thetaproton_copl_taps_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_inv_taps_collerated = new HistoManu("thetaproton_inv_taps_collerated", "Polarwinkel (Cut: invariante Masse(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);


//............................histograms for  Missing Mass................................

missingmass_collerated = new HistoManu("missingmass_collerated", "Missing Mass (Signal);m_{mm} [MeV]", 500,0,2200);
missingmass_collerated_proton = new HistoManu("missingmass_collerated", "Missing Mass (Signal);m_{mm} [MeV]", 500,0,2200);

missingmass_inv_collerated = new HistoManu("missingmass_inv_collerated", "Missing Mass (Cut: invariante Masse(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_inv_collerated_proton = new HistoManu("missingmass_inv_collerated_proton", "Missing Mass (Cut: invariante Masse(Signal));m_{mm} [MeV]", 500,0,2200);

missingmass_theta_copl_collerated = new HistoManu("missingmass_theta_copl_collerated", "Missing Mass (Cut: Polar+Azimut-Winkel(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_theta_copl_inv_collerated = new HistoManu("missingmass_theta_copl_inv_collerated", "Missing Mass (Cut: All(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_theta_collerated = new HistoManu("missingmass_theta_collerated", "Missing Mass (Cut: Polarwinkel(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_copl_collerated = new HistoManu("missingmass_copl_collerated", "Missing Mass (Cut: Azimutalwinkel(Signal));m_{mm} [MeV]", 500,0,2200);

//.......................histograms for  invariante Masse.....................................

massesumme_mass_theta_copl_inv_beam_collerated = new HistoManu2("massesumme_mass_theta_copl_inv_beam_collerated","inv Mass in dependency of E_{beam}; m_{#gamma #gamma} [MeV]; E^{rec}_{#gamma} [MeV]",400,0,600,37,200,1421);
massesumme_collerated = new HistoManu("massesumme_collerated", "Invariante Masse von der Summe (Signal); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_collerated_proton = new HistoManu("massesumme_collerated_proton", "Invariante Masse von der Summe (Signal); m_{#gamma #gamma} [MeV]", 400,0,600);

massesumme_mass_theta_copl_collerated = new HistoManu("massesumme_mass_theta_copl_collerated", "Invariante Masse von der Summe (Cut: All(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_theta_copl_inv1_collerated = new HistoManu("massesumme_mass_theta_copl_inv1_collerated", "Invariante Masse von der Summe (Cut: All+Cut auf invariante Masse 1(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_theta_copl_inv_collerated = new HistoManu("massesumme_mass_theta_copl_inv_collerated", "Invariante Masse von der Summe (Cut: All+Cut auf invariante Masse 1(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_copl_collerated = new HistoManu("massesumme_mass_copl_collerated", "Invariante Masse von der Summe (Cut: Azimutwinkel+Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_collerated = new HistoManu("massesumme_mass_collerated", "Invariante Masse von der Summe (Cut: Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_copl_collerated = new HistoManu("massesumme_copl_collerated", "Invariante Masse von der Summe (Cut: Azimutwinkel(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_theta_collerated = new HistoManu("massesumme_theta_collerated", "Invariante Masse von der Summe (Cut: Polarwinkel(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_collerated_proton = new HistoManu("massesumme_mass_collerated_proton", "Invariante Masse von der Summe (Cut: Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);


massegegenmasse_collerated = new HistoManu2("massegegenmasse_collerated","invariante Massen gegeneinander (Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

massegegenmasse_collerated_proton = new HistoManu2("massegegenmasse_collerated_proton","invariante Massen gegeneinander (Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

massegegenmasse_mass_copl_collerated = new HistoManu2("massegegenmasse_mass_copl_collerated","invariante Massen gegeneinander(Cut: Missing Mass+Azimutalwinkel(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_mass_theta_copl_collerated = new HistoManu2("massegegenmasse_mass_theta_copl_collerated","invariante Massen gegeneinander(Cut: All(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_copl_collerated = new HistoManu2("massegegenmasse_copl_collerated","invariante Massen gegeneinander(Cut: Azimutalwinkel(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_theta_collerated = new HistoManu2("massegegenmasse_theta_collerated","invariante Massen gegeneinander(Cut: Polarwinkel(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

massegegenmasse_mass_collerated = new HistoManu2("massegegenmasse_mass_collerated","invariante Massen gegeneinander(Cut: Missing Mass(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_mass_collerated_proton = new HistoManu2("massegegenmasse_mass_collerated_proton","invariante Massen gegeneinander(Cut: Missing Mass(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

massegegenmasse_allcuts_collerated = new HistoManu2("massegegenmasse_allcuts_collerated","invariante Massen gegeneinander(Cut: all); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

massegegenmasse_allcutsandinv_collerated = new HistoManu2("massegegenmasse_allcutsandinv_collerated","invariante Massen gegeneinander(Cut: all+invmass); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);


//ANNNNNNNNNAAAAAAALLLLLLLLLLYYYYYYYYYYYSSSSSSSSSSSSIIIIIIIIISSSSSSSSSSSSSSSSSSS

cosverteilung_collerated = new HistoManu2("cosverteilung_collerated", "Cos-Verteilung; cos(#theta_{2#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);
cosppi0verteilung_collerated = new HistoManu2("cosppi0verteilung_collerated", "Cos-Verteilung; cos(#theta_{p#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);

invverteilung_collerated = new HistoManu2("invverteilung_collerated", "Cos-Verteilung; m_{2#pi}[MeV];E_{beam} [MeV]", 100,300,800,37,200,1421);
invppi0verteilung_collerated = new HistoManu2("invppi0verteilung_collerated", "Cos-Verteilung; m_{p#pi}[MeV];E_{beam} [MeV]", 100,1100,1600,37,200,1421);

cosverteilung_thetaproton_collerated = new HistoManu2("cosverteilung_thetaproton_collerated", "Cos-Verteilung; cos(#theta_{2#pi});#theta_{proton} [deg]", 72,-1,1,180,0,180);

thetaverteilung_collerated_taps = new HistoManu2("thetaverteilung_collerated_taps", "Polar-Verteilung;#theta_{rek}-#theta_{meas} [deg]", 400,-200,200,37,200,1421);
thetaverteilung_collerated = new HistoManu2("thetaverteilung_collerated", "Polar-Verteilung;#theta_{rek}-#theta_{meas} [deg]", 400,-200,200,37,200,1421);

coplanarityverteilung_collerated = new HistoManu2("coplanarityverteilung_collerated", "Polar-Verteilung;#phi_{2#pi}-#phi_{p}[deg]", 400,0,360,37,200,1421);

missingmassverteilung_collerated = new HistoManu2("missingmassverteilung_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,37,200,1421);
missingmassverteilung_collerated_proton = new HistoManu2("missingmassverteilung_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,37,200,1421);

//histograms for  2pi0 System with Cos-Verteilung___________________________________________________
/*
kristallminus_targetplus_collerated = new HistoManu3("kristallminus_targetplus_collerated","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated = new HistoManu3("kristallminus_targetminus_collerated","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated = new HistoManu3("kristallplus_targetplus_collerated","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated = new HistoManu3("kristallplus_targetminus_collerated","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_collerated_pb = new HistoManu3("kristallminus_targetplus_collerated_pb","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pb = new HistoManu3("kristallminus_targetminus_collerated_pb","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pb = new HistoManu3("kristallplus_targetplus_collerated_pb","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pb = new HistoManu3("kristallplus_targetminus_collerated_pb","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_collerated_pt = new HistoManu3("kristallminus_targetplus_collerated_pt","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pt = new HistoManu3("kristallminus_targetminus_collerated_pt","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pt = new HistoManu3("kristallplus_targetplus_collerated_pt","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pt = new HistoManu3("kristallplus_targetminus_collerated_pt","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms for  ppi0 System  with Cos-Verteilung__________________________________________________________

kristallminus_targetplus_ppi0_collerated = new HistoManu3("kristallminus_targetplus_ppi0_collerated","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_ppi0_collerated = new HistoManu3("kristallminus_targetminus_ppi0_collerated","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_ppi0_collerated = new HistoManu3("kristallplus_targetplus_ppi0_collerated","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_ppi0_collerated = new HistoManu3("kristallplus_targetminus_ppi0_collerated","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_ppi0_collerated_pb = new HistoManu3("kristallminus_targetplus_ppi0_collerated_pb","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_ppi0_collerated_pb = new HistoManu3("kristallminus_targetminus_ppi0_collerated_pb","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_ppi0_collerated_pb = new HistoManu3("kristallplus_targetplus_ppi0_collerated_pb","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_ppi0_collerated_pb = new HistoManu3("kristallplus_targetminus_ppi0_collerated_pb","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_ppi0_collerated_pt = new HistoManu3("kristallminus_targetplus_ppi0_collerated_pt","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_ppi0_collerated_pt = new HistoManu3("kristallminus_targetminus_ppi0_collerated_pt","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_ppi0_collerated_pt = new HistoManu3("kristallplus_targetplus_ppi0_collerated_pt","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_ppi0_collerated_pt = new HistoManu3("kristallplus_targetminus_ppi0_collerated_pt","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms for  2pi0 System with invM-Verteilung___________________________________________________

kristallminus_targetplus_im_collerated = new HistoManu3("kristallminus_targetplus_im_collerated","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallminus_targetminus_im_collerated = new HistoManu3("kristallminus_targetminus_im_collerated","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetplus_im_collerated = new HistoManu3("kristallplus_targetplus_im_collerated","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetminus_im_collerated = new HistoManu3("kristallplus_targetminus_im_collerated","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_im_collerated_pb = new HistoManu3("kristallminus_targetplus_im_collerated_pb","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallminus_targetminus_im_collerated_pb = new HistoManu3("kristallminus_targetminus_im_collerated_pb","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetplus_im_collerated_pb = new HistoManu3("kristallplus_targetplus_im_collerated_pb","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetminus_im_collerated_pb = new HistoManu3("kristallplus_targetminus_im_collerated_pb","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_im_collerated_pt = new HistoManu3("kristallminus_targetplus_im_collerated_pt","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallminus_targetminus_im_collerated_pt = new HistoManu3("kristallminus_targetminus_im_collerated_pt","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetplus_im_collerated_pt = new HistoManu3("kristallplus_targetplus_im_collerated_pt","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetminus_im_collerated_pt = new HistoManu3("kristallplus_targetminus_im_collerated_pt","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

//histograms for  ppi0 System with invM-distribution___________________________________________________________________

kristallminus_targetplus_im_ppi0_collerated = new HistoManu3("kristallminus_targetplus_im_ppi0_collerated","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallminus_targetminus_im_ppi0_collerated = new HistoManu3("kristallminus_targetminus_im_ppi0_collerated","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetplus_im_ppi0_collerated = new HistoManu3("kristallplus_targetplus_im_ppi0_collerated","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetminus_im_ppi0_collerated = new HistoManu3("kristallplus_targetminus_im_ppi0_collerated","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_im_ppi0_collerated_pb = new HistoManu3("kristallminus_targetplus_im_ppi0_collerated_pb","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallminus_targetminus_im_ppi0_collerated_pb = new HistoManu3("kristallminus_targetminus_im_ppi0_collerated_pb","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetplus_im_ppi0_collerated_pb = new HistoManu3("kristallplus_targetplus_im_ppi0_collerated_pb","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetminus_im_ppi0_collerated_pb = new HistoManu3("kristallplus_targetminus_im_ppi0_collerated_pb","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_im_ppi0_collerated_pt = new HistoManu3("kristallminus_targetplus_im_ppi0_collerated_pt","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallminus_targetminus_im_ppi0_collerated_pt = new HistoManu3("kristallminus_targetminus_im_ppi0_collerated_pt","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetplus_im_ppi0_collerated_pt = new HistoManu3("kristallplus_targetplus_im_ppi0_collerated_pt","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetminus_im_ppi0_collerated_pt = new HistoManu3("kristallplus_targetminus_im_ppi0_collerated_pt","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);
*/
//PROTON IDENTIFIED

cosverteilung_collerated_proton = new HistoManu2("cosverteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{2#pi});E_{beam} [MeV]", 18,-1,1,30,220,1420);
/*
cosppi0verteilung_collerated_proton = new HistoManu2("cosppi0verteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{p#pi});E_{beam} [MeV]", 72,-1,1,37,200,1421);

invverteilung_collerated_proton = new HistoManu2("invverteilung_collerated_proton", "Cos-Verteilung; m_{2#pi}[MeV];E_{beam} [MeV]", 100,300,800,37,200,1421);

invppi0verteilung_collerated_proton = new HistoManu2("invppi0verteilung_collerated_proton", "Cos-Verteilung; m_{p#pi}[MeV];E_{beam} [MeV]", 100,1100,1600,37,200,1421);

cosverteilung_thetaproton_collerated_proton = new HistoManu2("cosverteilung_thetaproton_collerated_proton", "Cos-Verteilung; cos(#theta_{2#pi});#theta_{proton} [deg]", 72,-1,1,180,0,180);

//histograms for  2pi0 System with Cos-distribution___________________________________________________

kristallminus_targetplus_collerated_proton = new HistoManu3("kristallminus_targetplus_collerated_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_proton = new HistoManu3("kristallminus_targetminus_collerated_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_proton = new HistoManu3("kristallplus_targetplus_collerated_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_proton = new HistoManu3("kristallplus_targetminus_collerated_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_collerated_pb_proton = new HistoManu3("kristallminus_targetplus_collerated_pb_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pb_proton = new HistoManu3("kristallminus_targetminus_collerated_pb_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pb_proton = new HistoManu3("kristallplus_targetplus_collerated_pb_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pb_proton = new HistoManu3("kristallplus_targetminus_collerated_pb_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_collerated_pt_proton = new HistoManu3("kristallminus_targetplus_collerated_pt_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_collerated_pt_proton = new HistoManu3("kristallminus_targetminus_collerated_pt_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_collerated_pt_proton = new HistoManu3("kristallplus_targetplus_collerated_pt_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_collerated_pt_proton = new HistoManu3("kristallplus_targetminus_collerated_pt_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);


//histograms for  ppi0 System  with Cos-distribution__________________________________________________________

kristallminus_targetplus_ppi0_collerated_proton = new HistoManu3("kristallminus_targetplus_ppi0_collerated_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_ppi0_collerated_proton = new HistoManu3("kristallminus_targetminus_ppi0_collerated_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_ppi0_collerated_proton = new HistoManu3("kristallplus_targetplus_ppi0_collerated_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_ppi0_collerated_proton = new HistoManu3("kristallplus_targetminus_ppi0_collerated_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_ppi0_collerated_pb_proton = new HistoManu3("kristallminus_targetplus_ppi0_collerated_pb_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_ppi0_collerated_pb_proton = new HistoManu3("kristallminus_targetminus_ppi0_collerated_pb_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_ppi0_collerated_pb_proton = new HistoManu3("kristallplus_targetplus_ppi0_collerated_pb_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_ppi0_collerated_pb_proton = new HistoManu3("kristallplus_targetminus_ppi0_collerated_pb_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_ppi0_collerated_pt_proton = new HistoManu3("kristallminus_targetplus_ppi0_collerated_pt_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallminus_targetminus_ppi0_collerated_pt_proton = new HistoManu3("kristallminus_targetminus_ppi0_collerated_pt_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetplus_ppi0_collerated_pt_proton = new HistoManu3("kristallplus_targetplus_ppi0_collerated_pt_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

kristallplus_targetminus_ppi0_collerated_pt_proton = new HistoManu3("kristallplus_targetminus_ppi0_collerated_pt_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,37,200,1421);

//histograms for  2pi0 System mit invM-distribution___________________________________________________

kristallminus_targetplus_im_collerated_proton = new HistoManu3("kristallminus_targetplus_im_collerated_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallminus_targetminus_im_collerated_proton = new HistoManu3("kristallminus_targetminus_im_collerated_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetplus_im_collerated_proton = new HistoManu3("kristallplus_targetplus_im_collerated_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetminus_im_collerated_proton = new HistoManu3("kristallplus_targetminus_im_collerated_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_im_collerated_pb_proton = new HistoManu3("kristallminus_targetplus_im_collerated_pb_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallminus_targetminus_im_collerated_pb_proton = new HistoManu3("kristallminus_targetminus_im_collerated_pb_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetplus_im_collerated_pb_proton = new HistoManu3("kristallplus_targetplus_im_collerated_pb_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetminus_im_collerated_pb_proton = new HistoManu3("kristallplus_targetminus_im_collerated_pb_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_im_collerated_pt_proton = new HistoManu3("kristallminus_targetplus_im_collerated_pt_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallminus_targetminus_im_collerated_pt_proton = new HistoManu3("kristallminus_targetminus_im_collerated_pt_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetplus_im_collerated_pt_proton = new HistoManu3("kristallplus_targetplus_im_collerated_pt_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

kristallplus_targetminus_im_collerated_pt_proton = new HistoManu3("kristallplus_targetminus_im_collerated_pt_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,37,200,1421);

//histograms for  ppi0 System with invM-distribution___________________________________________________________________

kristallminus_targetplus_im_ppi0_collerated_proton = new HistoManu3("kristallminus_targetplus_im_ppi0_collerated_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallminus_targetminus_im_ppi0_collerated_proton = new HistoManu3("kristallminus_targetminus_im_ppi0_collerated_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetplus_im_ppi0_collerated_proton = new HistoManu3("kristallplus_targetplus_im_ppi0_collerated_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetminus_im_ppi0_collerated_proton = new HistoManu3("kristallplus_targetminus_im_ppi0_collerated_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

//histograms weighted with beam polarisation
kristallminus_targetplus_im_ppi0_collerated_pb_proton = new HistoManu3("kristallminus_targetplus_im_ppi0_collerated_pb_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallminus_targetminus_im_ppi0_collerated_pb_proton = new HistoManu3("kristallminus_targetminus_im_ppi0_collerated_pb_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetplus_im_ppi0_collerated_pb_proton = new HistoManu3("kristallplus_targetplus_im_ppi0_collerated_pb_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetminus_im_ppi0_collerated_pb_proton = new HistoManu3("kristallplus_targetminus_im_ppi0_collerated_pb_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

//histograms weighted with target polarisation
kristallminus_targetplus_im_ppi0_collerated_pt_proton = new HistoManu3("kristallminus_targetplus_im_ppi0_collerated_pt_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallminus_targetminus_im_ppi0_collerated_pt_proton = new HistoManu3("kristallminus_targetminus_im_ppi0_collerated_pt_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetplus_im_ppi0_collerated_pt_proton = new HistoManu3("kristallplus_targetplus_im_ppi0_collerated_pt_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);

kristallplus_targetminus_im_ppi0_collerated_pt_proton = new HistoManu3("kristallplus_targetminus_im_ppi0_collerated_pt_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,37,200,1421);
*/
h_energy_sum_2pi0_5ped = new HistoManu("h_energy_sum_2pi0_5ped ", "h_energy_sum_2pi0_5ped;E_{sum} [MeV]",  300, 0, 1557);

h_energy_taps_2pi0_5ped = new HistoManu2("h_energy_taps_2pi0_5ped ", "h_energy_taps_2pi0_5ped;E_{sum} [MeV] ",  900, 0, 900,300,0,1557);

h_energy_sum_2pi0_5ped_weightcos2pi0 = new HistoManu("h_energy_sum_2pi0_5ped_weightcos2pi0 ", "h_energy_sum_2pi0_5ped_weightcos2pi0;E_{sum} [MeV] ",  300, 0, 1557);

h_energy_sum_2pi0_5ped_weightcosppi0 = new HistoManu("h_energy_sum_2pi0_5ped_weightcosppi0 ", "h_energy_sum_2pi0_5ped_weightcosppi0;E_{sum} [MeV] ",  300, 0, 1557);


//generated
cos2pi0_beamphoton_monte = new TH2F("cos2pi0_beamphoton_monte", "cos2pi0_beamphoton_monte;cos #theta_{2#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420);

cosppi0_beamphoton_energysum_monte = new TH3F("cosppi0_beamphoton_energysum_monte", "cosppi0_beamphoton_energysum_monte;cos #theta_{p#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);

m2pi0_beamphoton_energysum_monte = new TH3F("m2pi0_beamphoton_energysum_monte", "m2pi0_beamphoton_energysum_monte;m_{2#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  70, 0, 700,30, 220, 1420 ,200, 0, 2000);

mppi0_beamphoton_energysum_monte = new TH3F("mppi0_beamphoton_energysum_monte", "mppi0_beamphoton_energysum_monte;m_{p#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  150, 0, 1500,30, 220, 1420 ,200, 0, 2000);

//reconstructed
cos2pi0_beamphoton_rek = new HistoManu2("cos2pi0_beamphoton_rek", "cos2pi0_beamphoton_rek;cos #theta_{2#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420);

cosppi0_beamphoton_energysum_rek = new HistoManu3("cosppi0_beamphoton_energysum_rek", "cosppi0_beamphoton_energysum_rek;cos #theta_{p#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);

m2pi0_beamphoton_energysum_rek = new HistoManu3("m2pi0_beamphoton_energysum_rek", "m2pi0_beamphoton_energysum_rek;m_{2#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  70, 0, 700,30, 220, 1420 ,200, 0, 2000);

mppi0_beamphoton_energysum_rek = new HistoManu3("mppi0_beamphoton_energysum_rek", "mppi0_beamphoton_energysum_rek;m_{p#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  150, 0, 1500,30, 220, 1420 ,200, 0, 2000);

h_bph_e_gen  = new TH1F("h_bph_e_gen ", "h_bph_e_gen;E_{#gamma} [MeV]",  300, 0, 1557);

triggertest = new TH1F("triggertest","triggerspattern",34,0,34);


Check_CBdE_E= new HistoManu2("Check_CBdE_E", "dE_E (all CB clusters compared to PID hits)", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E= new HistoManu2("Check_TAPSdE_E", "dE_E (all TAPS clusters compared to Veto hits)", 	400, 0, 400, 100, 0, 10);

//Initialize the fitter and the time background subtraction

fitter.AddConstraintsTotEnergy();
fitter.AddConstraintsTotMomentum();
fitter.AddConstraintsIM();
// // fitter.AddConstraintsMM();
// 
fitter1.AddConstraintsTotEnergy();
fitter1.AddConstraintsTotMomentum();
fitter1.AddConstraintsIM();
// // fitter1.AddConstraintsMM();
// 
fitter2.AddConstraintsTotEnergy();
fitter2.AddConstraintsTotMomentum();
fitter2.AddConstraintsIM();
fitter2.AddConstraintsMM();


HistoManu::InitCuts(-8, 8,100, 300);

}

P2Pi0Analyse::~P2Pi0Analyse()
{
}



Bool_t	P2Pi0Analyse::Start()
{
pt=P2Pi0Analyse::targetpol(inputFile);

//filenumber
TString* filename1 = new TString(inputFile->GetPath());
path11= new TString(filename1->Tokenize("_")->At(filename1->Tokenize("_")->GetEntries()-1)->GetName());
path11->Resize(path11->Length()-7);

//runnumber
stringstream ss(path11->Data());
ss >> runNumber;
if(runNumber>=1363 && runNumber<=3374)planesetting=P2Pi0Analyse::polplane(inputFile);
else{
if(GetLinpol()->GetPolarizationPlane()==1){planesetting="PERP";}
if(GetLinpol()->GetPolarizationPlane()==0){planesetting="PARA";}
}
calfile.open(Form("/disk/user/spieker/GoAT/txtfiles/liste_MWPC_prompt_2pi0_%i.txt",runNumber));
calfile1.open(Form("/disk/user/spieker/GoAT/txtfiles/liste_MWPC_side_2pi0_%i.txt",runNumber));


    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();


	return kTRUE;

}



void P2Pi0Analyse::fOnBeforeEventProcessing() {

}


void	P2Pi0Analyse::ProcessEvent()	
{

//for monte carlo
if(GetScalers()->GetNEntries()==0){

	TLorentzVector beam_mc = GetGeant()->GetBeam();
	Double_t e_beam_mc = beam_mc.E()*1000.;
	h_bph_e_gen->Fill(e_beam_mc);

	TLorentzVector proton_mc = GetGeant()->GetTrueVector(0);
	proton_mc.SetPxPyPzE(proton_mc.Px()*1000,proton_mc.Py()*1000,proton_mc.Pz()*1000,proton_mc.E()*1000);

	TLorentzVector g1_mc = GetGeant()->GetTrueVector(1);
	TLorentzVector g2_mc = GetGeant()->GetTrueVector(2);	
	TLorentzVector g3_mc = GetGeant()->GetTrueVector(3);
	TLorentzVector g4_mc = GetGeant()->GetTrueVector(4);
	
	TLorentzVector meson1_mc = g1_mc+g2_mc;
	TLorentzVector meson2_mc = g3_mc+g4_mc;

	meson1_mc.SetPxPyPzE(meson1_mc.Px()*1000,meson1_mc.Py()*1000,meson1_mc.Pz()*1000,meson1_mc.E()*1000);
	meson2_mc.SetPxPyPzE(meson2_mc.Px()*1000,meson2_mc.Py()*1000,meson2_mc.Pz()*1000,meson2_mc.E()*1000);

// 	Pi0Pi0-System
	TLorentzVector pi0pi0_4vektor_mc=meson1_mc+meson2_mc;
	TLorentzVector pi0pi0_4vektor_mc_boost = CMVector(pi0pi0_4vektor_mc, beam_mc, proton_mc);
	Double_t invariantmass_pi0pi0_mc = pi0pi0_4vektor_mc_boost.M();
	TVector3 pi0pi0_3vektor_boost_mc = pi0pi0_4vektor_mc_boost.Vect();
	Double_t cospi0pi0_mc = pi0pi0_3vektor_boost_mc.CosTheta();
	Double_t phimeson_mc = TMath::RadToDeg()*(pi0pi0_4vektor_mc_boost.Vect().Phi());
			
// 	pPi0-System for first Pion
	TLorentzVector ppi01_4vektor_mc = proton_mc+meson1_mc;
	TLorentzVector ppi01_4vektor_boost_mc = CMVector(ppi01_4vektor_mc, beam_mc, proton_mc);
	Double_t invariantmass_ppi01_mc = ppi01_4vektor_boost_mc.M();
	TVector3 ppi01_3vektor_boost_mc = ppi01_4vektor_boost_mc.Vect();
	Double_t cosppi01_mc = ppi01_3vektor_boost_mc.CosTheta();
	Double_t phippi01_mc = TMath::RadToDeg()*ppi01_3vektor_boost_mc.Phi();		
// 	pPi0-System for second Pion
	TLorentzVector ppi02_4vektor_mc = proton_mc+meson2_mc;
	TLorentzVector ppi02_4vektor_boost_mc = CMVector(ppi02_4vektor_mc, beam_mc, proton_mc);
	Double_t invariantmass_ppi02_mc = ppi02_4vektor_boost_mc.M();
	TVector3 ppi02_3vektor_boost_mc = ppi02_4vektor_boost_mc.Vect();
	Double_t cosppi02_mc = ppi02_3vektor_boost_mc.CosTheta();
	Double_t phippi02_mc = TMath::RadToDeg()*ppi02_3vektor_boost_mc.Phi();

	cos2pi0_beamphoton_monte->Fill(cospi0pi0_mc,e_beam_mc);
	m2pi0_beamphoton_energysum_monte->Fill(invariantmass_pi0pi0_mc,e_beam_mc,1000*GetGeant()->GetCBESum());

	cosppi0_beamphoton_energysum_monte->Fill(cosppi01_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
	mppi0_beamphoton_energysum_monte->Fill(invariantmass_ppi01_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
	cosppi0_beamphoton_energysum_monte->Fill(cosppi02_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
	mppi0_beamphoton_energysum_monte->Fill(invariantmass_ppi02_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
 }


if(GetPhotons()->GetNParticles()==4){
//::::::::::::::::::::::::::::::::::::::::::::::::::define photons::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
			
	TLorentzVector photon1_4vektor = GetPhotons()->Particle(0);
	TLorentzVector photon2_4vektor = GetPhotons()->Particle(1);
	TLorentzVector photon3_4vektor = GetPhotons()->Particle(2);
	TLorentzVector photon4_4vektor = GetPhotons()->Particle(3);

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::define 3 vector::::::::::::::::::::::::::::::::::::::::::::::::::
			
	TVector3 photon1_3vektor = photon1_4vektor.Vect();
	TVector3 photon2_3vektor = photon2_4vektor.Vect();
	TVector3 photon3_3vektor = photon3_4vektor.Vect();
	TVector3 photon4_3vektor = photon4_4vektor.Vect();
			
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::calculate inv mass:::::::::::::::::::::::::::::::::::::::::::::::::
			
	//Berechne invariante Masse von Photon 1 und 2			
	Double_t invariantemasse12 = invariantemasse(GetPhotons(),0,1);
	//Berechne invariante Masse von Photon 3 und 4	
	Double_t invariantemasse34 = invariantemasse(GetPhotons(),2,3);
		
	//Berechne invariante Masse von Photon 1 und 3
	Double_t invariantemasse13 = invariantemasse(GetPhotons(),0,2);
	//Berechne invariante Masse von Photon 2 und 4	
	Double_t invariantemasse24 = invariantemasse(GetPhotons(),1,3);

	//Berechne invariante Masse von Photon 1 und 4
	Double_t invariantemasse14 =invariantemasse(GetPhotons(),0,3);
	//Berechne invariante Masse von Photon 2 und 3	
	Double_t invariantemasse23 = invariantemasse(GetPhotons(),1,2);

// /::::::::::::::::::::::::::::::::::::::::::::::::::find correct combination:::::::::::::::::::::::::::::::::::::::::::::::::::::::	
// 
	Double_t invariantemasse1, invariantemasse2,Chi;

	Double_t pion1_E, pion2_E, pion1_m, pion2_m, pion1_pz, pion2_pz;
	Double_t photon1_E = photon1_4vektor.E();
	Double_t photon2_E = photon2_4vektor.E();
	Double_t photon3_E = photon3_4vektor.E();
	Double_t photon4_E = photon4_4vektor.E();
	Double_t photon1_pz = photon1_4vektor.Pz();
	Double_t photon2_pz = photon2_4vektor.Pz();
	Double_t photon3_pz = photon3_4vektor.Pz();
	Double_t photon4_pz = photon4_4vektor.Pz();

//_________________performe chi^2 for 2pi0 for all 4 gamma combination 

	Double_t Chi12 = ChiPionPion(GetPhotons(),invariantemasse12,0,1,invariantemasse34,2,3);
	Double_t Chi13 = ChiPionPion(GetPhotons(),invariantemasse13,0,2,invariantemasse24,1,3);
	Double_t Chi14 = ChiPionPion(GetPhotons(), invariantemasse14,0,3,invariantemasse23,1,2);

//____________________________set best combination as the pions 
	TLorentzVector pi01_4vect,pi02_4vect;

	if(Chi12 < Chi13 && Chi12 < Chi14){
		invariantemasse1 = invariantemasse12;
		pi01_4vect = GetPhotons()->Particle(0)+GetPhotons()->Particle(1);
		pion1_E = photon1_E + photon2_E ;
		pion1_pz = photon1_pz + photon2_pz;	
		invariantemasse2 = invariantemasse34;
		pi02_4vect = GetPhotons()->Particle(2)+GetPhotons()->Particle(3);
		pion2_E = photon3_E + photon4_E; 
		pion2_pz = photon3_pz + photon4_pz;
		Chi = Chi12;
	}
			
	if(Chi13 < Chi12 && Chi13 < Chi14){
		invariantemasse1 = invariantemasse13;
		pi01_4vect = GetPhotons()->Particle(0)+GetPhotons()->Particle(2);
		pion1_E = photon1_E + photon3_E ;
		pion1_pz = photon1_pz + photon3_pz;		
		invariantemasse2 = invariantemasse24;
		pi02_4vect = GetPhotons()->Particle(1)+GetPhotons()->Particle(3);
		pion2_E = photon2_E + photon4_E; 
		pion2_pz = photon2_pz + photon4_pz;
		Chi = Chi13;
	}
	
	if(Chi14 < Chi13 && Chi14 < Chi12){
		invariantemasse1 = invariantemasse14;
		pi01_4vect = GetPhotons()->Particle(0)+GetPhotons()->Particle(3);
		pion1_E = photon1_E + photon4_E ;
		pion1_pz = photon1_pz + photon4_pz;		
		invariantemasse2 = invariantemasse23;
		pi02_4vect = GetPhotons()->Particle(1)+GetPhotons()->Particle(2);
		pion2_E = photon2_E + photon3_E; 
		pion2_pz = photon2_pz + photon3_pz;
		Chi = Chi14;
	}	
// 
	pion1_m = invariantemasse1*invariantemasse1;
	pion2_m = invariantemasse2*invariantemasse2;

// 	printf("%.2f \t %.2f \t %.2f \t %.2f \t %.2f \n",Chi12,Chi13,Chi14,invariantemasse1,invariantemasse2);
// 	cout << Chi12 << "\t" << Chi13 << "\t" << Chi14 << invariantemasse1 << "\t" << invariantemasse2 << endl;

//___________________anti chi^2 test for pi0eta.. if this chi^2 is smaller than the chi^2 for 2pi0, then reject our event 
// 	Double_t Chi1[6];
// 				
// 	Chi1[0] = ChiPionEta(invariantemasse12,photon1_4vektor,photon2_4vektor,invariantemasse34,photon3_4vektor,photon3_4vektor);		
// 	Chi1[1] = ChiPionEta(invariantemasse34,photon3_4vektor,photon4_4vektor,invariantemasse12,photon1_4vektor,photon2_4vektor);
// 
// 	Chi1[2]= ChiPionEta(invariantemasse13,photon1_4vektor,photon3_4vektor,invariantemasse24,photon2_4vektor,photon4_4vektor);
// 	Chi1[3]= ChiPionEta(invariantemasse24,photon2_4vektor,photon4_4vektor,invariantemasse13,photon1_4vektor,photon3_4vektor);
// 
// 	Chi1[4]= ChiPionEta(invariantemasse14,photon1_4vektor,photon4_4vektor,invariantemasse23,photon2_4vektor,photon3_4vektor);
// 	Chi1[5]= ChiPionEta(invariantemasse23,photon2_4vektor,photon3_4vektor,invariantemasse14,photon1_4vektor,photon4_4vektor);
// 
// 	Double_t vergleich=Chi1[0];
// 
// 	for(Int_t i=1; i<6;i++){
// 		if(Chi1[i]<vergleich){vergleich = Chi1[i];}
// 	}

// 	cout << "FINALE: " << Chi << endl;		

		//TEST the target polarization (no entry with pt=0)
		test->Fill(pt);

	//Trigger setting
	for(Int_t i=0;i<GetTrigger()->GetNTriggerPattern();i++){	
// 	triggertest->Fill(GetTrigger()->GetTriggerPattern(i));
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



	for(Int_t j=0; j < GetTagger()->GetNTagged();j++)
	{
// 		if(Chi > vergleich){continue;} //reject the pi0eta events if chi^2 measured is used

		//reject the pbwo4
		if(GetPhotons()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;
		if(GetPhotons()->GetDetectors(1)==GTreeTrack::DETECTOR_PbWO4) continue;
		if(GetPhotons()->GetDetectors(2)==GTreeTrack::DETECTOR_PbWO4) continue;
		if(GetPhotons()->GetDetectors(3)==GTreeTrack::DETECTOR_PbWO4) continue;
		if(GetRootinos()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;

// 		if(pt==5){continue;} //for carbon measurements
// 		if(pt==5){continue;}
// 		if(GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j))==-1){continue;}//skip events with wrong polarisation!

		Double_t time= GetTagger()->GetTaggedTime(j) - 0.25*(GetPhotons()->GetTime(1)+GetPhotons()->GetTime(0)+GetPhotons()->GetTime(2)+GetPhotons()->GetTime(3));
		time1->Fill(time);

 		if(HistoManu::IsPrompt(time)){
			time_prompt->Fill(time);
		}

		if(HistoManu::IsRandom(time)){	
			time_side->Fill(time);
		}

		if(!(HistoManu::IsRandom(time) || HistoManu::IsPrompt(time))){continue;} //skip events which are not prompt or random

		poltable_energy->Fill(GetTagger()->GetTaggedEnergy(j));
		poltable_energy_weight->Fill(GetTagger()->GetTaggedEnergy(j),GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j)));

		//get target, beam and missng particle information
		TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
		TLorentzVector  beam_4vect = TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(j),GetTagger()->GetTaggedEnergy(j));
		TLorentzVector  missingp_4vect = beam_4vect + protonvektor_target - pi02_4vect-pi01_4vect;
		Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
		Double_t beamphoton1E=GetTagger()->GetTaggedEnergy(j);
		Double_t pb=GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j));
		//Missing Mass
		Double_t missingmass=missingp_4vect.M();


		//due to kinematic not possible -> kick them out
		if(anglethetaproton_rek > 90){continue;}
			
 		if(beamphoton1E<200 || beamphoton1E>1450){continue;}

		//energy dependent cuts
		Double_t oben_copl=180+15;
		Double_t unten_copl=180-15;
		Double_t oben_theta=12;
		Double_t unten_theta=-12;
		Double_t oben_mass=938+80;
		Double_t unten_mass=938-80;
		Double_t oben_mass_proton=938+80;
		Double_t unten_mass_proton=938-80;
			
		Double_t oben_inv=135+35;
		Double_t unten_inv=135-35;
		Double_t oben_inv_proton=135+35;
		Double_t unten_inv_proton=135-35;

		//Boost-System

		//Pi0Pi0-System
		TLorentzVector pi0pi0_4vektor=pi01_4vect+pi02_4vect;
		TLorentzVector pi0pi0_4vektor_boost = CMVector(pi0pi0_4vektor, beam_4vect, protonvektor_target);
		Double_t invariantmass_pi0pi0 = pi0pi0_4vektor_boost.M();
		TVector3 pi0pi0_3vektor_boost = pi0pi0_4vektor_boost.Vect();
		Double_t cospi0pi0 = pi0pi0_3vektor_boost.CosTheta();
		Double_t phimeson = TMath::RadToDeg()*(pi0pi0_4vektor.Vect().Phi());
			
		//pPi0-System for first Pion
		TLorentzVector ppi01_4vektor = missingp_4vect+pi01_4vect;
		TLorentzVector ppi01_4vektor_boost = CMVector(ppi01_4vektor, beam_4vect, protonvektor_target);
		Double_t invariantmass_ppi01 = ppi01_4vektor_boost.M();
		TVector3 ppi01_3vektor_boost = ppi01_4vektor_boost.Vect();
		Double_t cosppi01 = ppi01_3vektor_boost.CosTheta();
		Double_t phippi01 = TMath::RadToDeg()*ppi01_3vektor_boost.Phi();
			
		//pPi0-System for second Pion
		TLorentzVector ppi02_4vektor = missingp_4vect+pi02_4vect;
		TLorentzVector ppi02_4vektor_boost = CMVector(ppi02_4vektor, beam_4vect, protonvektor_target);
		Double_t invariantmass_ppi02 = ppi02_4vektor_boost.M();
		TVector3 ppi02_3vektor_boost = ppi02_4vektor_boost.Vect();
		Double_t cosppi02 = ppi02_3vektor_boost.CosTheta();
		Double_t phippi02 = TMath::RadToDeg()*ppi02_3vektor_boost.Phi();
			
			
		Double_t scalsignalpb_ppi0 = 0.5*pb;
		Double_t scalsignalpt_ppi0 = 0.5*pt;
			
		Double_t anglephiproton;
		Double_t phi;
		Double_t anglethetaproton_meas;
		Double_t thetadiff;

//_____________________________only one charged particle_________________________________

		if(GetRootinos()->GetNParticles()==1){

			//get charged particle information
			TLorentzVector proton_4vect_meas;
			proton_4vect_meas=GetRootinos()->Particle(0);

			//Phi-Difference
			anglephiproton = TMath::RadToDeg()*proton_4vect_meas.Vect().Phi();
			phi = phimeson - anglephiproton;
			if(phi < 0){phi = phi + 360;}
			
			//theta differenz
			anglethetaproton_meas = TMath::RadToDeg()*proton_4vect_meas.Vect().Theta();
			thetadiff=anglethetaproton_rek-anglethetaproton_meas;

			//due to kinematic not possible -> kick them out
			if(anglethetaproton_meas > 70){continue;}

			//...........................Histogramme vor jeglichen Cuts erstellen (3PED)..................................

			coplanarity_collerated->Fill(phi,time);
			thetaproton_collerated->Fill(thetadiff,time);
		
			missingmass_collerated_proton->Fill(missingmass,time);
		
			massesumme_collerated_proton->Fill(invariantemasse34,time);
			massesumme_collerated_proton->Fill(invariantemasse24,time);
			massesumme_collerated_proton->Fill(invariantemasse23,time);
			massesumme_collerated_proton->Fill(invariantemasse12,time);
			massesumme_collerated_proton->Fill(invariantemasse13,time);
			massesumme_collerated_proton->Fill(invariantemasse14,time);
			
			massegegenmasse_collerated_proton->Fill(invariantemasse12,invariantemasse34,time);
			massegegenmasse_collerated_proton->Fill(invariantemasse34,invariantemasse12,time);
			massegegenmasse_collerated_proton->Fill(invariantemasse24,invariantemasse13,time);
			massegegenmasse_collerated_proton->Fill(invariantemasse13,invariantemasse24,time);
			massegegenmasse_collerated_proton->Fill(invariantemasse14,invariantemasse23,time);
			massegegenmasse_collerated_proton->Fill(invariantemasse23,invariantemasse14,time);


			//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::single cuts:::::::::::::::::::::::::::::::::::::::::
					
			//__________________Polar angle______________________________________________________

			if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {
					
				coplanarity_theta_collerated->Fill(phi,time);
				missingmass_theta_collerated->Fill(missingmass,time);

				massesumme_theta_collerated->Fill(invariantemasse12,time);
				massesumme_theta_collerated->Fill(invariantemasse34,time);
				massesumme_theta_collerated->Fill(invariantemasse13,time);
				massesumme_theta_collerated->Fill(invariantemasse24,time);
				massesumme_theta_collerated->Fill(invariantemasse14,time);
				massesumme_theta_collerated->Fill(invariantemasse23,time);
					

				massegegenmasse_theta_collerated->Fill(invariantemasse12,invariantemasse34,time);
				massegegenmasse_theta_collerated->Fill(invariantemasse34,invariantemasse12,time);
				massegegenmasse_theta_collerated->Fill(invariantemasse13,invariantemasse24,time);
				massegegenmasse_theta_collerated->Fill(invariantemasse24,invariantemasse13,time);
				massegegenmasse_theta_collerated->Fill(invariantemasse23,invariantemasse14,time);
				massegegenmasse_theta_collerated->Fill(invariantemasse14,invariantemasse23,time);
	
			}

			//__________________Coplanarity______________________________________________________
					
			if((phi > ((unten_copl)) && phi < ((oben_copl)))){

				missingmass_copl_collerated->Fill(missingmass,time);

				massesumme_copl_collerated->Fill(invariantemasse12,time);
				massesumme_copl_collerated->Fill(invariantemasse34,time);
				massesumme_copl_collerated->Fill(invariantemasse13,time);
				massesumme_copl_collerated->Fill(invariantemasse24,time);
				massesumme_copl_collerated->Fill(invariantemasse14,time);
				massesumme_copl_collerated->Fill(invariantemasse23,time);

				massegegenmasse_copl_collerated->Fill(invariantemasse12,invariantemasse34,time);
				massegegenmasse_copl_collerated->Fill(invariantemasse34,invariantemasse12,time);
				massegegenmasse_copl_collerated->Fill(invariantemasse13,invariantemasse24,time);
				massegegenmasse_copl_collerated->Fill(invariantemasse24,invariantemasse13,time);
				massegegenmasse_copl_collerated->Fill(invariantemasse23,invariantemasse14,time);
				massegegenmasse_copl_collerated->Fill(invariantemasse14,invariantemasse23,time);
					
						
				if(GetRootinos()->HasTAPS(0)){
					thetaproton_copl_taps_collerated->Fill(thetadiff,time);
				}
				else
				{
					thetaproton_copl_collerated->Fill(thetadiff,time);
				}
					
			}
										
			//__________________invariant mass______________________________________________________
					
			if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){

				coplanarity_inv_collerated->Fill(phi,time);
				missingmass_inv_collerated_proton->Fill(missingmass,time);
					
				if(GetRootinos()->HasTAPS(0)){
					thetaproton_inv_taps_collerated->Fill(thetadiff,time);
				}
				else
				{
					thetaproton_inv_collerated->Fill(thetadiff,time);
				}

			}

			//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::aufeinanderfolgende Cuts:::::::::::::::::::::::::::::::::::::::::
					

//--------------------------------------------------------CopLanaritY-Cuts----------------------
					
		// 	__________________Missing Mass______________________________________________________
					
			if(missingmass>(unten_mass) && missingmass<(oben_mass)){
					
				coplanarity_mass_collerated->Fill(phi,time);


		// 	__________________Polarwinkel______________________________________________________
					
				if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {

					coplanarity_mass_theta_collerated->Fill(phi,time);
					
						
		// 			__________________invariante Masse______________________________________________________
					if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
							//if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
						coplanarity_mass_theta_inv_collerated->Fill(phi,time);
						coplanarityverteilung_collerated->Fill(phi,beamphoton1E,time);
									
					}
				}
			}
										
// --------------------------------------------------------Missing Mass Cuts----------------------
															
		// 	__________________Polarwinkel______________________________________________________
			
			if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {
					
				missingmass_theta_collerated->Fill(missingmass,time);

		// 	__________________Coplanarity______________________________________________________
					
				if((phi > ( (unten_copl)) && phi < ((oben_copl)))){

					missingmass_theta_copl_collerated->Fill(missingmass,time);
					
					
		// 	__________________invariante Masse______________________________________________________
					if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
							//if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
						missingmass_theta_copl_inv_collerated->Fill(missingmass,time);
						missingmassverteilung_collerated_proton->Fill(missingmass,beamphoton1E,time);
					}
				}
			}
					
// --------------------------------------------------------PolarWinkel Cuts----------------------			
					
		// 	__________________Missing Mass______________________________________________________
					
			if(missingmass>(unten_mass) && missingmass<(oben_mass)){
					
				if(GetRootinos()->HasTAPS(0)){
					thetaproton_mass_taps_collerated->Fill(thetadiff,time);
				}
				else{
					thetaproton_mass_collerated->Fill(thetadiff,time);		
				}		
					
		// 	__________________Coplanarity______________________________________________________
					
				if((phi > ((unten_copl)) && phi < ((oben_copl)))){

					if(GetRootinos()->HasTAPS(0)){
						thetaproton_mass_copl_taps_collerated->Fill(thetadiff,time);
					}
					else{
						thetaproton_mass_copl_collerated->Fill(thetadiff,time);		
					}	
					
					
		// __________________invariante Masse______________________________________________________
					if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
// 							if(((invariantemasse1 > (pionmasse + unten_inv) && invariantemasse1 < (pionmasse + oben_inv)) && (invariantemasse2 > (pionmasse + unten_inv) && invariantemasse2 < (pionmasse + oben_inv)))){
					
						if(GetRootinos()->HasTAPS(0)){
							thetaproton_mass_copl_inv_taps_collerated->Fill(thetadiff,time);
							thetaverteilung_collerated_taps->Fill(thetadiff,beamphoton1E,time);
						}
						else{
							thetaproton_mass_copl_inv_collerated->Fill(thetadiff,time);	
							thetaverteilung_collerated->Fill(thetadiff,beamphoton1E,time);	
						}
					}
				}	
			} 
// ................................invarianteMasse Cut...................................................................................
										
		// 	__________________Missing Mass______________________________________________________
					
			if(missingmass>(unten_mass) && missingmass<(oben_mass)){

				massesumme_mass_collerated_proton->Fill(invariantemasse12,time);
				massesumme_mass_collerated_proton->Fill(invariantemasse34,time);
				massesumme_mass_collerated_proton->Fill(invariantemasse13,time);
				massesumme_mass_collerated_proton->Fill(invariantemasse24,time);
				massesumme_mass_collerated_proton->Fill(invariantemasse14,time);
				massesumme_mass_collerated_proton->Fill(invariantemasse23,time);

				massegegenmasse_mass_collerated_proton->Fill(invariantemasse12,invariantemasse34,time);
				massegegenmasse_mass_collerated_proton->Fill(invariantemasse34,invariantemasse12,time);
				massegegenmasse_mass_collerated_proton->Fill(invariantemasse13,invariantemasse24,time);
				massegegenmasse_mass_collerated_proton->Fill(invariantemasse24,invariantemasse13,time);
				massegegenmasse_mass_collerated_proton->Fill(invariantemasse23,invariantemasse14,time);
				massegegenmasse_mass_collerated_proton->Fill(invariantemasse14,invariantemasse23,time);
									
					
		// 	__________________Polarwinkel______________________________________________________
					
				if((phi > ((unten_copl)) && phi < (oben_copl))){

					massesumme_mass_copl_collerated->Fill(invariantemasse12,time);
					massesumme_mass_copl_collerated->Fill(invariantemasse34,time);
					massesumme_mass_copl_collerated->Fill(invariantemasse13,time);
					massesumme_mass_copl_collerated->Fill(invariantemasse24,time);
					massesumme_mass_copl_collerated->Fill(invariantemasse14,time);
					massesumme_mass_copl_collerated->Fill(invariantemasse23,time);

					massegegenmasse_mass_copl_collerated->Fill(invariantemasse12,invariantemasse34,time);
					massegegenmasse_mass_copl_collerated->Fill(invariantemasse34,invariantemasse12,time);
					massegegenmasse_mass_copl_collerated->Fill(invariantemasse13,invariantemasse24,time);
					massegegenmasse_mass_copl_collerated->Fill(invariantemasse24,invariantemasse13,time);
					massegegenmasse_mass_copl_collerated->Fill(invariantemasse23,invariantemasse14,time);
					massegegenmasse_mass_copl_collerated->Fill(invariantemasse14,invariantemasse23,time);		
						
		// 	__________________Coplanarity______________________________________________________
																							
					if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {

						massesumme_mass_theta_copl_collerated->Fill(invariantemasse12,time);
						massesumme_mass_theta_copl_collerated->Fill(invariantemasse34,time);
						massesumme_mass_theta_copl_collerated->Fill(invariantemasse13,time);
						massesumme_mass_theta_copl_collerated->Fill(invariantemasse24,time);
						massesumme_mass_theta_copl_collerated->Fill(invariantemasse14,time);
						massesumme_mass_theta_copl_collerated->Fill(invariantemasse23,time);

						massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse12,invariantemasse34,time);
						massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse34,invariantemasse12,time);
						massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse13,invariantemasse24,time);
						massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse24,invariantemasse13,time);
						massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse23,invariantemasse14,time);
						massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse14,invariantemasse23,time);
					
					}							    				
				}		
			}
					//stuff after all cuts and creation of needed analysis histograms
	
			if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && (((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < oben_inv) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv)))){

// 			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& (((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < oben_inv) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv)))){

				//set the different fitter
				fitter.Set(GetPhotons(), 0, 1,2,3);
				fitter.SetProton(GetRootinos(),0);
				fitter.SetBeam(beamphoton1E);	
	
				fitter1.Set(GetPhotons(), 0, 2,1,3);
				fitter1.SetProton(GetRootinos(),0);
				fitter1.SetBeam(beamphoton1E);

				fitter2.Set(GetPhotons(), 0, 3,1,2);
				fitter2.SetProton(GetRootinos(),0);
				fitter2.SetBeam(beamphoton1E);

				//solve the different fitter
				Bool_t fit1,fit2,fit3;

				fit1=fitter.Solve();
				fit2=fitter1.Solve();
				fit3=fitter2.Solve();

				Double_t cl=0;
				Double_t    cl1 = fitter.GetFittedResult().Probability;
 				Double_t    cl2 = fitter1.GetFittedResult().Probability;
				Double_t    cl3 = fitter2.GetFittedResult().Probability;

				if(!(fit1||fit2||fit3))continue; //if no fitter solved the problem!


				//Find the best result

				//only one fitter solved
				if(fit1 && !(fit2) && !(fit3))cl1=cl;
				if(fit2 && !(fit1) && !(fit3))cl2=cl;
				if(fit3 && !(fit1) && !(fit2))cl3=cl;

				//two fitter solved

				//1 and 2
				if(fit1 && (fit2) && !(fit3) && (cl1>cl2))cl=cl1;
				if(fit1 && (fit2) && !(fit3) && (cl2>cl1))cl=cl2;

				//1 and 3
				if(fit1 && (!fit2) && (fit3) && (cl1>cl3))cl=cl1;
				if(fit1 && (!fit2) && (fit3) && (cl3>cl1))cl=cl3;

				//2 and 3
				if(!(fit1) && (fit2) && (fit3) && (cl2>cl3))cl=cl2;
				if(!(fit1) && (fit2) && (fit3) && (cl3>cl2))cl=cl3;

				//three fitter solved
				if((fit1) && (fit2) && (fit3) && (cl1>cl2 && cl1>cl3))cl=cl1;
				if((fit1) && (fit2) && (fit3) && (cl2>cl1 && cl2>cl3))cl=cl2;
				if((fit1) && (fit2) && (fit3) && (cl3>cl2 &&cl3>cl1))cl=cl3;

				//define the invariant masses and the subparticles
				Int_t pi01_sub1, pi01_sub2, pi02_sub1,pi02_sub2;
				if(cl==cl1){
					pi01_sub1=0;
					pi01_sub2=1;
					pi02_sub1=2;
					pi02_sub2=3;
					invariantemasse1=invariantemasse12;
					invariantemasse2=invariantemasse34;
				}

				if(cl==cl2){
					pi01_sub1=0;
					pi01_sub2=2;
					pi02_sub1=1;
					pi02_sub2=3;
					invariantemasse1=invariantemasse13;
					invariantemasse2=invariantemasse24;
				}
				if(cl==cl3){
					pi01_sub1=0;
					pi01_sub2=3;
					pi02_sub1=1;
					pi02_sub2=2;
					invariantemasse1=invariantemasse14;
					invariantemasse2=invariantemasse23;
				}
// 				cout << "Eventnumber: " <<GetEventParameters()->GetEventNumber() << "\t" << "Ntagged: " << GetTagger()->GetNTagged() << "\t" << j <<"\t" << "FitSolve: " << fit1 << "\t" << fit2 << "\t" << fit3  <<"\t" << "CL: " << cl1 << "\t" << cl2 << "\t" << cl3 << "\t" << "BEST CL: " <<cl <<  endl;
				CL->Fill(cl,time);

				//PULLS
				if(cl<0.01 || cl>1.01)continue;

				for(int i=0; i<4; i++)
				{

					if(!GetPhotons()->HasTAPS(i))continue; //reject if they are not detected in TAPS

					for(int p=0; p<3; p++)
					{
						std::stringstream s;
						s << "Ph" << i;
						s << "[" << p << "]";

						if(cl==cl1)pulls_4g_TAPS->Fill(fitter.GetFittedResult().Variables.at(s.str()).Pull, p, time);
						if(cl==cl2)pulls_4g_TAPS->Fill(fitter1.GetFittedResult().Variables.at(s.str()).Pull, p, time);
						if(cl==cl3)pulls_4g_TAPS->Fill(fitter2.GetFittedResult().Variables.at(s.str()).Pull, p, time);


					}
				}

				for(int i=0; i<4; i++)
				{

					if(!GetPhotons()->HasCB(i))continue;// reject if they are not detected in CB

					for(int p=0; p<3; p++)
					{
						std::stringstream s;
						s << "Ph" << i;
						s << "[" << p << "]";
// 
						if(cl==cl1)pulls_4g_CB->Fill(fitter.GetFittedResult().Variables.at(s.str()).Pull, p, time);
						if(cl==cl2)pulls_4g_CB->Fill(fitter1.GetFittedResult().Variables.at(s.str()).Pull, p, time);
						if(cl==cl3)pulls_4g_CB->Fill(fitter2.GetFittedResult().Variables.at(s.str()).Pull, p, time);

					}
				}


				for(Int_t p=0;p<3;p++){

					std::stringstream s;
					s << "Be[" << p << "]";;

					if(cl==cl1)pulls_beam->Fill(fitter.GetFittedResult().Variables.at(s.str()).Pull, p, time);
					if(cl==cl2)pulls_beam->Fill(fitter1.GetFittedResult().Variables.at(s.str()).Pull, p, time);
					if(cl==cl3)pulls_beam->Fill(fitter2.GetFittedResult().Variables.at(s.str()).Pull, p, time);
				}

				for(Int_t p=0;p<2;p++){

					std::stringstream s;
					s << "PA[" << p << "]";

					if(cl==cl1&& GetRootinos()->HasTAPS(0))pulls_proton_TAPS->Fill(fitter.GetFittedResult().Variables.at(s.str()).Pull,p,time);
					if(cl==cl2&& GetRootinos()->HasTAPS(0))pulls_proton_TAPS->Fill(fitter1.GetFittedResult().Variables.at(s.str()).Pull,p,time);
					if(cl==cl3&& GetRootinos()->HasTAPS(0))pulls_proton_TAPS->Fill(fitter2.GetFittedResult().Variables.at(s.str()).Pull,p,time);

					if(cl==cl1&& GetRootinos()->HasCB(0))pulls_proton_CB->Fill(fitter.GetFittedResult().Variables.at(s.str()).Pull,p,time);
					if(cl==cl2&& GetRootinos()->HasCB(0))pulls_proton_CB->Fill(fitter1.GetFittedResult().Variables.at(s.str()).Pull,p,time);
					if(cl==cl3&& GetRootinos()->HasCB(0))pulls_proton_CB->Fill(fitter2.GetFittedResult().Variables.at(s.str()).Pull,p,time);

				}

				//Check of fit quality
				Double_t anglethetaproton_fit;
				Double_t anglephiproton_fit; 
				Double_t theta_fit,phi_fit;
				if(thetadiff > (unten_theta) && thetadiff < (oben_theta)){

					//Phi-Difference
					if(cl==cl1)anglephiproton_fit= TMath::RadToDeg()*fitter.GetFittedResult().Variables.at("PA[1]").Value.After;
					if(cl==cl2)anglephiproton_fit= TMath::RadToDeg()*fitter1.GetFittedResult().Variables.at("PA[1]").Value.After;
					if(cl==cl3)anglephiproton_fit= TMath::RadToDeg()*fitter2.GetFittedResult().Variables.at("PA[1]").Value.After;

					phi_fit = phimeson - anglephiproton_fit;
					if(phi_fit < 0){phi_fit = phi_fit + 360;}
					coplanarity_fit_vs_rek->Fill(phi_fit,time);
				}

				if(phi > (unten_copl) && phi < (oben_copl)){

					//theta differenz
					if(cl==cl1)anglethetaproton_fit = TMath::RadToDeg()*fitter.GetFittedResult().Variables.at("PA[0]").Value.After;
					if(cl==cl2)anglethetaproton_fit = TMath::RadToDeg()*fitter1.GetFittedResult().Variables.at("PA[0]").Value.After;
					if(cl==cl3)anglethetaproton_fit = TMath::RadToDeg()*fitter2.GetFittedResult().Variables.at("PA[0]").Value.After;

					theta_fit = anglethetaproton_meas-anglethetaproton_fit;
					if(GetRootinos()->HasCB(0))thetaproton_fit_vs_rek_CB->Fill(theta_fit,time);
					if(GetRootinos()->HasTAPS(0))thetaproton_fit_vs_rek_TAPS->Fill(theta_fit,time);

				}

				if(cl==cl1)invmass_diff_rec_kin->Fill(invariantemasse1-(fitter.GetFittedGamma(0)+fitter.GetFittedGamma(1)).M(),time);
				if(cl==cl1)invmass_diff_rec_kin->Fill(invariantemasse2-(fitter.GetFittedGamma(2)+fitter.GetFittedGamma(3)).M(),time);

				if(cl==cl2)invmass_diff_rec_kin->Fill(invariantemasse1-(fitter1.GetFittedGamma(0)+fitter1.GetFittedGamma(1)).M(),time);
				if(cl==cl2)invmass_diff_rec_kin->Fill(invariantemasse2-(fitter1.GetFittedGamma(2)+fitter1.GetFittedGamma(3)).M(),time);

				if(cl==cl3)invmass_diff_rec_kin->Fill(invariantemasse1-(fitter2.GetFittedGamma(0)+fitter2.GetFittedGamma(1)).M(),time);
				if(cl==cl3)invmass_diff_rec_kin->Fill(invariantemasse2-(fitter2.GetFittedGamma(2)+fitter2.GetFittedGamma(3)).M(),time);
// 
//    				if(HistoManu::IsPrompt(time)) calfile << anglethetaproton_meas << "\t" << GetEventParameters()->GetEventNumber() << "\t" <<GetRootinos()->GetDetectors(0)<< "\n";
// 
//      			if(HistoManu::IsRandom(time)) calfile1 << anglethetaproton_meas<< "\t" << GetEventParameters()->GetEventNumber() << "\t" << GetRootinos()->GetDetectors(0) <<"\n";
// 
				Double_t openinganle_pi01=TMath::ACos(TMath::Cos(GetPhotons()->Particle(pi01_sub1).Phi()-(GetPhotons()->Particle(pi01_sub2)).Phi())*TMath::Sin(GetPhotons()->Particle(pi01_sub1).Theta())*TMath::Sin((GetPhotons()->Particle(pi01_sub2)).Theta())+ TMath::Cos(GetPhotons()->Particle(pi01_sub1).Theta())*TMath::Cos((GetPhotons()->Particle(pi01_sub2)).Theta()) )*TMath::RadToDeg();

				Double_t openinganle_pi02=TMath::ACos(TMath::Cos(GetPhotons()->Particle(pi02_sub1).Phi()-(GetPhotons()->Particle(pi02_sub2)).Phi())*TMath::Sin(GetPhotons()->Particle(pi02_sub1).Theta())*TMath::Sin((GetPhotons()->Particle(pi02_sub2)).Theta())+ TMath::Cos(GetPhotons()->Particle(pi02_sub1).Theta())*TMath::Cos((GetPhotons()->Particle(pi02_sub2)).Theta()))*TMath::RadToDeg();

				openangle_pi01_to_pi02_energy->Fill(openinganle_pi01,openinganle_pi02, beamphoton1E,time);

				massegegenmasse_allcuts_collerated->Fill(invariantemasse12,invariantemasse34,time);
				massegegenmasse_allcuts_collerated->Fill(invariantemasse34,invariantemasse12,time);
				massegegenmasse_allcuts_collerated->Fill(invariantemasse13,invariantemasse24,time);
				massegegenmasse_allcuts_collerated->Fill(invariantemasse24,invariantemasse13,time);
				massegegenmasse_allcuts_collerated->Fill(invariantemasse23,invariantemasse14,time);
				massegegenmasse_allcuts_collerated->Fill(invariantemasse14,invariantemasse23,time);
// 
				massegegenmasse_allcutsandinv_collerated->Fill(invariantemasse2,invariantemasse1,time);
				massegegenmasse_allcutsandinv_collerated->Fill(invariantemasse1,invariantemasse2,time);
			


				//Some Cout Checks.

// 				if(cl==cl1)cout << "Start BEAME: "<< beamphoton1E << "\t" << "FIT BEAME: " << fitter.GetFittedResult().Variables.at("Be[0]").Value.After << "\t" << "THETA: " << fitter.GetFittedResult().Variables.at("Be[1]").Value.After << "\t" << "PHI: " << fitter.GetFittedResult().Variables.at("Be[2]").Value.After << endl;
// 
// 				if(cl==cl2)cout << "Start BEAME: "<< beamphoton1E << "\t" << "FIT BEAME: " << fitter1.GetFittedResult().Variables.at("Be[0]").Value.After << "\t" << "THETA: " << fitter1.GetFittedResult().Variables.at("Be[1]").Value.After << "\t" << "PHI: " << fitter1.GetFittedResult().Variables.at("Be[2]").Value.After << endl;
// 
// 				if(cl==cl3)cout << "Start BEAME: "<< beamphoton1E << "\t" << "FIT BEAME: " << fitter2.GetFittedResult().Variables.at("Be[0]").Value.After << "\t" << "THETA: " << fitter2.GetFittedResult().Variables.at("Be[1]").Value.After << "\t" << "PHI: " << fitter2.GetFittedResult().Variables.at("Be[2]").Value.After << endl;

// 				cout << "\t" << GetRootinos()->GetNParticles() << "\t" << fitter.Solve() << "\t" << anglethetaproton_fit << "\t" << anglethetaproton_meas << "\t" << TMath::RadToDeg()*fitter.GetFittedResult().Variables.at("PA[0]").Value.Before << endl;
// 				cout << beamphoton1E << "\t" << fitter.GetFittedResult().Variables.at("Be[0]").Value.After << "\t" << fitter.GetFittedResult().Variables.at("Be[0]").Value.Before << "\t" << cl << endl;
/*
				cout << phi_fit << "\t" << anglephiproton << "\t" << anglephiproton_fit << endl;
				cout << "m12: " << (fitter.GetFittedGamma(0)+fitter.GetFittedGamma(1)).M() << "\t"<< invariantemasse12<< endl;
				cout << "m34:  " <<(fitter.GetFittedGamma(2)+fitter.GetFittedGamma(3)).M() << "\t" << invariantemasse34 << endl;
				cout << "m13: " << (fitter1.GetFittedGamma(0)+fitter1.GetFittedGamma(1)).M() << "\t" << invariantemasse13 << endl;
				cout << "m24: " << (fitter1.GetFittedGamma(2)+fitter1.GetFittedGamma(3)).M() << "\t" << invariantemasse24 << endl;
				cout << "m14: " << (fitter2.GetFittedGamma(0)+fitter2.GetFittedGamma(1)).M() << "\t" << invariantemasse14 <<  endl;
				cout << "m23: " << (fitter2.GetFittedGamma(2)+fitter2.GetFittedGamma(3)).M() << "\t" << invariantemasse23 << endl;
				cout << "_________________________________________________________________________________________" << endl;*/

				
				if(missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)  && thetadiff > (unten_theta) && thetadiff < (oben_theta) && phi > (unten_copl) && phi < (oben_copl) && invariantemasse2>unten_inv && invariantemasse2<oben_inv){
	
					massesumme_mass_theta_copl_inv_beam_collerated->Fill(invariantemasse1,beamphoton1E,time);

					if(invariantemasse1>unten_inv && invariantemasse1<oben_inv){
					
						if(GetRootinos()->HasCB(0)){
							Check_CBdE_E->Fill(proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
						}
		
						if(GetRootinos()->HasTAPS(0)){
							Check_TAPSdE_E->Fill(proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
						}
		
						if(GetScalers()->GetNEntries()==0){
							h_energy_sum_2pi0_5ped->Fill(GetTrigger()->GetEnergySum(),time);
							h_energy_sum_2pi0_5ped_weightcos2pi0->Fill(GetTrigger()->GetEnergySum(),time,dwq(beamphoton1E, cospi0pi0, "2pi0",""));
							h_energy_sum_2pi0_5ped_weightcosppi0->Fill(GetTrigger()->GetEnergySum(),time,dwq(beamphoton1E,  0.5*(cosppi01+cosppi02), "2pi0","cosppi0"));
			
							cos2pi0_beamphoton_rek->Fill(cospi0pi0,beamphoton1E,time);
							m2pi0_beamphoton_energysum_rek->Fill(invariantmass_pi0pi0,beamphoton1E,GetTrigger()->GetEnergySum(),time);
							cosppi0_beamphoton_energysum_rek->Fill(cosppi01,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);
							cosppi0_beamphoton_energysum_rek->Fill(cosppi02,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);
							mppi0_beamphoton_energysum_rek->Fill(invariantmass_ppi01,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);
							mppi0_beamphoton_energysum_rek->Fill(invariantmass_ppi02,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);
		
							for(Int_t i=0; i<GetTracks()->GetNTracks() ;i++){//loop over all GetTracks()
								if(GetTracks()->HasTAPS(i)){//if track is from taps
									for(Int_t j=0; j<GetDetectorHits()->GetNBaF2Hits() ;j++){//loop over all hits
										if(GetDetectorHits()->GetBaF2Hits(j)!=GetTracks()->GetCentralCrystal(i)){continue;}//if central crystal index is unequal baf2 crystal index skip
											h_energy_taps_2pi0_5ped->Fill(GetDetectorHits()->GetBaF2Energy(j),GetTrigger()->GetEnergySum(),time);
									}
							
								}
							}
						}
						cosverteilung_collerated_proton->Fill(cospi0pi0,beamphoton1E,time);
	
					}

				}
	
/*
				invverteilung_collerated_proton->Fill(invariantmass_pi0pi0,beamphoton1E,time);
				cosppi0verteilung_collerated_proton->Fill(cosppi01,beamphoton1E,time,0.5);
				cosppi0verteilung_collerated_proton->Fill(cosppi02,beamphoton1E,time,0.5);
				invppi0verteilung_collerated_proton->Fill(invariantmass_ppi01,beamphoton1E,time,0.5);
				invppi0verteilung_collerated_proton->Fill(invariantmass_ppi02,beamphoton1E,time,0.5);
				cosverteilung_thetaproton_collerated_proton->Fill(cospi0pi0,anglethetaproton_meas, time);;					
				if(pt>0 || pt==5){
				
 					if(planesetting=="PARA"){

						kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time);
						kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
						kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						kristallminus_targetplus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						kristallminus_targetplus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						kristallminus_targetplus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						kristallminus_targetplus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
						kristallminus_targetplus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetplus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallminus_targetplus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
						kristallminus_targetplus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetplus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						kristallminus_targetplus_im_ppi0_collerated_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						kristallminus_targetplus_im_ppi0_collerated_pb_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetplus_im_ppi0_collerated_pt_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallminus_targetplus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						kristallminus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

					}

					if(planesetting=="PERP"){

						kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time);
						kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
						kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						kristallplus_targetplus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						kristallplus_targetplus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						kristallplus_targetplus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						kristallplus_targetplus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
						kristallplus_targetplus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetplus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallplus_targetplus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
						kristallplus_targetplus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetplus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						kristallplus_targetplus_im_ppi0_collerated_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						kristallplus_targetplus_im_ppi0_collerated_pb_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetplus_im_ppi0_collerated_pt_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallplus_targetplus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						kristallplus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);
					}
				}
		
				if(pt<0 || pt==5){

 					if(planesetting=="PARA"){

						kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time);
						kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
						kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						kristallminus_targetminus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						kristallminus_targetminus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						kristallminus_targetminus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						kristallminus_targetminus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
						kristallminus_targetminus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetminus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallminus_targetminus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
						kristallminus_targetminus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetminus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						kristallminus_targetminus_im_ppi0_collerated_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						kristallminus_targetminus_im_ppi0_collerated_pb_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetminus_im_ppi0_collerated_pt_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallminus_targetminus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						kristallminus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallminus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

					}

 					if(planesetting=="PERP"){

						kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time);
						kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
						kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						kristallplus_targetminus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						kristallplus_targetminus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						kristallplus_targetminus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						kristallplus_targetminus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
						kristallplus_targetminus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetminus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallplus_targetminus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
						kristallplus_targetminus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetminus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						kristallplus_targetminus_im_ppi0_collerated_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						kristallplus_targetminus_im_ppi0_collerated_pb_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetminus_im_ppi0_collerated_pt_proton->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						kristallplus_targetminus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						kristallplus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						kristallplus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

					}
				}*/
			}



// 			}

		}
//_________________________________________________ end of all events with proton identified!!!!!_________

	
/*
//...........................Histogramme vor jeglichen Cuts erstellen(2+3ped)..................................
		missingmass_collerated->Fill(missingmass,time);

		massesumme_collerated->Fill(invariantemasse34,time);
		massesumme_collerated->Fill(invariantemasse24,time);
		massesumme_collerated->Fill(invariantemasse23,time);
		massesumme_collerated->Fill(invariantemasse12,time);
		massesumme_collerated->Fill(invariantemasse13,time);
		massesumme_collerated->Fill(invariantemasse14,time);

		massegegenmasse_collerated->Fill(invariantemasse12,invariantemasse34,time);
		massegegenmasse_collerated->Fill(invariantemasse34,invariantemasse12,time);
		massegegenmasse_collerated->Fill(invariantemasse24,invariantemasse13,time);
		massegegenmasse_collerated->Fill(invariantemasse13,invariantemasse24,time);
		massegegenmasse_collerated->Fill(invariantemasse14,invariantemasse23,time);
		massegegenmasse_collerated->Fill(invariantemasse23,invariantemasse14,time);

		if(GetRootinos()->GetNParticles()==1){

				if(!(thetadiff > (unten_theta) && thetadiff < (oben_theta)) ||  !(phi > (unten_copl) && phi < (oben_copl))){continue;} 

		}					
			
// 			::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::einzelne Cuts:::::::::::::::::::::::::::::::::::::::::
			
			
// 			__________________invariante Masse______________________________________________________
		if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < oben_inv) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
			//if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){

			missingmass_inv_collerated->Fill(missingmass,time);
			missingmassverteilung_collerated->Fill(missingmass,beamphoton1E,time);
		}

		if(missingmass>(unten_mass) && missingmass<(oben_mass)){
									
			massesumme_mass_collerated->Fill(invariantemasse34,time);
			massesumme_mass_collerated->Fill(invariantemasse24,time);
			massesumme_mass_collerated->Fill(invariantemasse23,time);
			massesumme_mass_collerated->Fill(invariantemasse12,time);
			massesumme_mass_collerated->Fill(invariantemasse13,time);
			massesumme_mass_collerated->Fill(invariantemasse14,time);
	
			massegegenmasse_mass_collerated->Fill(invariantemasse12,invariantemasse34,time);
			massegegenmasse_mass_collerated->Fill(invariantemasse34,invariantemasse12,time);
			massegegenmasse_mass_collerated->Fill(invariantemasse24,invariantemasse13,time);
			massegegenmasse_mass_collerated->Fill(invariantemasse13,invariantemasse24,time);
			massegegenmasse_mass_collerated->Fill(invariantemasse14,invariantemasse23,time);
			massegegenmasse_mass_collerated->Fill(invariantemasse23,invariantemasse14,time);
		}
					
						

		if((missingmass>(unten_mass) && missingmass<(oben_mass))&& (((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < oben_inv) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv)))&&GetTrigger()->GetNErrors()!=0){
			events_witherror->Fill(phimeson,cospi0pi0, beamphoton1E);
		}

		if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)))&& ((invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
			events_all->Fill(phimeson,cospi0pi0, beamphoton1E);
		}

		if((missingmass>(unten_mass) && missingmass<(oben_mass))&& (((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < oben_inv) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv)))){

			cosverteilung_collerated->Fill(cospi0pi0,beamphoton1E,time);
			invverteilung_collerated->Fill(invariantmass_pi0pi0,beamphoton1E,time);
			cosppi0verteilung_collerated->Fill(cosppi01,beamphoton1E,time,0.5);
			cosppi0verteilung_collerated->Fill(cosppi02,beamphoton1E,time,0.5);
			invppi0verteilung_collerated->Fill(invariantmass_ppi01,beamphoton1E,time,0.5);
			invppi0verteilung_collerated->Fill(invariantmass_ppi02,beamphoton1E,time,0.5);
			
 			if(pt>0 || pt==5){
				
 				if(planesetting=="PARA"){

					kristallminus_targetplus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E,time);
					kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
					kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
					kristallminus_targetplus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
					kristallminus_targetplus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
					kristallminus_targetplus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

					kristallminus_targetplus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
					kristallminus_targetplus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetplus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallminus_targetplus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
					kristallminus_targetplus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetplus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

					kristallminus_targetplus_im_ppi0_collerated->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
					kristallminus_targetplus_im_ppi0_collerated_pb->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetplus_im_ppi0_collerated_pt->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallminus_targetplus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
					kristallminus_targetplus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetplus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

				}

				if(planesetting=="PERP"){

					kristallplus_targetplus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E,time);
					kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
					kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
					kristallplus_targetplus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
					kristallplus_targetplus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
					kristallplus_targetplus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

					kristallplus_targetplus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
					kristallplus_targetplus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetplus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallplus_targetplus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
					kristallplus_targetplus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetplus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

					kristallplus_targetplus_im_ppi0_collerated->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
					kristallplus_targetplus_im_ppi0_collerated_pb->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetplus_im_ppi0_collerated_pt->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallplus_targetplus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
					kristallplus_targetplus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetplus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);
				}
			}
		
			if(pt<0 || pt==5){

				if(planesetting=="PARA"){

					kristallminus_targetminus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E,time);
					kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
					kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
					kristallminus_targetminus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
					kristallminus_targetminus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
					kristallminus_targetminus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

					kristallminus_targetminus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
					kristallminus_targetminus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetminus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallminus_targetminus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
					kristallminus_targetminus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetminus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

					kristallminus_targetminus_im_ppi0_collerated->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
					kristallminus_targetminus_im_ppi0_collerated_pb->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetminus_im_ppi0_collerated_pt->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallminus_targetminus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
					kristallminus_targetminus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallminus_targetminus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

				}

 				if(planesetting=="PERP"){

					kristallplus_targetminus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E,time);
					kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,time,pb);
					kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,time,pt);
	
					kristallplus_targetminus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time);
					kristallplus_targetminus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
					kristallplus_targetminus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

					kristallplus_targetminus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,time,0.5);
					kristallplus_targetminus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetminus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallplus_targetminus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,time,0.5);
					kristallplus_targetminus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetminus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

					kristallplus_targetminus_im_ppi0_collerated->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
					kristallplus_targetminus_im_ppi0_collerated_pb->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetminus_im_ppi0_collerated_pt->Fill(invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
					kristallplus_targetminus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
					kristallplus_targetminus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
					kristallplus_targetminus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

				}
			}
					
		}//end of analysis histos*/
	}//GetTagger()
}//4 gammas
}


void P2Pi0Analyse::fOnEndProcessing() {

}

Bool_t	P2Pi0Analyse::Write(){
//FOR CARBON
// polsetting=P2Pi0Analyse::poledge(inputFile);
// TString *data=new TString(polsetting);
// TDirectory* curDir1  = outputFile->mkdir(Form("%s",data->Data()));

TDirectory* curDir1  = outputFile->mkdir(Form("Butanol%i",(Int_t)GetLinpol()->GetEdgeSetting()));


TString* filename1 = new TString(inputFile->GetPath());
TString* path11 = new TString(filename1->Tokenize("_")->At(filename1->Tokenize("_")->GetEntries()-1)->GetName());
path11->Resize(path11->Length()-7);
TNamed filenumber=TNamed(path11->Data(), "Filenumber");

h_energy_sum_2pi0_5ped->Write();
h_energy_sum_2pi0_5ped_weightcos2pi0->Write();
h_energy_sum_2pi0_5ped_weightcosppi0->Write();
h_energy_taps_2pi0_5ped->Write();
test->Write();
cos2pi0_beamphoton_monte->Write();
cosppi0_beamphoton_energysum_monte->Write();

m2pi0_beamphoton_energysum_monte->Write();
mppi0_beamphoton_energysum_monte->Write();

//reconstructed
cos2pi0_beamphoton_rek->Write();
cosppi0_beamphoton_energysum_rek->Write();

m2pi0_beamphoton_energysum_rek->Write();
mppi0_beamphoton_energysum_rek->Write();

curDir1->cd();
poltable_energy->Write();
poltable_energy_weight->Write();
filenumber.Write();
// TDirectory* curDir7  = curDir1->mkdir("KinFit");
// curDir7->cd();
// CL->Write();
// pulls_4g_CB->Write();
// pulls_4g_TAPS->Write();
// pulls_beam->Write();
// pulls_proton_CB->Write();
// pulls_proton_TAPS->Write();
// coplanarity_fit_vs_rek->Write();
// thetaproton_fit_vs_rek_TAPS->Write();
// thetaproton_fit_vs_rek_CB->Write();
// openangle_pi01_to_pi02_energy->Write();
// invmass_diff_rec_kin->Write();
// curDir1->cd();
// TDirectory* curDir3  = curDir1->mkdir("Selektion");
// curDir3->cd();
// time_prompt->Write();
// time1->Write();
// time_side->Write();
// 
// missingmass_collerated->Write();
// missingmass_inv_collerated->Write();
// 
// massesumme_collerated->Write();
// massesumme_mass_collerated->Write();
// massegegenmasse_collerated->Write();
// massegegenmasse_mass_collerated->Write();


curDir1->cd();
TDirectory* curDir5  = curDir1->mkdir("Selektion_withProton");
curDir5->cd();
coplanarity_collerated->Write();
coplanarity_mass_theta_collerated->Write();
coplanarity_mass_theta_inv_collerated->Write();
coplanarity_mass_collerated->Write();
coplanarity_theta_collerated->Write();
coplanarity_inv_collerated->Write();
// 
missingmass_theta_copl_collerated->Write();
missingmass_theta_copl_inv_collerated->Write();
missingmass_theta_collerated->Write();
missingmass_copl_collerated->Write();
missingmass_inv_collerated->Write();
missingmass_inv_collerated_proton->Write();
missingmass_collerated_proton->Write();
//
massesumme_mass_theta_copl_inv_beam_collerated->Write();
massesumme_collerated_proton->Write();
massesumme_mass_collerated_proton->Write();
massesumme_mass_theta_copl_collerated->Write();
massesumme_mass_theta_copl_inv1_collerated->Write();
massesumme_mass_theta_copl_inv_collerated->Write();
massesumme_mass_copl_collerated->Write();
massesumme_copl_collerated->Write();
massesumme_theta_collerated->Write();

massegegenmasse_collerated_proton->Write();
massegegenmasse_mass_copl_collerated->Write();
massegegenmasse_mass_theta_copl_collerated->Write();
massegegenmasse_copl_collerated->Write();
massegegenmasse_theta_collerated->Write();
massegegenmasse_mass_collerated_proton->Write();
massegegenmasse_allcutsandinv_collerated->Write();
massegegenmasse_allcuts_collerated->Write();

//
thetaproton_collerated->Write();
thetaproton_mass_copl_collerated->Write();
thetaproton_mass_copl_inv_collerated->Write();
thetaproton_mass_collerated->Write();
thetaproton_copl_collerated->Write();
thetaproton_inv_collerated->Write();
thetaproton_taps_collerated->Write();
thetaproton_mass_copl_taps_collerated->Write();
thetaproton_mass_copl_inv_taps_collerated->Write();
thetaproton_mass_taps_collerated->Write();
thetaproton_copl_taps_collerated->Write();
thetaproton_inv_taps_collerated->Write();
// Check_CBdE_E->Write();
// Check_TAPSdE_E->Write();
// 
// 
// curDir1->cd();
// TDirectory* curDir2  = curDir1->mkdir("Allgemein");
// curDir2->cd();
// cosverteilung_collerated->Write();
// cosppi0verteilung_collerated->Write(); 
// invverteilung_collerated->Write();
// invppi0verteilung_collerated->Write();
// 
// kristallminus_targetplus_collerated->Write();
// kristallminus_targetminus_collerated->Write();
// kristallplus_targetplus_collerated->Write();
// kristallplus_targetminus_collerated->Write();
// 
// kristallminus_targetplus_collerated_pb->Write();
// kristallminus_targetminus_collerated_pb->Write();
// kristallplus_targetplus_collerated_pb->Write();
// kristallplus_targetminus_collerated_pb->Write();
// 
// kristallminus_targetplus_collerated_pt->Write();
// kristallminus_targetminus_collerated_pt->Write();
// kristallplus_targetplus_collerated_pt->Write();
// kristallplus_targetminus_collerated_pt->Write();
// 
// kristallminus_targetplus_ppi0_collerated->Write();
// kristallminus_targetminus_ppi0_collerated->Write();
// kristallplus_targetplus_ppi0_collerated->Write();
// kristallplus_targetminus_ppi0_collerated->Write();
// 
// kristallminus_targetplus_ppi0_collerated_pb->Write();
// kristallminus_targetminus_ppi0_collerated_pb->Write();
// kristallplus_targetplus_ppi0_collerated_pb->Write();
// kristallplus_targetminus_ppi0_collerated_pb->Write();
// 
// kristallminus_targetplus_ppi0_collerated_pt->Write();
// kristallminus_targetminus_ppi0_collerated_pt->Write();
// kristallplus_targetplus_ppi0_collerated_pt->Write();
// kristallplus_targetminus_ppi0_collerated_pt->Write();
// 
// kristallminus_targetplus_im_collerated->Write();
// kristallminus_targetminus_im_collerated->Write();
// kristallplus_targetplus_im_collerated->Write();
// kristallplus_targetminus_im_collerated->Write();
// 
// kristallminus_targetplus_im_collerated_pb->Write();
// kristallminus_targetminus_im_collerated_pb->Write();
// kristallplus_targetplus_im_collerated_pb->Write();
// kristallplus_targetminus_im_collerated_pb->Write();
// 
// kristallminus_targetplus_im_collerated_pt->Write();
// kristallminus_targetminus_im_collerated_pt->Write();
// kristallplus_targetplus_im_collerated_pt->Write();
// kristallplus_targetminus_im_collerated_pt->Write();
// 
// kristallminus_targetplus_im_ppi0_collerated->Write();
// kristallminus_targetminus_im_ppi0_collerated->Write();
// kristallplus_targetplus_im_ppi0_collerated->Write();
// kristallplus_targetminus_im_ppi0_collerated->Write();
// 
// kristallminus_targetplus_im_ppi0_collerated_pb->Write();
// kristallminus_targetminus_im_ppi0_collerated_pb->Write();
// kristallplus_targetplus_im_ppi0_collerated_pb->Write();
// kristallplus_targetminus_im_ppi0_collerated_pb->Write();
// 
// kristallminus_targetplus_im_ppi0_collerated_pt->Write();
// kristallminus_targetminus_im_ppi0_collerated_pt->Write();
// kristallplus_targetplus_im_ppi0_collerated_pt->Write();
// kristallplus_targetminus_im_ppi0_collerated_pt->Write();
// 
// missingmassverteilung_collerated->Write();
// events_witherror->Write();
// events_all->Write();
// 
// curDir1->cd();
// TDirectory* curDir4  = curDir1->mkdir("Allgemein_withProton");
// curDir4->cd();
cosverteilung_collerated_proton->Write();
// cosppi0verteilung_collerated_proton->Write(); 
// invverteilung_collerated_proton->Write();
// invppi0verteilung_collerated_proton->Write();
// 
// kristallminus_targetplus_collerated_proton->Write();
// kristallminus_targetminus_collerated_proton->Write();
// kristallplus_targetplus_collerated_proton->Write();
// kristallplus_targetminus_collerated_proton->Write();
// 
// kristallminus_targetplus_collerated_pb_proton->Write();
// kristallminus_targetminus_collerated_pb_proton->Write();
// kristallplus_targetplus_collerated_pb_proton->Write();
// kristallplus_targetminus_collerated_pb_proton->Write();
// 
// kristallminus_targetplus_collerated_pt_proton->Write();
// kristallminus_targetminus_collerated_pt_proton->Write();
// kristallplus_targetplus_collerated_pt_proton->Write();
// kristallplus_targetminus_collerated_pt_proton->Write();
// 
// kristallminus_targetplus_ppi0_collerated->Write();
// kristallminus_targetminus_ppi0_collerated->Write();
// kristallplus_targetplus_ppi0_collerated->Write();
// kristallplus_targetminus_ppi0_collerated->Write();
// 
// kristallminus_targetplus_ppi0_collerated_pb->Write();
// kristallminus_targetminus_ppi0_collerated_pb->Write();
// kristallplus_targetplus_ppi0_collerated_pb->Write();
// kristallplus_targetminus_ppi0_collerated_pb->Write();
// 
// kristallminus_targetplus_ppi0_collerated_pt->Write();
// kristallminus_targetminus_ppi0_collerated_pt->Write();
// kristallplus_targetplus_ppi0_collerated_pt->Write();
// kristallplus_targetminus_ppi0_collerated_pt->Write();
// 
// kristallminus_targetplus_im_collerated_proton->Write();
// kristallminus_targetminus_im_collerated_proton->Write();
// kristallplus_targetplus_im_collerated_proton->Write();
// kristallplus_targetminus_im_collerated_proton->Write();
// 
// kristallminus_targetplus_im_collerated_pb_proton->Write();
// kristallminus_targetminus_im_collerated_pb_proton->Write();
// kristallplus_targetplus_im_collerated_pb_proton->Write();
// kristallplus_targetminus_im_collerated_pb_proton->Write();
// 
// kristallminus_targetplus_im_collerated_pt_proton->Write();
// kristallminus_targetminus_im_collerated_pt_proton->Write();
// kristallplus_targetplus_im_collerated_pt_proton->Write();
// kristallplus_targetminus_im_collerated_pt_proton->Write();
// 
// kristallminus_targetplus_im_ppi0_collerated_proton->Write();
// kristallminus_targetminus_im_ppi0_collerated_proton->Write();
// kristallplus_targetplus_im_ppi0_collerated_proton->Write();
// kristallplus_targetminus_im_ppi0_collerated_proton->Write();
// 
// kristallminus_targetplus_im_ppi0_collerated_pb_proton->Write();
// kristallminus_targetminus_im_ppi0_collerated_pb_proton->Write();
// kristallplus_targetplus_im_ppi0_collerated_pb_proton->Write();
// kristallplus_targetminus_im_ppi0_collerated_pb_proton->Write();
// 
// kristallminus_targetplus_im_ppi0_collerated_pt_proton->Write();
// kristallminus_targetminus_im_ppi0_collerated_pt_proton->Write();
// kristallplus_targetplus_im_ppi0_collerated_pt_proton->Write();
// kristallplus_targetminus_im_ppi0_collerated_pt_proton->Write();
// 
// thetaverteilung_collerated->Write();
// thetaverteilung_collerated_taps->Write();
// coplanarityverteilung_collerated->Write();
// missingmassverteilung_collerated_proton->Write();

outputFile->Close();
}


void	P2Pi0Analyse::ProcessScalerRead()
{
    //time.ScalerReadCorrection(5);
}
