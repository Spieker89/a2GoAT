#include "P2Pi0Analyse.h"

P2Pi0Analyse::P2Pi0Analyse()
{ 

time_prompt = new TH1F("time_prompt", 	"time_prompt", 	1400, -700, 700);
time1= new TH1F("time1", 	"time1", 	1400, -700, 700);
time_side 	= new TH1F("time_side", 	"time_side", 	1400, -700, 700);
   
 
poltable_energy          = new TH1F("poltable_energy",         "poltable_energy",           352,   0, 1448);
poltable_energy_weight          = new TH1F("poltable_energy_weight",         "poltable_energy_weight",           352,   0, 1448);

test = new TH1F("test","test",1100,-5.5,5.5);



events_all = new TH3F("events_all","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
events_witherror = new TH3F("events_witherror","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);

//........................................Histogramme für Coplanarity.............................

coplanarity_collerated = new TH1F("coplanarity_collerated", "Azimutwinkel (Signal); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_mass_theta_collerated = new TH1F("coplanarity_mass_theta_collerated", "Azimutwinkel (Cut: PolarWinkel+MissingMass(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_mass_theta_inv_collerated = new TH1F("coplanarity_mass_theta_inv_collerated", "Azimutwinkel (Cut: All(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_mass_collerated = new TH1F("coplanarity_mass_collerated", "Azimutwinkel (Cut: MissingMass(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_theta_collerated = new TH1F("coplanarity_theta_collerated", "Azimutwinkel (Cut: PolarWinkel+MissingMass(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);
coplanarity_inv_collerated = new TH1F("coplanarity_inv_collerated", "Azimutwinkel (Cut: invariante Masse(Signal)); #phi_{2#pi}-#phi_{p}[deg]", 400,0,360);

//..........................Histogramme für ThetaProton Polarwinkel...........................

thetaproton_collerated = new TH1F("thetaproton_collerated", "Polarwinkel (Signal); #theta_{rek}-#theta_{meas} [deg]", 400,-200,200);
thetaproton_mass_copl_collerated = new TH1F("thetaproton_mass_copl_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_copl_inv_collerated = new TH1F("thetaproton_mass_copl_inv_collerated", "Polarwinkel (Cut: All(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_collerated = new TH1F("thetaproton_mass_collerated", "Polarwinkel (Cut: Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_copl_collerated = new TH1F("thetaproton_copl_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_inv_collerated = new TH1F("thetaproton_inv_collerated", "Polarwinkel (Cut: invariante Masse(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);

//detektiert in TAPS
thetaproton_taps_collerated = new TH1F("thetaproton_taps_collerated", "Polarwinkel (Signal); #theta_{rek}-#theta_{meas} [deg]", 400,-200,200);
thetaproton_mass_copl_taps_collerated = new TH1F("thetaproton_mass_copl_taps_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_copl_inv_taps_collerated = new TH1F("thetaproton_mass_copl_inv_taps_collerated", "Polarwinkel (Cut: All(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_mass_taps_collerated = new TH1F("thetaproton_mass_taps_collerated", "Polarwinkel (Cut: Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_copl_taps_collerated = new TH1F("thetaproton_copl_taps_collerated", "Polarwinkel (Cut: Azimutwinkel+Missing Mass(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);
thetaproton_inv_taps_collerated = new TH1F("thetaproton_inv_taps_collerated", "Polarwinkel (Cut: invariante Masse(Signal)); #theta_{rek}-#theta_{meas}", 400,-200,200);


//.............................Histogramme für Missing Mass................................

missingmass_collerated = new TH1F("missingmass_collerated", "Missing Mass (Signal);m_{mm} [MeV]", 500,0,2200);
missingmass_theta_copl_collerated = new TH1F("missingmass_theta_copl_collerated", "Missing Mass (Cut: Polar+Azimut-Winkel(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_theta_copl_inv_collerated = new TH1F("missingmass_theta_copl_inv_collerated", "Missing Mass (Cut: All(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_theta_collerated = new TH1F("missingmass_theta_collerated", "Missing Mass (Cut: Polarwinkel(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_copl_collerated = new TH1F("missingmass_copl_collerated", "Missing Mass (Cut: Azimutalwinkel(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_inv_collerated = new TH1F("missingmass_inv_collerated", "Missing Mass (Cut: invariante Masse(Signal));m_{mm} [MeV]", 500,0,2200);

//........................Histogramme für invariante Masse.....................................
massesumme_mass_theta_copl_inv_beam_collerated = new TH2F("massesumme_mass_theta_copl_inv_beam_collerated","inv Mass in dependency of E_{beam}; m_{#gamma #gamma} [MeV]; E^{rec}_{#gamma} [MeV]",400,0,600,200,200,800);
massesumme_mass_theta_copl_inv_beam_collerated->Sumw2();
massesumme_collerated = new TH1F("massesumme_collerated", "Invariante Masse von der Summe (Signal); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_theta_copl_collerated = new TH1F("massesumme_mass_theta_copl_collerated", "Invariante Masse von der Summe (Cut: All(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_theta_copl_inv1_collerated = new TH1F("massesumme_mass_theta_copl_inv1_collerated", "Invariante Masse von der Summe (Cut: All+Cut auf invariante Masse 1(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_theta_copl_inv_collerated = new TH1F("massesumme_mass_theta_copl_inv_collerated", "Invariante Masse von der Summe (Cut: All+Cut auf invariante Masse 1(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_copl_collerated = new TH1F("massesumme_mass_copl_collerated", "Invariante Masse von der Summe (Cut: Azimutwinkel+Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_collerated = new TH1F("massesumme_mass_collerated", "Invariante Masse von der Summe (Cut: Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_copl_collerated = new TH1F("massesumme_copl_collerated", "Invariante Masse von der Summe (Cut: Azimutwinkel(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_theta_collerated = new TH1F("massesumme_theta_collerated", "Invariante Masse von der Summe (Cut: Polarwinkel(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_collerated_proton = new TH1F("massesumme_mass_collerated_proton", "Invariante Masse von der Summe (Cut: Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);


massegegenmasse_collerated = new TH2F("massegegenmasse_collerated","invariante Massen gegeneinander (Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_mass_copl_collerated = new TH2F("massegegenmasse_mass_copl_collerated","invariante Massen gegeneinander(Cut: Missing Mass+Azimutalwinkel(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_mass_theta_copl_collerated = new TH2F("massegegenmasse_mass_theta_copl_collerated","invariante Massen gegeneinander(Cut: All(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_copl_collerated = new TH2F("massegegenmasse_copl_collerated","invariante Massen gegeneinander(Cut: Azimutalwinkel(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_theta_collerated = new TH2F("massegegenmasse_theta_collerated","invariante Massen gegeneinander(Cut: Polarwinkel(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_mass_collerated = new TH2F("massegegenmasse_mass_collerated","invariante Massen gegeneinander(Cut: Missing Mass(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);
massegegenmasse_mass_collerated_proton = new TH2F("massegegenmasse_mass_collerated_proton","invariante Massen gegeneinander(Cut: Missing Mass(Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

//ANNNNNNNNNAAAAAAALLLLLLLLLLYYYYYYYYYYYSSSSSSSSSSSSIIIIIIIIISSSSSSSSSSSSSSSSSSS

cosverteilung_collerated = new TH2F("cosverteilung_collerated", "Cos-Verteilung; cos(#theta_{2#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosverteilung_collerated->Sumw2();

cosppi0verteilung_collerated = new TH2F("cosppi0verteilung_collerated", "Cos-Verteilung; cos(#theta_{p#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosppi0verteilung_collerated->Sumw2();

invverteilung_collerated = new TH2F("invverteilung_collerated", "Cos-Verteilung; m_{2#pi}[MeV];E_{beam} [MeV]", 100,300,800,200,200,800);
invverteilung_collerated->Sumw2();

invppi0verteilung_collerated = new TH2F("invppi0verteilung_collerated", "Cos-Verteilung; m_{p#pi}[MeV];E_{beam} [MeV]", 100,1100,1600,200,200,800);
invppi0verteilung_collerated->Sumw2();

cosverteilung_thetaproton_collerated = new TH2F("cosverteilung_thetaproton_collerated", "Cos-Verteilung; cos(#theta_{2#pi});#theta_{proton} [deg]", 72,-1,1,180,0,180);
cosverteilung_thetaproton_collerated->Sumw2();

thetaverteilung_collerated_taps = new TH2F("thetaverteilung_collerated_taps", "Polar-Verteilung;#theta_{rek}-#theta_{meas} [deg]", 400,-200,200,200,200,800);
thetaverteilung_collerated_taps->Sumw2();

thetaverteilung_collerated = new TH2F("thetaverteilung_collerated", "Polar-Verteilung;#theta_{rek}-#theta_{meas} [deg]", 400,-200,200,200,200,800);
thetaverteilung_collerated->Sumw2();

coplanarityverteilung_collerated = new TH2F("coplanarityverteilung_collerated", "Polar-Verteilung;#phi_{2#pi}-#phi_{p}[deg]", 400,0,360,200,200,800);
coplanarityverteilung_collerated->Sumw2();

missingmassverteilung_collerated = new TH2F("missingmassverteilung_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,200,200,800);
missingmassverteilung_collerated->Sumw2();

missingmassverteilung_collerated_proton = new TH2F("missingmassverteilung_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,200,200,800);
missingmassverteilung_collerated_proton->Sumw2();

//Histogramme für 2pi0 System mit Cos-Verteilung___________________________________________________

kristallminus_targetplus_collerated = new TH3F("kristallminus_targetplus_collerated","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated->Sumw2();

kristallminus_targetminus_collerated = new TH3F("kristallminus_targetminus_collerated","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated->Sumw2();

kristallplus_targetplus_collerated = new TH3F("kristallplus_targetplus_collerated","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated->Sumw2();

kristallplus_targetminus_collerated = new TH3F("kristallplus_targetminus_collerated","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb = new TH3F("kristallminus_targetplus_collerated_pb","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pb->Sumw2();

kristallminus_targetminus_collerated_pb = new TH3F("kristallminus_targetminus_collerated_pb","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pb->Sumw2();

kristallplus_targetplus_collerated_pb = new TH3F("kristallplus_targetplus_collerated_pb","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pb->Sumw2();

kristallplus_targetminus_collerated_pb = new TH3F("kristallplus_targetminus_collerated_pb","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pb->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt = new TH3F("kristallminus_targetplus_collerated_pt","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pt->Sumw2();

kristallminus_targetminus_collerated_pt = new TH3F("kristallminus_targetminus_collerated_pt","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pt->Sumw2();

kristallplus_targetplus_collerated_pt = new TH3F("kristallplus_targetplus_collerated_pt","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pt->Sumw2();

kristallplus_targetminus_collerated_pt = new TH3F("kristallplus_targetminus_collerated_pt","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pt->Sumw2();

//Histogramme für ppi0 System  mit Cos-Verteilung__________________________________________________________

kristallminus_targetplus_ppi0_collerated = new TH3F("kristallminus_targetplus_ppi0_collerated","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_ppi0_collerated->Sumw2();

kristallminus_targetminus_ppi0_collerated = new TH3F("kristallminus_targetminus_ppi0_collerated","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_ppi0_collerated->Sumw2();

kristallplus_targetplus_ppi0_collerated = new TH3F("kristallplus_targetplus_ppi0_collerated","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_ppi0_collerated->Sumw2();

kristallplus_targetminus_ppi0_collerated = new TH3F("kristallplus_targetminus_ppi0_collerated","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_ppi0_collerated->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_ppi0_collerated_pb = new TH3F("kristallminus_targetplus_ppi0_collerated_pb","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_ppi0_collerated_pb->Sumw2();

kristallminus_targetminus_ppi0_collerated_pb = new TH3F("kristallminus_targetminus_ppi0_collerated_pb","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_ppi0_collerated_pb->Sumw2();

kristallplus_targetplus_ppi0_collerated_pb = new TH3F("kristallplus_targetplus_ppi0_collerated_pb","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_ppi0_collerated_pb->Sumw2();

kristallplus_targetminus_ppi0_collerated_pb = new TH3F("kristallplus_targetminus_ppi0_collerated_pb","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_ppi0_collerated_pb->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_ppi0_collerated_pt = new TH3F("kristallminus_targetplus_ppi0_collerated_pt","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_ppi0_collerated_pt->Sumw2();

kristallminus_targetminus_ppi0_collerated_pt = new TH3F("kristallminus_targetminus_ppi0_collerated_pt","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_ppi0_collerated_pt->Sumw2();

kristallplus_targetplus_ppi0_collerated_pt = new TH3F("kristallplus_targetplus_ppi0_collerated_pt","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_ppi0_collerated_pt->Sumw2();

kristallplus_targetminus_ppi0_collerated_pt = new TH3F("kristallplus_targetminus_ppi0_collerated_pt","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_ppi0_collerated_pt->Sumw2();

//Histogramme für 2pi0 System mit invM-Verteilung___________________________________________________

kristallminus_targetplus_im_collerated = new TH3F("kristallminus_targetplus_im_collerated","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetplus_im_collerated->Sumw2();

kristallminus_targetminus_im_collerated = new TH3F("kristallminus_targetminus_im_collerated","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetminus_im_collerated->Sumw2();

kristallplus_targetplus_im_collerated = new TH3F("kristallplus_targetplus_im_collerated","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetplus_im_collerated->Sumw2();

kristallplus_targetminus_im_collerated = new TH3F("kristallplus_targetminus_im_collerated","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetminus_im_collerated->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_im_collerated_pb = new TH3F("kristallminus_targetplus_im_collerated_pb","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetplus_im_collerated_pb->Sumw2();

kristallminus_targetminus_im_collerated_pb = new TH3F("kristallminus_targetminus_im_collerated_pb","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetminus_im_collerated_pb->Sumw2();

kristallplus_targetplus_im_collerated_pb = new TH3F("kristallplus_targetplus_im_collerated_pb","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetplus_im_collerated_pb->Sumw2();

kristallplus_targetminus_im_collerated_pb = new TH3F("kristallplus_targetminus_im_collerated_pb","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetminus_im_collerated_pb->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_im_collerated_pt = new TH3F("kristallminus_targetplus_im_collerated_pt","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetplus_im_collerated_pt->Sumw2();

kristallminus_targetminus_im_collerated_pt = new TH3F("kristallminus_targetminus_im_collerated_pt","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetminus_im_collerated_pt->Sumw2();

kristallplus_targetplus_im_collerated_pt = new TH3F("kristallplus_targetplus_im_collerated_pt","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetplus_im_collerated_pt->Sumw2();

kristallplus_targetminus_im_collerated_pt = new TH3F("kristallplus_targetminus_im_collerated_pt","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetminus_im_collerated_pt->Sumw2();

//Histogramme für ppi0 System mit invM-Verteilung___________________________________________________________________

kristallminus_targetplus_im_ppi0_collerated = new TH3F("kristallminus_targetplus_im_ppi0_collerated","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetplus_im_ppi0_collerated->Sumw2();

kristallminus_targetminus_im_ppi0_collerated = new TH3F("kristallminus_targetminus_im_ppi0_collerated","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetminus_im_ppi0_collerated->Sumw2();

kristallplus_targetplus_im_ppi0_collerated = new TH3F("kristallplus_targetplus_im_ppi0_collerated","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetplus_im_ppi0_collerated->Sumw2();

kristallplus_targetminus_im_ppi0_collerated = new TH3F("kristallplus_targetminus_im_ppi0_collerated","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetminus_im_ppi0_collerated->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_im_ppi0_collerated_pb = new TH3F("kristallminus_targetplus_im_ppi0_collerated_pb","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetplus_im_ppi0_collerated_pb->Sumw2();

kristallminus_targetminus_im_ppi0_collerated_pb = new TH3F("kristallminus_targetminus_im_ppi0_collerated_pb","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetminus_im_ppi0_collerated_pb->Sumw2();

kristallplus_targetplus_im_ppi0_collerated_pb = new TH3F("kristallplus_targetplus_im_ppi0_collerated_pb","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetplus_im_ppi0_collerated_pb->Sumw2();

kristallplus_targetminus_im_ppi0_collerated_pb = new TH3F("kristallplus_targetminus_im_ppi0_collerated_pb","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetminus_im_ppi0_collerated_pb->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_im_ppi0_collerated_pt = new TH3F("kristallminus_targetplus_im_ppi0_collerated_pt","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetplus_im_ppi0_collerated_pt->Sumw2();

kristallminus_targetminus_im_ppi0_collerated_pt = new TH3F("kristallminus_targetminus_im_ppi0_collerated_pt","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetminus_im_ppi0_collerated_pt->Sumw2();

kristallplus_targetplus_im_ppi0_collerated_pt = new TH3F("kristallplus_targetplus_im_ppi0_collerated_pt","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetplus_im_ppi0_collerated_pt->Sumw2();

kristallplus_targetminus_im_ppi0_collerated_pt = new TH3F("kristallplus_targetminus_im_ppi0_collerated_pt","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetminus_im_ppi0_collerated_pt->Sumw2();

//PROTON IDENTIFIED

cosverteilung_collerated_proton = new TH2F("cosverteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{2#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosverteilung_collerated_proton->Sumw2();

cosppi0verteilung_collerated_proton = new TH2F("cosppi0verteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{p#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosppi0verteilung_collerated_proton->Sumw2();

invverteilung_collerated_proton = new TH2F("invverteilung_collerated_proton", "Cos-Verteilung; m_{2#pi}[MeV];E_{beam} [MeV]", 100,300,800,200,200,800);
invverteilung_collerated_proton->Sumw2();

invppi0verteilung_collerated_proton = new TH2F("invppi0verteilung_collerated_proton", "Cos-Verteilung; m_{p#pi}[MeV];E_{beam} [MeV]", 100,1100,1600,200,200,800);
invppi0verteilung_collerated_proton->Sumw2();

cosverteilung_thetaproton_collerated_proton = new TH2F("cosverteilung_thetaproton_collerated_proton", "Cos-Verteilung; cos(#theta_{2#pi});#theta_{proton} [deg]", 72,-1,1,180,0,180);
cosverteilung_thetaproton_collerated_proton->Sumw2();

//Histogramme für 2pi0 System mit Cos-Verteilung___________________________________________________

kristallminus_targetplus_collerated_proton = new TH3F("kristallminus_targetplus_collerated_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_proton->Sumw2();

kristallminus_targetminus_collerated_proton = new TH3F("kristallminus_targetminus_collerated_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_proton->Sumw2();

kristallplus_targetplus_collerated_proton = new TH3F("kristallplus_targetplus_collerated_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_proton->Sumw2();

kristallplus_targetminus_collerated_proton = new TH3F("kristallplus_targetminus_collerated_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_proton->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb_proton = new TH3F("kristallminus_targetplus_collerated_pb_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pb_proton->Sumw2();

kristallminus_targetminus_collerated_pb_proton = new TH3F("kristallminus_targetminus_collerated_pb_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pb_proton->Sumw2();

kristallplus_targetplus_collerated_pb_proton = new TH3F("kristallplus_targetplus_collerated_pb_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pb_proton->Sumw2();

kristallplus_targetminus_collerated_pb_proton = new TH3F("kristallplus_targetminus_collerated_pb_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pb_proton->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt_proton = new TH3F("kristallminus_targetplus_collerated_pt_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pt_proton->Sumw2();

kristallminus_targetminus_collerated_pt_proton = new TH3F("kristallminus_targetminus_collerated_pt_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pt_proton->Sumw2();

kristallplus_targetplus_collerated_pt_proton = new TH3F("kristallplus_targetplus_collerated_pt_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pt_proton->Sumw2();

kristallplus_targetminus_collerated_pt_proton = new TH3F("kristallplus_targetminus_collerated_pt_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; cos(#theta)_{2#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pt_proton->Sumw2();


//Histogramme für ppi0 System  mit Cos-Verteilung__________________________________________________________

kristallminus_targetplus_ppi0_collerated_proton = new TH3F("kristallminus_targetplus_ppi0_collerated_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_ppi0_collerated_proton->Sumw2();

kristallminus_targetminus_ppi0_collerated_proton = new TH3F("kristallminus_targetminus_ppi0_collerated_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_ppi0_collerated_proton->Sumw2();

kristallplus_targetplus_ppi0_collerated_proton = new TH3F("kristallplus_targetplus_ppi0_collerated_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_ppi0_collerated_proton->Sumw2();

kristallplus_targetminus_ppi0_collerated_proton = new TH3F("kristallplus_targetminus_ppi0_collerated_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_ppi0_collerated_proton->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_ppi0_collerated_pb_proton = new TH3F("kristallminus_targetplus_ppi0_collerated_pb_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_ppi0_collerated_pb_proton->Sumw2();

kristallminus_targetminus_ppi0_collerated_pb_proton = new TH3F("kristallminus_targetminus_ppi0_collerated_pb_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_ppi0_collerated_pb_proton->Sumw2();

kristallplus_targetplus_ppi0_collerated_pb_proton = new TH3F("kristallplus_targetplus_ppi0_collerated_pb_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_ppi0_collerated_pb_proton->Sumw2();

kristallplus_targetminus_ppi0_collerated_pb_proton = new TH3F("kristallplus_targetminus_ppi0_collerated_pb_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_ppi0_collerated_pb_proton->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_ppi0_collerated_pt_proton = new TH3F("kristallminus_targetplus_ppi0_collerated_pt_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_ppi0_collerated_pt_proton->Sumw2();

kristallminus_targetminus_ppi0_collerated_pt_proton = new TH3F("kristallminus_targetminus_ppi0_collerated_pt_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_ppi0_collerated_pt_proton->Sumw2();

kristallplus_targetplus_ppi0_collerated_pt_proton = new TH3F("kristallplus_targetplus_ppi0_collerated_pt_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_ppi0_collerated_pt_proton->Sumw2();

kristallplus_targetminus_ppi0_collerated_pt_proton = new TH3F("kristallplus_targetminus_ppi0_collerated_pt_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; cos(#theta)_{p#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_ppi0_collerated_pt_proton->Sumw2();

//Histogramme für 2pi0 System mit invM-Verteilung___________________________________________________

kristallminus_targetplus_im_collerated_proton = new TH3F("kristallminus_targetplus_im_collerated_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetplus_im_collerated_proton->Sumw2();

kristallminus_targetminus_im_collerated_proton = new TH3F("kristallminus_targetminus_im_collerated_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetminus_im_collerated_proton->Sumw2();

kristallplus_targetplus_im_collerated_proton = new TH3F("kristallplus_targetplus_im_collerated_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetplus_im_collerated_proton->Sumw2();

kristallplus_targetminus_im_collerated_proton = new TH3F("kristallplus_targetminus_im_collerated_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetminus_im_collerated_proton->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_im_collerated_pb_proton = new TH3F("kristallminus_targetplus_im_collerated_pb_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetplus_im_collerated_pb_proton->Sumw2();

kristallminus_targetminus_im_collerated_pb_proton = new TH3F("kristallminus_targetminus_im_collerated_pb_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetminus_im_collerated_pb_proton->Sumw2();

kristallplus_targetplus_im_collerated_pb_proton = new TH3F("kristallplus_targetplus_im_collerated_pb_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetplus_im_collerated_pb_proton->Sumw2();

kristallplus_targetminus_im_collerated_pb_proton = new TH3F("kristallplus_targetminus_im_collerated_pb_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetminus_im_collerated_pb_proton->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_im_collerated_pt_proton = new TH3F("kristallminus_targetplus_im_collerated_pt_proton","Kristall-45 Target positiv; #phi_{2#pi}[deg];m_{2#pi}[MeV];E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetplus_im_collerated_pt_proton->Sumw2();

kristallminus_targetminus_im_collerated_pt_proton = new TH3F("kristallminus_targetminus_im_collerated_pt_proton","Kristall-45 Target negativ; #phi_{2#pi}[deg];m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallminus_targetminus_im_collerated_pt_proton->Sumw2();

kristallplus_targetplus_im_collerated_pt_proton = new TH3F("kristallplus_targetplus_im_collerated_pt_proton","Kristall+45 Target positiv; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetplus_im_collerated_pt_proton->Sumw2();

kristallplus_targetminus_im_collerated_pt_proton = new TH3F("kristallplus_targetminus_im_collerated_pt_proton","Kristall+45 Target negativ; #phi_{2#pi}[deg]; m_{2#pi}[MeV]; E_{beam} [MeV]", 96,-180,180,100,300,800,200,200,800);
kristallplus_targetminus_im_collerated_pt_proton->Sumw2();

//Histogramme für ppi0 System mit invM-Verteilung___________________________________________________________________

kristallminus_targetplus_im_ppi0_collerated_proton = new TH3F("kristallminus_targetplus_im_ppi0_collerated_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetplus_im_ppi0_collerated_proton->Sumw2();

kristallminus_targetminus_im_ppi0_collerated_proton = new TH3F("kristallminus_targetminus_im_ppi0_collerated_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetminus_im_ppi0_collerated_proton->Sumw2();

kristallplus_targetplus_im_ppi0_collerated_proton = new TH3F("kristallplus_targetplus_im_ppi0_collerated_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetplus_im_ppi0_collerated_proton->Sumw2();

kristallplus_targetminus_im_ppi0_collerated_proton = new TH3F("kristallplus_targetminus_im_ppi0_collerated_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetminus_im_ppi0_collerated_proton->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_im_ppi0_collerated_pb_proton = new TH3F("kristallminus_targetplus_im_ppi0_collerated_pb_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetplus_im_ppi0_collerated_pb_proton->Sumw2();

kristallminus_targetminus_im_ppi0_collerated_pb_proton = new TH3F("kristallminus_targetminus_im_ppi0_collerated_pb_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetminus_im_ppi0_collerated_pb_proton->Sumw2();

kristallplus_targetplus_im_ppi0_collerated_pb_proton = new TH3F("kristallplus_targetplus_im_ppi0_collerated_pb_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetplus_im_ppi0_collerated_pb_proton->Sumw2();

kristallplus_targetminus_im_ppi0_collerated_pb_proton = new TH3F("kristallplus_targetminus_im_ppi0_collerated_pb_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetminus_im_ppi0_collerated_pb_proton->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_im_ppi0_collerated_pt_proton = new TH3F("kristallminus_targetplus_im_ppi0_collerated_pt_proton","Kristall-45 Target positiv; #phi_{p#pi}[deg];m_{p#pi};E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetplus_im_ppi0_collerated_pt_proton->Sumw2();

kristallminus_targetminus_im_ppi0_collerated_pt_proton = new TH3F("kristallminus_targetminus_im_ppi0_collerated_pt_proton","Kristall-45 Target negativ; #phi_{p#pi}[deg];m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallminus_targetminus_im_ppi0_collerated_pt_proton->Sumw2();

kristallplus_targetplus_im_ppi0_collerated_pt_proton = new TH3F("kristallplus_targetplus_im_ppi0_collerated_pt_proton","Kristall+45 Target positiv; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetplus_im_ppi0_collerated_pt_proton->Sumw2();

kristallplus_targetminus_im_ppi0_collerated_pt_proton = new TH3F("kristallplus_targetminus_im_ppi0_collerated_pt_proton","Kristall+45 Target negativ; #phi_{p#pi}[deg]; m_{p#pi}; E_{beam} [MeV]", 96,-180,180,100,1100,1600,200,200,800);
kristallplus_targetminus_im_ppi0_collerated_pt_proton->Sumw2();


Check_CBdE_E= new TH2F("Check_CBdE_E", "dE_E (all CB clusters compared to PID hits)", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E= new TH2F("Check_TAPSdE_E", "dE_E (all TAPS clusters compared to Veto hits)", 	400, 0, 400, 100, 0, 10);

}

P2Pi0Analyse::~P2Pi0Analyse()
{
}



Bool_t	P2Pi0Analyse::Start()
{
pt=P2Pi0Analyse::targetpol(inputFile);

// planesetting=P2Pi0Analyse::polplane(inputFile);
if(linpol->GetPolarizationPlane()==1){planesetting="PERP";}
if(linpol->GetPolarizationPlane()==0){planesetting="PARA";}
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

//:::::::::::::::::::::::::::::::::::::::::::::::::::Photonen festlegen::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
			
	TLorentzVector photon1_4vektor = photons->Particle(0);
	TLorentzVector photon2_4vektor = photons->Particle(1);
	TLorentzVector photon3_4vektor = photons->Particle(2);
	TLorentzVector photon4_4vektor = photons->Particle(3);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Vektoren definieren:::::::::::::::::::::::::::::::::::::::::::::::::::
			
	TVector3 photon1_3vektor = photon1_4vektor.Vect();
	TVector3 photon2_3vektor = photon2_4vektor.Vect();
	TVector3 photon3_3vektor = photon3_4vektor.Vect();
	TVector3 photon4_3vektor = photon4_4vektor.Vect();
			
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::Berechne invariante Masse:::::::::::::::::::::::::::::::::::::::::::::::::
			
	//Berechne invariante Masse von Photon 1 und 2			
	Double_t invariantemasse12 = invariantemasse(photon1_4vektor,photon2_4vektor);
	//Berechne invariante Masse von Photon 3 und 4	
	Double_t invariantemasse34 = invariantemasse(photon3_4vektor,photon4_4vektor);
		
	//Berechne invariante Masse von Photon 1 und 3
	Double_t invariantemasse13 = invariantemasse(photon1_4vektor,photon3_4vektor);
	//Berechne invariante Masse von Photon 2 und 4	
	Double_t invariantemasse24 = invariantemasse(photon2_4vektor,photon4_4vektor);

	//Berechne invariante Masse von Photon 1 und 4
	Double_t invariantemasse14 =invariantemasse(photon1_4vektor,photon4_4vektor);
	//Berechne invariante Masse von Photon 2 und 3	
	Double_t invariantemasse23 = invariantemasse(photon2_4vektor,photon3_4vektor);

//_______________________________passende Kombination finden____________________________________________	

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

	Double_t Chi12 = ChiPionPion(invariantemasse12,photon1_4vektor,photon2_4vektor,invariantemasse34,photon3_4vektor,photon4_4vektor);
	Double_t Chi13 = ChiPionPion(invariantemasse13,photon1_4vektor,photon3_4vektor,invariantemasse24,photon2_4vektor,photon4_4vektor);
	Double_t Chi14 = ChiPionPion(invariantemasse14,photon1_4vektor,photon4_4vektor,invariantemasse23,photon2_4vektor,photon3_4vektor);

	TLorentzVector pi01_4vect,pi02_4vect;

	if(Chi12 < Chi13 && Chi12 < Chi14){
		invariantemasse1 = invariantemasse12;
		pi01_4vect = photons->Particle(0)+photons->Particle(1);
		pion1_E = photon1_E + photon2_E ;
		pion1_pz = photon1_pz + photon2_pz;	
		invariantemasse2 = invariantemasse34;
		pi02_4vect = photons->Particle(2)+photons->Particle(3);
		pion2_E = photon3_E + photon4_E; 
		pion2_pz = photon3_pz + photon4_pz;
		Chi = Chi12;
	}
			
	if(Chi13 < Chi12 && Chi13 < Chi14){
		invariantemasse1 = invariantemasse13;
		pi01_4vect = photons->Particle(0)+photons->Particle(2);
		pion1_E = photon1_E + photon3_E ;
		pion1_pz = photon1_pz + photon3_pz;		
		invariantemasse2 = invariantemasse24;
		pi02_4vect = photons->Particle(1)+photons->Particle(3);
		pion2_E = photon2_E + photon4_E; 
		pion2_pz = photon2_pz + photon4_pz;
		Chi = Chi13;
	}
	
	if(Chi14 < Chi13 && Chi14 < Chi12){
		invariantemasse1 = invariantemasse14;
		pi01_4vect = photons->Particle(0)+photons->Particle(3);
		pion1_E = photon1_E + photon4_E ;
		pion1_pz = photon1_pz + photon4_pz;		
		invariantemasse2 = invariantemasse23;
		pi02_4vect = photons->Particle(1)+photons->Particle(2);;
		pion2_E = photon2_E + photon3_E; 
		pion2_pz = photon2_pz + photon3_pz;
		Chi = Chi14;
	}	

	pion1_m = invariantemasse1*invariantemasse1;
	pion2_m = invariantemasse2*invariantemasse2;

	Double_t Chi1[6];
				
	Chi1[0] = ChiPionEta(invariantemasse12,photon1_4vektor,photon2_4vektor,invariantemasse34,photon3_4vektor,photon3_4vektor);		
	Chi1[1] = ChiPionEta(invariantemasse34,photon3_4vektor,photon4_4vektor,invariantemasse12,photon1_4vektor,photon2_4vektor);

	Chi1[2]= ChiPionEta(invariantemasse13,photon1_4vektor,photon3_4vektor,invariantemasse24,photon2_4vektor,photon4_4vektor);
	Chi1[3]= ChiPionEta(invariantemasse24,photon2_4vektor,photon4_4vektor,invariantemasse13,photon1_4vektor,photon3_4vektor);

	Chi1[4]= ChiPionEta(invariantemasse14,photon1_4vektor,photon4_4vektor,invariantemasse23,photon2_4vektor,photon3_4vektor);
	Chi1[5]= ChiPionEta(invariantemasse23,photon2_4vektor,photon3_4vektor,invariantemasse14,photon1_4vektor,photon4_4vektor);

	Double_t vergleich=Chi1[0];

	for(Int_t i=1; i<6;i++){
		if(Chi1[i]<vergleich){vergleich = Chi1[i];}
	}

// 	cout << "FINALE: " << Chi << endl;		

		//TEST the target polarization (no entry with pt=0)
		test->Fill(pt);



		for(Int_t j=0; j < tagger->GetNTagged();j++)
		{
		if(Chi > vergleich){continue;}

		if(pt==5){continue;}

			Double_t time= tagger->GetTaggedTime(j) - 0.25*(photons->GetTime(1)+photons->GetTime(0)+photons->GetTime(2)+photons->GetTime(3));
			time1->Fill(time);

			if(time > -20 && time < 5){
				time_prompt->Fill(time);
			}

				if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
					time_side->Fill(time);
				}

			poltable_energy->Fill(tagger->GetTaggedEnergy(j));
			poltable_energy_weight->Fill(tagger->GetTaggedEnergy(j),linpol->GetPolarizationDegree(tagger->GetTaggedChannel(j)));

			//get target, beam and missng particle information
			TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
			TLorentzVector  beam_4vect = TLorentzVector(0.,0.,tagger->GetTaggedEnergy(j),tagger->GetTaggedEnergy(j));
			TLorentzVector  missingp_4vect = beam_4vect + protonvektor_target - pi02_4vect-pi01_4vect;
			Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
			Double_t beamphoton1E=tagger->GetTaggedEnergy(j);
			Double_t pb=linpol->GetPolarizationDegree(tagger->GetTaggedChannel(j));

			//due to kinematic not possible -> kick them out
			if(anglethetaproton_rek > 90){continue;}
			
			if(beamphoton1E<200 || beamphoton1E>800){continue;}

			//energy dependent cuts
			Double_t oben_copl=270.135+(-0.37891)*beamphoton1E+(0.000637701)*beamphoton1E*beamphoton1E+(-3.62165e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_copl=89.4344+(0.360771)*beamphoton1E+(-0.000593373)*beamphoton1E*beamphoton1E+(3.32871e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_theta=-52.7257+(0.352941)*beamphoton1E+(-0.000625629)*beamphoton1E*beamphoton1E+(3.63602e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_theta=38.2131+(-0.302196)*beamphoton1E+(0.00058914)*beamphoton1E*beamphoton1E+(-3.67153e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_mass=910.25+(0.375124)*beamphoton1E+(-0.000611761)*beamphoton1E*beamphoton1E+(3.88824e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_mass=989.045+(-0.466069)*beamphoton1E+(0.000689688)*beamphoton1E*beamphoton1E+(-4.82698e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_mass_proton=969.521+(0.0658678)*beamphoton1E+(-9.43064e-05)*beamphoton1E*beamphoton1E+(8.86854e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_mass_proton=912.747+(0.00410493)*beamphoton1E+(-0.000197333)*beamphoton1E*beamphoton1E+(6.67124e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			
			Double_t oben_inv=142.348+(0.0886059)*beamphoton1E+(-0.000129399)*beamphoton1E*beamphoton1E+(7.21352e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_inv=112.242+(-0.0110781)*beamphoton1E+(1.20014e-06)*beamphoton1E*beamphoton1E+(8.16401e-09)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t oben_inv_proton=152.662+(0.00115463)*beamphoton1E+(6.3428e-05)*beamphoton1E*beamphoton1E+(-5.08229e-08)*beamphoton1E*beamphoton1E*beamphoton1E;
			Double_t unten_inv_proton=109.617+(0.0221158)*beamphoton1E+(-6.68222e-05)*beamphoton1E*beamphoton1E+(4.07149e-08)*beamphoton1E*beamphoton1E*beamphoton1E;

			//Missing Mass
			Double_t missingmass=missingp_4vect.M();

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
			
			Double_t scalbackground_ppi0=-0.0625*0.5;
			
			Double_t scalbackgroundpb = -0.0625*pb;
			Double_t scalbackgroundpt = -0.0625*pt;
			
			Double_t scalbackgroundpb_ppi0 = -0.5*0.0625*pb;
			Double_t scalbackgroundpt_ppi0 = -0.5*0.0625*pt;
			
			Double_t scalsignalpb_ppi0 = 0.5*pb;
			Double_t scalsignalpt_ppi0 = 0.5*pt;

//...........................Histogramme vor jeglichen Cuts erstellen..................................
			
			
			if(time > -20 && time < 5 ){

				missingmass_collerated->Fill(missingmass);
				
				massesumme_collerated->Fill(invariantemasse34);
				massesumme_collerated->Fill(invariantemasse24);
				massesumme_collerated->Fill(invariantemasse23);	
				massesumme_collerated->Fill(invariantemasse12);
				massesumme_collerated->Fill(invariantemasse13);
				massesumme_collerated->Fill(invariantemasse14);
				
				massegegenmasse_collerated->Fill(invariantemasse12,invariantemasse34);
				massegegenmasse_collerated->Fill(invariantemasse34,invariantemasse12);
				massegegenmasse_collerated->Fill(invariantemasse24,invariantemasse13);
				massegegenmasse_collerated->Fill(invariantemasse13,invariantemasse24);
				massegegenmasse_collerated->Fill(invariantemasse14,invariantemasse23);
				massegegenmasse_collerated->Fill(invariantemasse23,invariantemasse14);
			}
			
			if((time > -300 && time < -100) || (time > 100 && time < 300)){	
				missingmass_collerated->Fill(missingmass,-0.0625);
				
				massesumme_collerated->Fill(invariantemasse34,-0.0625);
				massesumme_collerated->Fill(invariantemasse24,-0.0625);
				massesumme_collerated->Fill(invariantemasse23,-0.0625);	
				massesumme_collerated->Fill(invariantemasse12,-0.0625);
				massesumme_collerated->Fill(invariantemasse13,-0.0625);
				massesumme_collerated->Fill(invariantemasse14,-0.0625);
				
				massegegenmasse_collerated->Fill(invariantemasse12,invariantemasse34,-0.0625);
				massegegenmasse_collerated->Fill(invariantemasse34,invariantemasse12,-0.0625);
				massegegenmasse_collerated->Fill(invariantemasse24,invariantemasse13,-0.0625);
				massegegenmasse_collerated->Fill(invariantemasse13,invariantemasse24,-0.0625);
				massegegenmasse_collerated->Fill(invariantemasse14,invariantemasse23,-0.0625);
				massegegenmasse_collerated->Fill(invariantemasse23,invariantemasse14,-0.0625);
			}								
									
			
// 			::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::einzelne Cuts:::::::::::::::::::::::::::::::::::::::::
			
			
// 			__________________invariante Masse______________________________________________________
			
			if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
			
				if(time > -20 && time < 5){
			
					missingmass_inv_collerated->Fill(missingmass); 
					missingmassverteilung_collerated->Fill(missingmass, beamphoton1E); 

				}
					if((time > -300 && time < -100) ||(time > 100 && time < 300)){
	
						missingmass_inv_collerated->Fill(missingmass,-0.0625); 
						missingmassverteilung_collerated->Fill(missingmass, beamphoton1E,-0.0625); 

					}
			}

			if(missingmass>(unten_mass) && missingmass<(oben_mass)){
									
				if(time > -20 && time < 5){
					massesumme_mass_collerated->Fill(invariantemasse13);
					massesumme_mass_collerated->Fill(invariantemasse12);				
					massesumme_mass_collerated->Fill(invariantemasse14);
					massesumme_mass_collerated->Fill(invariantemasse34);
					massesumme_mass_collerated->Fill(invariantemasse24);
					massesumme_mass_collerated->Fill(invariantemasse23);
					massegegenmasse_mass_collerated->Fill(invariantemasse12,invariantemasse34);
					massegegenmasse_mass_collerated->Fill(invariantemasse13,invariantemasse24);
					massegegenmasse_mass_collerated->Fill(invariantemasse14,invariantemasse23);
					massegegenmasse_mass_collerated->Fill(invariantemasse34,invariantemasse12);
					massegegenmasse_mass_collerated->Fill(invariantemasse24,invariantemasse13);
					massegegenmasse_mass_collerated->Fill(invariantemasse23,invariantemasse14);
				}
					
						if((time > -300 && time < -100) ||(time > 100 && time < 300)){						
							massesumme_mass_collerated->Fill(invariantemasse13,-0.0625);
							massesumme_mass_collerated->Fill(invariantemasse12,-0.0625);				
							massesumme_mass_collerated->Fill(invariantemasse14,-0.0625);
							massesumme_mass_collerated->Fill(invariantemasse34,-0.0625);
							massesumme_mass_collerated->Fill(invariantemasse24,-0.0625);
							massesumme_mass_collerated->Fill(invariantemasse23,-0.0625);
							massegegenmasse_mass_collerated->Fill(invariantemasse12,invariantemasse34,-0.0625);
							massegegenmasse_mass_collerated->Fill(invariantemasse13,invariantemasse24,-0.0625);
							massegegenmasse_mass_collerated->Fill(invariantemasse14,invariantemasse23,-0.0625);
							massegegenmasse_mass_collerated->Fill(invariantemasse34,invariantemasse12,-0.0625);
							massegegenmasse_mass_collerated->Fill(invariantemasse24,invariantemasse13,-0.0625);
							massegegenmasse_mass_collerated->Fill(invariantemasse23,invariantemasse14,-0.0625);
						}
			}
			

//_____________________________only one charged particle_________________________________

				if(protons->GetNParticles()==1 || electrons->GetNParticles()==1 || chargedPions->GetNParticles()==1 || rootinos->GetNParticles()==1){

					//get charged particle information
					TLorentzVector proton_4vect_meas;
					if(electrons->GetNParticles()==1){proton_4vect_meas=electrons->Particle(0);}
					if(protons->GetNParticles()==1){proton_4vect_meas=protons->Particle(0);}
					if(chargedPions->GetNParticles()==1){proton_4vect_meas=chargedPions->Particle(0);}
					if(rootinos->GetNParticles()==1){proton_4vect_meas=rootinos->Particle(0);}

					//Phi-Difference
					Double_t anglephiproton = TMath::RadToDeg()*proton_4vect_meas.Vect().Phi();
					Double_t phi = phimeson - anglephiproton;
					if(phi < 0){phi = phi + 360;}
			
					//theta differenz
					Double_t anglethetaproton_meas = TMath::RadToDeg()*proton_4vect_meas.Vect().Theta();
					Double_t thetadiff=anglethetaproton_rek-anglethetaproton_meas;

					//due to kinematic not possible -> kick them out
					if(anglethetaproton_meas > 90){continue;}

					if(time > -20 && time < 5 ){
					
						coplanarity_collerated->Fill(phi);
						thetaproton_collerated->Fill(thetadiff);
				
					}
					
					if((time > -300 && time < -100) || (time > 100 && time < 300)){	
					
						coplanarity_collerated->Fill(phi,-0.0625);
						thetaproton_collerated->Fill(thetadiff,-0.0625);
					
					}	

		// 			::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::aufeinanderfolgende Cuts:::::::::::::::::::::::::::::::::::::::::
					
					
// 			--------------------------------------------------------CopLanaritY-Cuts----------------------
					
		// 			__________________Missing Mass______________________________________________________
					
					if(missingmass>(unten_mass) && missingmass<(oben_mass)){
					
						if(time > -20 && time < 5){
					
							coplanarity_mass_collerated->Fill(phi);
						}
							if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
									
								coplanarity_mass_collerated->Fill(phi,-0.0625);
							}	
						
		// 			__________________Polarwinkel______________________________________________________
					
						if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {
					
							if(time > -20 && time < 5){
					
								coplanarity_mass_theta_collerated->Fill(phi);
							}
								if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
						
									coplanarity_mass_theta_collerated->Fill(phi,-0.0625);
								}
						
		// 			__________________invariante Masse______________________________________________________
					
							if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
					
								if(time > -20 && time < 5){
					
									coplanarity_mass_theta_inv_collerated->Fill(phi);
									coplanarityverteilung_collerated->Fill(phi, beamphoton1E);
					
								}
									if((time > -300 && time < -100) ||(time > 100 && time < 300)){					
										coplanarity_mass_theta_inv_collerated->Fill(phi,-0.0625);
										coplanarityverteilung_collerated->Fill(phi, beamphoton1E,-0.0625);
					
									}
							}
						}
					}
					
// 			--------------------------------------------------------Missing Mass Cuts----------------------
															
		// 			__________________Polarwinkel______________________________________________________
					
					if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {
					
						if(time > -20 && time < 5){
					
							missingmass_theta_collerated->Fill(missingmass); 
						}
							if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
								
								missingmass_theta_collerated->Fill(missingmass,-0.0625); 
							}
								
					
		// 			__________________Coplanarity______________________________________________________
					
						if((phi > ( (unten_copl)) && phi < ((oben_copl)))){
					
							if(time > -20 && time < 5){
						
								missingmass_theta_copl_collerated->Fill(missingmass); 
							}
								if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
									
									missingmass_theta_copl_collerated->Fill(missingmass,-0.0625);
								}
					
		// 			__________________invariante Masse______________________________________________________
					
							if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
					
								if(time > -20 && time < 5){
					
									missingmass_theta_copl_inv_collerated->Fill(missingmass);
									missingmassverteilung_collerated_proton->Fill(missingmass, beamphoton1E); 
								}
									if((time > -300 && time < -100) ||(time > 100 && time < 300)){
									
										missingmass_theta_copl_inv_collerated->Fill(missingmass,-0.0625);
										missingmassverteilung_collerated_proton->Fill(missingmass, beamphoton1E,-0.0625);
									}								
							}
						}
					}
					
// 			--------------------------------------------------------PolarWinkel Cuts----------------------			
					
		// 			__________________Missing Mass______________________________________________________
					
					if(missingmass>(unten_mass) && missingmass<(oben_mass)){
					
						if(time > -20 && time < 5){
					
								if(rootinos->HasTAPS(0)){
								thetaproton_mass_taps_collerated->Fill(thetadiff); 
								}
								else
								{thetaproton_mass_collerated->Fill(thetadiff); 
								}	
						}
					
							if((time > -300 && time < -100) ||(time > 100 && time < 300)){							
								if(rootinos->HasTAPS(0)){
								thetaproton_mass_taps_collerated->Fill(thetadiff,-0.0625); 
								}
								else
								{thetaproton_mass_collerated->Fill(thetadiff,-0.0625); 
								}
							}		
					
		// 			__________________Coplanarity______________________________________________________
					
						if((phi > ((unten_copl)) && phi < ((oben_copl)))){
					
							if(time > -20 && time < 5){
					
								if(rootinos->HasTAPS(0)){
								thetaproton_mass_copl_taps_collerated->Fill(thetadiff); 
								}
								else
								{thetaproton_mass_copl_collerated->Fill(thetadiff); 
								}		
							}
									
								if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
									
									if(rootinos->HasTAPS(0)){
									thetaproton_mass_copl_taps_collerated->Fill(thetadiff,-0.0625); 
									}
									else
									{thetaproton_mass_copl_collerated->Fill(thetadiff,-0.0625); 
									}	
								}
					
		// 			__________________invariante Masse______________________________________________________
					
							if(((invariantemasse1 > (pionmasse + unten_inv) && invariantemasse1 < (pionmasse + oben_inv)) && (invariantemasse2 > (pionmasse + unten_inv) && invariantemasse2 < (pionmasse + oben_inv)))){
					
								if(time > -20 && time < 5){
					
									if(rootinos->HasTAPS(0)){
									thetaproton_mass_copl_inv_taps_collerated->Fill(thetadiff); 
									thetaverteilung_collerated_taps->Fill(thetadiff, beamphoton1E);
									}
									else
									{thetaproton_mass_copl_inv_collerated->Fill(thetadiff); 
									}				
									thetaverteilung_collerated->Fill(thetadiff, beamphoton1E);
								}
					
									if((time > -300 && time < -100) ||(time > 100 && time < 300)){
					
										if(rootinos->HasTAPS(0)){
										thetaproton_mass_copl_inv_taps_collerated->Fill(thetadiff,-0.0625); 
										thetaverteilung_collerated_taps->Fill(thetadiff, beamphoton1E,-0.0625);

										}
										else
										{thetaproton_mass_copl_inv_collerated->Fill(thetadiff,-0.0625); 
										}	
										thetaverteilung_collerated->Fill(thetadiff, beamphoton1E,-0.0625);
									}
							}
						}	
					}
// 			................................invarianteMasse Cut...................................................................................
										
		// 			__________________Missing Mass______________________________________________________
					
					if(missingmass>(unten_mass) && missingmass<(oben_mass)){
									
						if(time > -20 && time < 5){
							massesumme_mass_collerated_proton->Fill(invariantemasse13);
							massesumme_mass_collerated_proton->Fill(invariantemasse12);				
							massesumme_mass_collerated_proton->Fill(invariantemasse14);
							massesumme_mass_collerated_proton->Fill(invariantemasse34);
							massesumme_mass_collerated_proton->Fill(invariantemasse24);
							massesumme_mass_collerated_proton->Fill(invariantemasse23);
							massegegenmasse_mass_collerated_proton->Fill(invariantemasse12,invariantemasse34);
							massegegenmasse_mass_collerated_proton->Fill(invariantemasse13,invariantemasse24);
							massegegenmasse_mass_collerated_proton->Fill(invariantemasse14,invariantemasse23);
							massegegenmasse_mass_collerated_proton->Fill(invariantemasse34,invariantemasse12);
							massegegenmasse_mass_collerated_proton->Fill(invariantemasse24,invariantemasse13);
							massegegenmasse_mass_collerated_proton->Fill(invariantemasse23,invariantemasse14);
					}
					
							if((time > -300 && time < -100) ||(time > 100 && time < 300)){						
								massesumme_mass_collerated_proton->Fill(invariantemasse13,-0.0625);
								massesumme_mass_collerated_proton->Fill(invariantemasse12,-0.0625);				
								massesumme_mass_collerated_proton->Fill(invariantemasse14,-0.0625);
								massesumme_mass_collerated_proton->Fill(invariantemasse34,-0.0625);
								massesumme_mass_collerated_proton->Fill(invariantemasse24,-0.0625);
								massesumme_mass_collerated_proton->Fill(invariantemasse23,-0.0625);
								massegegenmasse_mass_collerated_proton->Fill(invariantemasse12,invariantemasse34,-0.0625);
								massegegenmasse_mass_collerated_proton->Fill(invariantemasse13,invariantemasse24,-0.0625);
								massegegenmasse_mass_collerated_proton->Fill(invariantemasse14,invariantemasse23,-0.0625);
								massegegenmasse_mass_collerated_proton->Fill(invariantemasse34,invariantemasse12,-0.0625);
								massegegenmasse_mass_collerated_proton->Fill(invariantemasse24,invariantemasse13,-0.0625);
								massegegenmasse_mass_collerated_proton->Fill(invariantemasse23,invariantemasse14,-0.0625);
							}
					
		// 			__________________Polarwinkel______________________________________________________
					
						if((phi > ((unten_copl)) && phi < (oben_copl))){		
									
							if(time > -20 && time < 5){
								massesumme_mass_copl_collerated->Fill(invariantemasse13);
								massesumme_mass_copl_collerated->Fill(invariantemasse12);				
								massesumme_mass_copl_collerated->Fill(invariantemasse14);
								massesumme_mass_copl_collerated->Fill(invariantemasse34);
								massesumme_mass_copl_collerated->Fill(invariantemasse24);
								massesumme_mass_copl_collerated->Fill(invariantemasse23);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse12,invariantemasse34);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse13,invariantemasse24);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse14,invariantemasse23);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse34,invariantemasse12);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse24,invariantemasse13);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse23,invariantemasse14);
							}
					
								if((time > -300 && time < -100) ||(time > 100 && time < 300)){					
								massesumme_mass_copl_collerated->Fill(invariantemasse13,-0.0625);
								massesumme_mass_copl_collerated->Fill(invariantemasse12,-0.0625);				
								massesumme_mass_copl_collerated->Fill(invariantemasse14,-0.0625);
								massesumme_mass_copl_collerated->Fill(invariantemasse34,-0.0625);
								massesumme_mass_copl_collerated->Fill(invariantemasse24,-0.0625);
								massesumme_mass_copl_collerated->Fill(invariantemasse23,-0.0625);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse12,invariantemasse34,-0.0625);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse13,invariantemasse24,-0.0625);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse14,invariantemasse23,-0.0625);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse34,invariantemasse12,-0.0625);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse24,invariantemasse13,-0.0625);
								massegegenmasse_mass_copl_collerated->Fill(invariantemasse23,invariantemasse14,-0.0625);
								}	
									
		// 			__________________Coplanarity______________________________________________________
																							
							if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {
						
								if(time > -20 && time < 5){
									massesumme_mass_theta_copl_collerated->Fill(invariantemasse13);
									massesumme_mass_theta_copl_collerated->Fill(invariantemasse12);				
									massesumme_mass_theta_copl_collerated->Fill(invariantemasse14);
									massesumme_mass_theta_copl_collerated->Fill(invariantemasse34);
									massesumme_mass_theta_copl_collerated->Fill(invariantemasse24);
									massesumme_mass_theta_copl_collerated->Fill(invariantemasse23);
					
									massesumme_mass_theta_copl_inv_collerated->Fill(invariantemasse2);
									massesumme_mass_theta_copl_inv_collerated->Fill(invariantemasse1);
					
									massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse12,invariantemasse34);
									massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse13,invariantemasse24);
									massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse14,invariantemasse23);
									massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse34,invariantemasse12);
									massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse24,invariantemasse13);
									massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse23,invariantemasse14);
								}
					
									if((time > -300 && time < -100) ||(time > 100 && time < 300)){					
										massesumme_mass_theta_copl_collerated->Fill(invariantemasse13,-0.0625);
										massesumme_mass_theta_copl_collerated->Fill(invariantemasse12,-0.0625);				
										massesumme_mass_theta_copl_collerated->Fill(invariantemasse14,-0.0625);
										massesumme_mass_theta_copl_collerated->Fill(invariantemasse34,-0.0625);
										massesumme_mass_theta_copl_collerated->Fill(invariantemasse24,-0.0625);
										massesumme_mass_theta_copl_collerated->Fill(invariantemasse23,-0.0625);
					
										massesumme_mass_theta_copl_inv_collerated->Fill(invariantemasse2,-0.0625);
										massesumme_mass_theta_copl_inv_collerated->Fill(invariantemasse1,-0.0625);
					
										massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse12,invariantemasse34,-0.0625);
										massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse13,invariantemasse24,-0.0625);
										massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse14,invariantemasse23,-0.0625);
										massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse34,invariantemasse12,-0.0625);
										massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse24,invariantemasse13,-0.0625);
										massegegenmasse_mass_theta_copl_collerated->Fill(invariantemasse23,invariantemasse14,-0.0625);
									}
					
								if(invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)){
					
									if(time > -20 && time < 5 ){
									massesumme_mass_theta_copl_inv1_collerated->Fill(invariantemasse2);
									massesumme_mass_theta_copl_inv_beam_collerated->Fill(invariantemasse2,beamphoton1E);
					
									}
									
										if((time > -300 && time < -100) ||(time > 100 && time < 300)){				
										massesumme_mass_theta_copl_inv1_collerated->Fill(invariantemasse2,-0.0625);
										massesumme_mass_theta_copl_inv_beam_collerated->Fill(invariantemasse2,beamphoton1E,-0.0625);
					
										}
													
								}	
								
							}							    				
						}		
					}

					//stuff after all cuts and creation of needed analysis histograms
					if(time > -20 && time < 5){
	
						if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((invariantemasse1 > (unten_inv_proton) && invariantemasse1 < (oben_inv_proton)))&& ((invariantemasse2 > (unten_inv_proton) && invariantemasse2 < (oben_inv_proton)))){

						if(rootinos->HasCB(0)){
							Check_CBdE_E->Fill(proton_4vect_meas.E(),rootinos->GetVetoEnergy(0));
						}

						if(rootinos->HasTAPS(0)){
							Check_TAPSdE_E->Fill(proton_4vect_meas.E(),rootinos->GetVetoEnergy(0));
						}

						cosverteilung_collerated_proton->Fill(cospi0pi0,beamphoton1E);
						cosppi0verteilung_collerated_proton->Fill(cosppi01,beamphoton1E,0.5);
						cosppi0verteilung_collerated_proton->Fill(cosppi02,beamphoton1E,0.5);
				
						invverteilung_collerated_proton->Fill(invariantmass_pi0pi0,beamphoton1E);
						invppi0verteilung_collerated_proton->Fill(invariantmass_ppi01,beamphoton1E,0.5);
						invppi0verteilung_collerated_proton->Fill(invariantmass_ppi02,beamphoton1E,0.5);
				
						cosverteilung_thetaproton_collerated_proton->Fill(cospi0pi0,anglethetaproton_meas);
									
 							if(planesetting=="PARA"){

				
								if(pt>0 || pt==5){
									kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0pi0, beamphoton1E);
									kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0pi0, beamphoton1E,pb);
									kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0pi0, beamphoton1E,pt);
									kristallminus_targetplus_ppi0_collerated_proton->Fill(phippi01,cosppi01, beamphoton1E,0.5);
									kristallminus_targetplus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01, beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetplus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01, beamphoton1E,scalsignalpt_ppi0);
									kristallminus_targetplus_ppi0_collerated_proton->Fill(phippi02,cosppi02, beamphoton1E,0.5);
									kristallminus_targetplus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02, beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetplus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02, beamphoton1E,scalsignalpt_ppi0);
					
									kristallminus_targetplus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E);
									kristallminus_targetplus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,pb);
									kristallminus_targetplus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,pt);
									kristallminus_targetplus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01, beamphoton1E,0.5);
									kristallminus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalsignalpt_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02, beamphoton1E,0.5);
									kristallminus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalsignalpt_ppi0);
								}
		
								if(pt<0 || pt==5){
									kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E);
									kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,pb);
									kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,pt);
									kristallminus_targetminus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,0.5);
									kristallminus_targetminus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetminus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpt_ppi0);
									kristallminus_targetminus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,0.5);
									kristallminus_targetminus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetminus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpt_ppi0);
					
									kristallminus_targetminus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E);
									kristallminus_targetminus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pb);
									kristallminus_targetminus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pt);
									kristallminus_targetminus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,0.5);
									kristallminus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpt_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,0.5);
									kristallminus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpb_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpt_ppi0);
								}
							}

 							if(planesetting=="PERP"){

			
								if(pt>0 || pt==5){
									kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E);
									kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,pb);
									kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,pt);
									kristallplus_targetplus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,0.5);
									kristallplus_targetplus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetplus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpt_ppi0);
									kristallplus_targetplus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,0.5);
									kristallplus_targetplus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetplus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpt_ppi0);
					
									kristallplus_targetplus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E);
									kristallplus_targetplus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pb);
									kristallplus_targetplus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pt);
									kristallplus_targetplus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,0.5);
									kristallplus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpt_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,0.5);
									kristallplus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpt_ppi0);
								}
				
								if(pt<0 || pt==5){
									kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E);
									kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,pb);
									kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,pt);
									kristallplus_targetminus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,0.5);
									kristallplus_targetminus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetminus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpt_ppi0);
									kristallplus_targetminus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,0.5);
									kristallplus_targetminus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetminus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpt_ppi0);
					
									kristallplus_targetminus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E);
									kristallplus_targetminus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pb);
									kristallplus_targetminus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pt);
									kristallplus_targetminus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,0.5);
									kristallplus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpt_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,0.5);
									kristallplus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpb_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpt_ppi0);
								}
	
							}						
	
						}
	
					}//end prompt
	
					if((time > -300 && time < -100) ||(time > 100 && time < 300)){	
	
						if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((invariantemasse1 > (unten_inv_proton) && invariantemasse1 < (oben_inv_proton)))&& ((invariantemasse2 > (unten_inv_proton) && invariantemasse2 < (oben_inv_proton)))){

							if(rootinos->HasCB(0)){
								Check_CBdE_E->Fill(proton_4vect_meas.E(),rootinos->GetVetoEnergy(0),-0.0625);
							}
	
							if(rootinos->HasTAPS(0)){
								Check_TAPSdE_E->Fill(proton_4vect_meas.E(),rootinos->GetVetoEnergy(0),-0.0625);
							}
		
							cosverteilung_collerated_proton->Fill(cospi0pi0,beamphoton1E,-0.0625);
							cosppi0verteilung_collerated_proton->Fill(cosppi01,beamphoton1E,scalbackground_ppi0);
							cosppi0verteilung_collerated_proton->Fill(cosppi02,beamphoton1E,scalbackground_ppi0);
				
							invverteilung_collerated_proton->Fill(invariantmass_pi0pi0,beamphoton1E,-0.0625);
							invppi0verteilung_collerated_proton->Fill(invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
							invppi0verteilung_collerated_proton->Fill(invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
				
							cosverteilung_thetaproton_collerated_proton->Fill(cospi0pi0,anglethetaproton_meas,-0.0625);

							if(planesetting=="PARA"){

			
								if(pt>0 || pt==5){
									kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0pi0, beamphoton1E,-0.0625);
									kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0pi0, beamphoton1E,scalbackgroundpb);
									kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0pi0, beamphoton1E,scalbackgroundpt);
									kristallminus_targetplus_ppi0_collerated_proton->Fill(phippi01,cosppi01, beamphoton1E,scalbackground_ppi0);
									kristallminus_targetplus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01, beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetplus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01, beamphoton1E,scalbackgroundpt_ppi0);
									kristallminus_targetplus_ppi0_collerated_proton->Fill(phippi02,cosppi02, beamphoton1E,scalbackground_ppi0);
									kristallminus_targetplus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02, beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetplus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02, beamphoton1E,scalbackgroundpt_ppi0);
					
									kristallminus_targetplus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,-0.0625);
									kristallminus_targetplus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,scalbackgroundpb);
									kristallminus_targetplus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,scalbackgroundpt);
									kristallminus_targetplus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalbackground_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalbackgroundpt_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalbackground_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalbackgroundpt_ppi0);
								}
			
								if(pt<0 || pt==5){
									kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E,-0.0625);
									kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpb);
									kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpt);
									kristallminus_targetminus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackground_ppi0);
									kristallminus_targetminus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetminus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpt_ppi0);
									kristallminus_targetminus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackground_ppi0);
									kristallminus_targetminus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetminus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpt_ppi0);
					
									kristallminus_targetminus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,-0.0625);
									kristallminus_targetminus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpb);
									kristallminus_targetminus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpt);
									kristallminus_targetminus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01,scalbackgroundpt_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpb_ppi0);
									kristallminus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpt_ppi0);
								}
							}

							if(planesetting=="PERP"){

			
								if(pt>0 || pt==5){
									kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E,-0.0625);
									kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpb);
									kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpt);
									kristallplus_targetplus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetplus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetplus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpt_ppi0);
									kristallplus_targetplus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetplus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetplus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpt_ppi0);
					
									kristallplus_targetplus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,-0.0625);
									kristallplus_targetplus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpb);
									kristallplus_targetplus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpt);
									kristallplus_targetplus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpt_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetplus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpt_ppi0);
								}
				
								if(pt<0 || pt==5){
									kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0pi0,beamphoton1E,-0.0625);
									kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpb);
									kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpt);
									kristallplus_targetminus_ppi0_collerated_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetminus_ppi0_collerated_pb_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetminus_ppi0_collerated_pt_proton->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpt_ppi0);
									kristallplus_targetminus_ppi0_collerated_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetminus_ppi0_collerated_pb_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetminus_ppi0_collerated_pt_proton->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpt_ppi0);
					
									kristallplus_targetminus_im_collerated_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,-0.0625);
									kristallplus_targetminus_im_collerated_pb_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpb);
									kristallplus_targetminus_im_collerated_pt_proton->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpt);
									kristallplus_targetminus_im_ppi0_collerated_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpt_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_pb_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpb_ppi0);
									kristallplus_targetminus_im_ppi0_collerated_pt_proton->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpt_ppi0);
								}
							}						
						}
	
					}//end sideband
// 				
				}
//_________________________________________________ end of all events with proton identified

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)))&& ((invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))&&trigger->GetNErrors()!=0){
				events_witherror->Fill(phimeson,cospi0pi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)))&& ((invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
				events_all->Fill(phimeson,cospi0pi0, beamphoton1E);
			}

			if(time > -20 && time < 5){

				if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)))&& ((invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
					cosverteilung_collerated->Fill(cospi0pi0,beamphoton1E);
					cosppi0verteilung_collerated->Fill(cosppi01,beamphoton1E,0.5);
					cosppi0verteilung_collerated->Fill(cosppi02,beamphoton1E,0.5);
			
					invverteilung_collerated->Fill(invariantmass_pi0pi0,beamphoton1E);
					invppi0verteilung_collerated->Fill(invariantmass_ppi01,beamphoton1E,0.5);
					invppi0verteilung_collerated->Fill(invariantmass_ppi02,beamphoton1E,0.5);
			
 					if(planesetting=="PARA"){

		
						if(pt>0 || pt==5){
							kristallminus_targetplus_collerated->Fill(phimeson,cospi0pi0, beamphoton1E);
							kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0pi0, beamphoton1E,pb);
							kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0pi0, beamphoton1E,pt);
							kristallminus_targetplus_ppi0_collerated->Fill(phippi01,cosppi01, beamphoton1E,0.5);
							kristallminus_targetplus_ppi0_collerated_pb->Fill(phippi01,cosppi01, beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetplus_ppi0_collerated_pt->Fill(phippi01,cosppi01, beamphoton1E,scalsignalpt_ppi0);
							kristallminus_targetplus_ppi0_collerated->Fill(phippi02,cosppi02, beamphoton1E,0.5);
							kristallminus_targetplus_ppi0_collerated_pb->Fill(phippi02,cosppi02, beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetplus_ppi0_collerated_pt->Fill(phippi02,cosppi02, beamphoton1E,scalsignalpt_ppi0);
			
							kristallminus_targetplus_im_collerated->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E);
							kristallminus_targetplus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,pb);
							kristallminus_targetplus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,pt);
							kristallminus_targetplus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01, beamphoton1E,0.5);
							kristallminus_targetplus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetplus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalsignalpt_ppi0);
							kristallminus_targetplus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02, beamphoton1E,0.5);
							kristallminus_targetplus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetplus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalsignalpt_ppi0);
						}

						if(pt<0 || pt==5){
							kristallminus_targetminus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E);
							kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,pb);
							kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,pt);
							kristallminus_targetminus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,0.5);
							kristallminus_targetminus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetminus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpt_ppi0);
							kristallminus_targetminus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,0.5);
							kristallminus_targetminus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetminus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpt_ppi0);
			
							kristallminus_targetminus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E);
							kristallminus_targetminus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pb);
							kristallminus_targetminus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pt);
							kristallminus_targetminus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01,beamphoton1E,0.5);
							kristallminus_targetminus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetminus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpt_ppi0);
							kristallminus_targetminus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,0.5);
							kristallminus_targetminus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpb_ppi0);
							kristallminus_targetminus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpt_ppi0);
						}
					}

					if(planesetting=="PERP"){
		
						if(pt>0 || pt==5){
							kristallplus_targetplus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E);
							kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,pb);
							kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,pt);
							kristallplus_targetplus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,0.5);
							kristallplus_targetplus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetplus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpt_ppi0);
							kristallplus_targetplus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,0.5);
							kristallplus_targetplus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetplus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpt_ppi0);
			
							kristallplus_targetplus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E);
							kristallplus_targetplus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pb);
							kristallplus_targetplus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pt);
							kristallplus_targetplus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01,beamphoton1E,0.5);
							kristallplus_targetplus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetplus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpt_ppi0);
							kristallplus_targetplus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,0.5);
							kristallplus_targetplus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetplus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpt_ppi0);
						}
			
						if(pt<0 || pt==5){
							kristallplus_targetminus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E);
							kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,pb);
							kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,pt);
							kristallplus_targetminus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,0.5);
							kristallplus_targetminus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetminus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,scalsignalpt_ppi0);
							kristallplus_targetminus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,0.5);
							kristallplus_targetminus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetminus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,scalsignalpt_ppi0);
			
							kristallplus_targetminus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E);
							kristallplus_targetminus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pb);
							kristallplus_targetminus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,pt);
							kristallplus_targetminus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01,beamphoton1E,0.5);
							kristallplus_targetminus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetminus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalsignalpt_ppi0);
							kristallplus_targetminus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,0.5);
							kristallplus_targetminus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpb_ppi0);
							kristallplus_targetminus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalsignalpt_ppi0);
						}
					}						
				}
			}//end prompt

			if((time > -300 && time < -100) ||(time > 100 && time < 300)){	

				if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)))&& ((invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){

					cosverteilung_collerated->Fill(cospi0pi0,beamphoton1E,-0.0625);
					cosppi0verteilung_collerated->Fill(cosppi01,beamphoton1E,scalbackground_ppi0);
					cosppi0verteilung_collerated->Fill(cosppi02,beamphoton1E,scalbackground_ppi0);
		
					invverteilung_collerated->Fill(invariantmass_pi0pi0,beamphoton1E,-0.0625);
					invppi0verteilung_collerated->Fill(invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
					invppi0verteilung_collerated->Fill(invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
		
 					if(planesetting=="PARA"){

		
						if(pt>0  || pt==5){
							kristallminus_targetplus_collerated->Fill(phimeson,cospi0pi0, beamphoton1E,-0.0625);
							kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0pi0, beamphoton1E,scalbackgroundpb);
							kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0pi0, beamphoton1E,scalbackgroundpt);
							kristallminus_targetplus_ppi0_collerated->Fill(phippi01,cosppi01, beamphoton1E,scalbackground_ppi0);
							kristallminus_targetplus_ppi0_collerated_pb->Fill(phippi01,cosppi01, beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetplus_ppi0_collerated_pt->Fill(phippi01,cosppi01, beamphoton1E,scalbackgroundpt_ppi0);
							kristallminus_targetplus_ppi0_collerated->Fill(phippi02,cosppi02, beamphoton1E,scalbackground_ppi0);
							kristallminus_targetplus_ppi0_collerated_pb->Fill(phippi02,cosppi02, beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetplus_ppi0_collerated_pt->Fill(phippi02,cosppi02, beamphoton1E,scalbackgroundpt_ppi0);
			
							kristallminus_targetplus_im_collerated->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,-0.0625);
							kristallminus_targetplus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,scalbackgroundpb);
							kristallminus_targetplus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0, beamphoton1E,scalbackgroundpt);
							kristallminus_targetplus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalbackground_ppi0);
							kristallminus_targetplus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetplus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01, beamphoton1E,scalbackgroundpt_ppi0);
							kristallminus_targetplus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalbackground_ppi0);
							kristallminus_targetplus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetplus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02, beamphoton1E,scalbackgroundpt_ppi0);
						}
			
						if(pt<0  || pt==5){
							kristallminus_targetminus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E,-0.0625);
							kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpb);
							kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpt);
							kristallminus_targetminus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,scalbackground_ppi0);
							kristallminus_targetminus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetminus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpt_ppi0);
							kristallminus_targetminus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,scalbackground_ppi0);
							kristallminus_targetminus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetminus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpt_ppi0);
			
							kristallminus_targetminus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,-0.0625);
							kristallminus_targetminus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpb);
							kristallminus_targetminus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpt);
							kristallminus_targetminus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
							kristallminus_targetminus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetminus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01,scalbackgroundpt_ppi0);
							kristallminus_targetminus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
							kristallminus_targetminus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpb_ppi0);
							kristallminus_targetminus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpt_ppi0);
						}
					}

					if(planesetting=="PERP"){

		
						if(pt>0  || pt==5){
							kristallplus_targetplus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E,-0.0625);
							kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpb);
							kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpt);
							kristallplus_targetplus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetplus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetplus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpt_ppi0);
							kristallplus_targetplus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetplus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetplus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpt_ppi0);
			
							kristallplus_targetplus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,-0.0625);
							kristallplus_targetplus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpb);
							kristallplus_targetplus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpt);
							kristallplus_targetplus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetplus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetplus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpt_ppi0);
							kristallplus_targetplus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetplus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetplus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpt_ppi0);
						}

						if(pt<0  || pt==5){
							 kristallplus_targetminus_collerated->Fill(phimeson,cospi0pi0,beamphoton1E,-0.0625);
							kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpb);
							kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0pi0,beamphoton1E,scalbackgroundpt);
							kristallplus_targetminus_ppi0_collerated->Fill(phippi01,cosppi01,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetminus_ppi0_collerated_pb->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetminus_ppi0_collerated_pt->Fill(phippi01,cosppi01,beamphoton1E,scalbackgroundpt_ppi0);
							kristallplus_targetminus_ppi0_collerated->Fill(phippi02,cosppi02,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetminus_ppi0_collerated_pb->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetminus_ppi0_collerated_pt->Fill(phippi02,cosppi02,beamphoton1E,scalbackgroundpt_ppi0);
			
							kristallplus_targetminus_im_collerated->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,-0.0625);
							kristallplus_targetminus_im_collerated_pb->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpb);
							kristallplus_targetminus_im_collerated_pt->Fill(phimeson,invariantmass_pi0pi0,beamphoton1E,scalbackgroundpt);
							kristallplus_targetminus_im_ppi0_collerated->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetminus_im_ppi0_collerated_pb->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetminus_im_ppi0_collerated_pt->Fill(phippi01,invariantmass_ppi01,beamphoton1E,scalbackgroundpt_ppi0);
							kristallplus_targetminus_im_ppi0_collerated->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackground_ppi0);
							kristallplus_targetminus_im_ppi0_collerated_pb->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpb_ppi0);
							kristallplus_targetminus_im_ppi0_collerated_pt->Fill(phippi02,invariantmass_ppi02,beamphoton1E,scalbackgroundpt_ppi0);
						}
					}						
				}
			}//end sideband

		}//tagger
}


void P2Pi0Analyse::fOnEndProcessing() {

}

Bool_t	P2Pi0Analyse::Write(){
//FOR CARBON
// polsetting=P2Pi0Analyse::poledge(inputFile);
// TString *data=new TString(polsetting);
// TDirectory* curDir1  = outputFile->mkdir(Form("%s",data->Data()));

TDirectory* curDir1  = outputFile->mkdir(Form("Butanol%i",(Int_t)linpol->GetEdgeSetting()));
curDir1->cd();
poltable_energy->Write();
poltable_energy_weight->Write();
test->Write();
TDirectory* curDir3  = curDir1->mkdir("Selektion");
curDir3->cd();
time_prompt->Write();
time1->Write();
time_side->Write();
missingmass_collerated->Write();
missingmass_inv_collerated->Write();
massesumme_collerated->Write();
massegegenmasse_mass_collerated->Write();
massesumme_mass_collerated->Write();


curDir1->cd();
TDirectory* curDir5  = curDir1->mkdir("Selektion_withProton");
curDir5->cd();
coplanarity_collerated->Write();
coplanarity_mass_theta_collerated->Write();
coplanarity_mass_theta_inv_collerated->Write();
coplanarity_mass_collerated->Write();
coplanarity_theta_collerated->Write();
coplanarity_inv_collerated->Write();
massegegenmasse_mass_collerated_proton->Write();
massesumme_mass_collerated_proton->Write();
// 
missingmass_theta_copl_collerated->Write();
missingmass_theta_copl_inv_collerated->Write();
missingmass_theta_collerated->Write();
missingmass_copl_collerated->Write();
missingmass_inv_collerated->Write();
// 
massesumme_mass_theta_copl_inv_beam_collerated->Write();
massesumme_collerated->Write();
massesumme_mass_theta_copl_collerated->Write();
massesumme_mass_theta_copl_inv1_collerated->Write();
massesumme_mass_theta_copl_inv_collerated->Write();
massesumme_mass_copl_collerated->Write();
massesumme_copl_collerated->Write();
massesumme_theta_collerated->Write();
massegegenmasse_collerated->Write();
massegegenmasse_mass_copl_collerated->Write();
massegegenmasse_mass_theta_copl_collerated->Write();
massegegenmasse_copl_collerated->Write();
massegegenmasse_theta_collerated->Write();
massegegenmasse_mass_collerated->Write();
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
Check_CBdE_E->Write();
Check_TAPSdE_E->Write();


curDir1->cd();
TDirectory* curDir2  = curDir1->mkdir("Allgemein");
curDir2->cd();
cosverteilung_collerated->Write();
cosppi0verteilung_collerated->Write(); 
invverteilung_collerated->Write();
invppi0verteilung_collerated->Write();

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

kristallminus_targetplus_ppi0_collerated->Write();
kristallminus_targetminus_ppi0_collerated->Write();
kristallplus_targetplus_ppi0_collerated->Write();
kristallplus_targetminus_ppi0_collerated->Write();

kristallminus_targetplus_ppi0_collerated_pb->Write();
kristallminus_targetminus_ppi0_collerated_pb->Write();
kristallplus_targetplus_ppi0_collerated_pb->Write();
kristallplus_targetminus_ppi0_collerated_pb->Write();

kristallminus_targetplus_ppi0_collerated_pt->Write();
kristallminus_targetminus_ppi0_collerated_pt->Write();
kristallplus_targetplus_ppi0_collerated_pt->Write();
kristallplus_targetminus_ppi0_collerated_pt->Write();

kristallminus_targetplus_im_collerated->Write();
kristallminus_targetminus_im_collerated->Write();
kristallplus_targetplus_im_collerated->Write();
kristallplus_targetminus_im_collerated->Write();

kristallminus_targetplus_im_collerated_pb->Write();
kristallminus_targetminus_im_collerated_pb->Write();
kristallplus_targetplus_im_collerated_pb->Write();
kristallplus_targetminus_im_collerated_pb->Write();

kristallminus_targetplus_im_collerated_pt->Write();
kristallminus_targetminus_im_collerated_pt->Write();
kristallplus_targetplus_im_collerated_pt->Write();
kristallplus_targetminus_im_collerated_pt->Write();

kristallminus_targetplus_im_ppi0_collerated->Write();
kristallminus_targetminus_im_ppi0_collerated->Write();
kristallplus_targetplus_im_ppi0_collerated->Write();
kristallplus_targetminus_im_ppi0_collerated->Write();

kristallminus_targetplus_im_ppi0_collerated_pb->Write();
kristallminus_targetminus_im_ppi0_collerated_pb->Write();
kristallplus_targetplus_im_ppi0_collerated_pb->Write();
kristallplus_targetminus_im_ppi0_collerated_pb->Write();

kristallminus_targetplus_im_ppi0_collerated_pt->Write();
kristallminus_targetminus_im_ppi0_collerated_pt->Write();
kristallplus_targetplus_im_ppi0_collerated_pt->Write();
kristallplus_targetminus_im_ppi0_collerated_pt->Write();

missingmassverteilung_collerated->Write();
events_witherror->Write();
events_all->Write();

curDir1->cd();
TDirectory* curDir4  = curDir1->mkdir("Allgemein_withProton");
curDir4->cd();
cosverteilung_collerated_proton->Write();
cosppi0verteilung_collerated_proton->Write(); 
invverteilung_collerated_proton->Write();
invppi0verteilung_collerated_proton->Write();

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

kristallminus_targetplus_ppi0_collerated->Write();
kristallminus_targetminus_ppi0_collerated->Write();
kristallplus_targetplus_ppi0_collerated->Write();
kristallplus_targetminus_ppi0_collerated->Write();

kristallminus_targetplus_ppi0_collerated_pb->Write();
kristallminus_targetminus_ppi0_collerated_pb->Write();
kristallplus_targetplus_ppi0_collerated_pb->Write();
kristallplus_targetminus_ppi0_collerated_pb->Write();

kristallminus_targetplus_ppi0_collerated_pt->Write();
kristallminus_targetminus_ppi0_collerated_pt->Write();
kristallplus_targetplus_ppi0_collerated_pt->Write();
kristallplus_targetminus_ppi0_collerated_pt->Write();

kristallminus_targetplus_im_collerated_proton->Write();
kristallminus_targetminus_im_collerated_proton->Write();
kristallplus_targetplus_im_collerated_proton->Write();
kristallplus_targetminus_im_collerated_proton->Write();

kristallminus_targetplus_im_collerated_pb_proton->Write();
kristallminus_targetminus_im_collerated_pb_proton->Write();
kristallplus_targetplus_im_collerated_pb_proton->Write();
kristallplus_targetminus_im_collerated_pb_proton->Write();

kristallminus_targetplus_im_collerated_pt_proton->Write();
kristallminus_targetminus_im_collerated_pt_proton->Write();
kristallplus_targetplus_im_collerated_pt_proton->Write();
kristallplus_targetminus_im_collerated_pt_proton->Write();

kristallminus_targetplus_im_ppi0_collerated_proton->Write();
kristallminus_targetminus_im_ppi0_collerated_proton->Write();
kristallplus_targetplus_im_ppi0_collerated_proton->Write();
kristallplus_targetminus_im_ppi0_collerated_proton->Write();

kristallminus_targetplus_im_ppi0_collerated_pb_proton->Write();
kristallminus_targetminus_im_ppi0_collerated_pb_proton->Write();
kristallplus_targetplus_im_ppi0_collerated_pb_proton->Write();
kristallplus_targetminus_im_ppi0_collerated_pb_proton->Write();

kristallminus_targetplus_im_ppi0_collerated_pt_proton->Write();
kristallminus_targetminus_im_ppi0_collerated_pt_proton->Write();
kristallplus_targetplus_im_ppi0_collerated_pt_proton->Write();
kristallplus_targetminus_im_ppi0_collerated_pt_proton->Write();

thetaverteilung_collerated->Write();
thetaverteilung_collerated_taps->Write();
coplanarityverteilung_collerated->Write();
missingmassverteilung_collerated_proton->Write();

outputFile->Close();
}


void	P2Pi0Analyse::ProcessScalerRead()
{
    //time.ScalerReadCorrection(5);
}
