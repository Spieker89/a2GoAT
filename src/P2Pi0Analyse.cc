#include "P2Pi0Analyse.h"

P2Pi0Analyse::P2Pi0Analyse()
{ 
fitter.AddConstraintsTotEnergy();
fitter.AddConstraintsTotMomentum();
fitter.AddConstraintsIM();
fitter.AddConstraintsMM();
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
missingmass_collerated_proton = new TH1F("missingmass_collerated", "Missing Mass (Signal);m_{mm} [MeV]", 500,0,2200);

missingmass_inv_collerated = new TH1F("missingmass_inv_collerated", "Missing Mass (Cut: invariante Masse(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_inv_collerated_proton = new TH1F("missingmass_inv_collerated_proton", "Missing Mass (Cut: invariante Masse(Signal));m_{mm} [MeV]", 500,0,2200);

missingmass_theta_copl_collerated = new TH1F("missingmass_theta_copl_collerated", "Missing Mass (Cut: Polar+Azimut-Winkel(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_theta_copl_inv_collerated = new TH1F("missingmass_theta_copl_inv_collerated", "Missing Mass (Cut: All(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_theta_collerated = new TH1F("missingmass_theta_collerated", "Missing Mass (Cut: Polarwinkel(Signal));m_{mm} [MeV]", 500,0,2200);
missingmass_copl_collerated = new TH1F("missingmass_copl_collerated", "Missing Mass (Cut: Azimutalwinkel(Signal));m_{mm} [MeV]", 500,0,2200);

//........................Histogramme für invariante Masse.....................................
massesumme_mass_theta_copl_inv_beam_collerated = new TH2F("massesumme_mass_theta_copl_inv_beam_collerated","inv Mass in dependency of E_{beam}; m_{#gamma #gamma} [MeV]; E^{rec}_{#gamma} [MeV]",400,0,600,200,200,800);
massesumme_mass_theta_copl_inv_beam_collerated->Sumw2();
massesumme_collerated = new TH1F("massesumme_collerated", "Invariante Masse von der Summe (Signal); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_collerated_proton = new TH1F("massesumme_collerated_proton", "Invariante Masse von der Summe (Signal); m_{#gamma #gamma} [MeV]", 400,0,600);

massesumme_mass_theta_copl_collerated = new TH1F("massesumme_mass_theta_copl_collerated", "Invariante Masse von der Summe (Cut: All(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_theta_copl_inv1_collerated = new TH1F("massesumme_mass_theta_copl_inv1_collerated", "Invariante Masse von der Summe (Cut: All+Cut auf invariante Masse 1(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_theta_copl_inv_collerated = new TH1F("massesumme_mass_theta_copl_inv_collerated", "Invariante Masse von der Summe (Cut: All+Cut auf invariante Masse 1(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_copl_collerated = new TH1F("massesumme_mass_copl_collerated", "Invariante Masse von der Summe (Cut: Azimutwinkel+Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_collerated = new TH1F("massesumme_mass_collerated", "Invariante Masse von der Summe (Cut: Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_copl_collerated = new TH1F("massesumme_copl_collerated", "Invariante Masse von der Summe (Cut: Azimutwinkel(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_theta_collerated = new TH1F("massesumme_theta_collerated", "Invariante Masse von der Summe (Cut: Polarwinkel(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);
massesumme_mass_collerated_proton = new TH1F("massesumme_mass_collerated_proton", "Invariante Masse von der Summe (Cut: Missing Mass(Signal)); m_{#gamma #gamma} [MeV]", 400,0,600);


massegegenmasse_collerated = new TH2F("massegegenmasse_collerated","invariante Massen gegeneinander (Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

massegegenmasse_collerated_proton = new TH2F("massegegenmasse_collerated_proton","invariante Massen gegeneinander (Signal); m_{#gamma 1  #gamma 2} [MeV]; m_{#gamma 3  #gamma 4 } [MeV]",200,0,350,200,0,350);

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

h_energy_sum_2pi0_5ped = new TH1F("h_energy_sum_2pi0_5ped ", "h_energy_sum_2pi0_5ped ",  300, 0, 1557);
h_energy_sum_2pi0_5ped->Sumw2();

h_energy_taps_2pi0_5ped = new TH2F("h_energy_taps_2pi0_5ped ", "h_energy_taps_2pi0_5ped ",  900, 0, 900,300,0,1557);
h_energy_taps_2pi0_5ped->Sumw2();



h_energy_sum_2pi0_5ped_weightcos2pi0 = new TH1F("h_energy_sum_2pi0_5ped_weightcos2pi0 ", "h_energy_sum_2pi0_5ped_weightcos2pi0 ",  300, 0, 1557);
h_energy_sum_2pi0_5ped_weightcos2pi0->Sumw2();

h_energy_sum_2pi0_5ped_weightcosppi0 = new TH1F("h_energy_sum_2pi0_5ped_weightcosppi0 ", "h_energy_sum_2pi0_5ped_weightcosppi0 ",  300, 0, 1557);
h_energy_sum_2pi0_5ped_weightcosppi0->Sumw2();


//generated
cos2pi0_beamphoton_energysum_monte = new TH3F("cos2pi0_beamphoton_energysum_monte", "cos2pi0_beamphoton_energysum_monte;cos #theta_{2#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);
cos2pi0_beamphoton_energysum_monte->Sumw2();

cosppi0_beamphoton_energysum_monte = new TH3F("cosppi0_beamphoton_energysum_monte", "cosppi0_beamphoton_energysum_monte;cos #theta_{p#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);
cosppi0_beamphoton_energysum_monte->Sumw2();

m2pi0_beamphoton_energysum_monte = new TH3F("m2pi0_beamphoton_energysum_monte", "m2pi0_beamphoton_energysum_monte;m_{2#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  70, 0, 700,30, 220, 1420 ,200, 0, 2000);
m2pi0_beamphoton_energysum_monte->Sumw2();

mppi0_beamphoton_energysum_monte = new TH3F("mppi0_beamphoton_energysum_monte", "mppi0_beamphoton_energysum_monte;m_{p#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  150, 0, 1500,30, 220, 1420 ,200, 0, 2000);
mppi0_beamphoton_energysum_monte->Sumw2();

//reconstructed
cos2pi0_beamphoton_energysum_rek = new TH3F("cos2pi0_beamphoton_energysum_rek", "cos2pi0_beamphoton_energysum_rek;cos #theta_{2#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);
cos2pi0_beamphoton_energysum_rek->Sumw2();

cosppi0_beamphoton_energysum_rek = new TH3F("cosppi0_beamphoton_energysum_rek", "cosppi0_beamphoton_energysum_rek;cos #theta_{p#pi};E_{#gamma} [MeV];E_{sum} [MeV]",  18, -1, 1,30, 220, 1420 ,200, 0, 2000);
cosppi0_beamphoton_energysum_rek->Sumw2();

m2pi0_beamphoton_energysum_rek = new TH3F("m2pi0_beamphoton_energysum_rek", "m2pi0_beamphoton_energysum_rek;m_{2#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  70, 0, 700,30, 220, 1420 ,200, 0, 2000);
m2pi0_beamphoton_energysum_rek->Sumw2();

mppi0_beamphoton_energysum_rek = new TH3F("mppi0_beamphoton_energysum_rek", "mppi0_beamphoton_energysum_rek;m_{p#pi} [MeV];E_{#gamma} [MeV];E_{sum} [MeV]",  150, 0, 1500,30, 220, 1420 ,200, 0, 2000);
mppi0_beamphoton_energysum_rek->Sumw2();


h_bph_e_gen  = new TH1F("h_bph_e_gen ", "h_bph_e_gen ",  300, 0, 1557);

triggertest = new TH1F("triggertest","triggerspattern",34,0,34);


Check_CBdE_E= new TH2F("Check_CBdE_E", "dE_E (all CB clusters compared to PID hits)", 	400, 0, 400, 100, 0, 10);
Check_TAPSdE_E= new TH2F("Check_TAPSdE_E", "dE_E (all TAPS clusters compared to Veto hits)", 	400, 0, 400, 100, 0, 10);

HistoManu::InitCutss(-8, 8,100, 300);

}

P2Pi0Analyse::~P2Pi0Analyse()
{
}



Bool_t	P2Pi0Analyse::Start()
{
pt=P2Pi0Analyse::targetpol(inputFile);

// planesetting=P2Pi0Analyse::polplane(inputFile);
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

	//Pi0Pi0-System
	TLorentzVector pi0pi0_4vektor_mc=meson1_mc+meson2_mc;
	TLorentzVector pi0pi0_4vektor_mc_boost = CMVector(pi0pi0_4vektor_mc, beam_mc, proton_mc);
	Double_t invariantmass_pi0pi0_mc = pi0pi0_4vektor_mc_boost.M();
	TVector3 pi0pi0_3vektor_boost_mc = pi0pi0_4vektor_mc_boost.Vect();
	Double_t cospi0pi0_mc = pi0pi0_3vektor_boost_mc.CosTheta();
	Double_t phimeson_mc = TMath::RadToDeg()*(pi0pi0_4vektor_mc_boost.Vect().Phi());
			
	//pPi0-System for first Pion
	TLorentzVector ppi01_4vektor_mc = proton_mc+meson1_mc;
	TLorentzVector ppi01_4vektor_boost_mc = CMVector(ppi01_4vektor_mc, beam_mc, proton_mc);
	Double_t invariantmass_ppi01_mc = ppi01_4vektor_boost_mc.M();
	TVector3 ppi01_3vektor_boost_mc = ppi01_4vektor_boost_mc.Vect();
	Double_t cosppi01_mc = ppi01_3vektor_boost_mc.CosTheta();
	Double_t phippi01_mc = TMath::RadToDeg()*ppi01_3vektor_boost_mc.Phi();		
	//pPi0-System for second Pion
	TLorentzVector ppi02_4vektor_mc = proton_mc+meson2_mc;
	TLorentzVector ppi02_4vektor_boost_mc = CMVector(ppi02_4vektor_mc, beam_mc, proton_mc);
	Double_t invariantmass_ppi02_mc = ppi02_4vektor_boost_mc.M();
	TVector3 ppi02_3vektor_boost_mc = ppi02_4vektor_boost_mc.Vect();
	Double_t cosppi02_mc = ppi02_3vektor_boost_mc.CosTheta();
	Double_t phippi02_mc = TMath::RadToDeg()*ppi02_3vektor_boost_mc.Phi();

	cos2pi0_beamphoton_energysum_monte->Fill(cospi0pi0_mc,e_beam_mc,1000*GetGeant()->GetCBESum());
	m2pi0_beamphoton_energysum_monte->Fill(invariantmass_pi0pi0_mc,e_beam_mc,1000*GetGeant()->GetCBESum());

	cosppi0_beamphoton_energysum_monte->Fill(cosppi01_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
	mppi0_beamphoton_energysum_monte->Fill(invariantmass_ppi01_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
	cosppi0_beamphoton_energysum_monte->Fill(cosppi02_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
	mppi0_beamphoton_energysum_monte->Fill(invariantmass_ppi02_mc,e_beam_mc,1000*GetGeant()->GetCBESum(),0.5);
 }


if(GetPhotons()->GetNParticles()==4){
//:::::::::::::::::::::::::::::::::::::::::::::::::::Photonen festlegen::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
			
	TLorentzVector photon1_4vektor = GetPhotons()->Particle(0);
	TLorentzVector photon2_4vektor = GetPhotons()->Particle(1);
	TLorentzVector photon3_4vektor = GetPhotons()->Particle(2);
	TLorentzVector photon4_4vektor = GetPhotons()->Particle(3);

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

// //_______________________________passende Kombination finden____________________________________________	
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
// 
	Double_t Chi12 = ChiPionPion(invariantemasse12,photon1_4vektor,photon2_4vektor,invariantemasse34,photon3_4vektor,photon4_4vektor);
	Double_t Chi13 = ChiPionPion(invariantemasse13,photon1_4vektor,photon3_4vektor,invariantemasse24,photon2_4vektor,photon4_4vektor);
	Double_t Chi14 = ChiPionPion(invariantemasse14,photon1_4vektor,photon4_4vektor,invariantemasse23,photon2_4vektor,photon3_4vektor);
// 
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
// 
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



		for(Int_t j=0; j < GetTagger()->GetNTagged();j++)
		{
// 		if(Chi > vergleich){continue;}

			//reject the pbwo4
			if(GetPhotons()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;
			if(GetPhotons()->GetDetectors(1)==GTreeTrack::DETECTOR_PbWO4) continue;
			if(GetRootinos()->GetDetectors(0)==GTreeTrack::DETECTOR_PbWO4) continue;

// 			if(pt==5){continue;} //for carbon measurements
			//if(pt==0){cotinue;}
			//if(GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j))==-1){continue;}//skip events with wrong polarisation!

			Double_t time= GetTagger()->GetTaggedTime(j) - 0.25*(GetPhotons()->GetTime(1)+GetPhotons()->GetTime(0)+GetPhotons()->GetTime(2)+GetPhotons()->GetTime(3));
			time1->Fill(time);



 			if(HistoManu::IsPromptt(time)){
				time_prompt->Fill(time);
			}

				if(HistoManu::IsRandomm(time)){	
					time_side->Fill(time);
				}
			poltable_energy->Fill(GetTagger()->GetTaggedEnergy(j));
			poltable_energy_weight->Fill(GetTagger()->GetTaggedEnergy(j),GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j)));

			//get target, beam and missng particle information
			TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
			TLorentzVector  beam_4vect = TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(j),GetTagger()->GetTaggedEnergy(j));
			TLorentzVector  missingp_4vect = beam_4vect + protonvektor_target - pi02_4vect-pi01_4vect;
			Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
			Double_t beamphoton1E=GetTagger()->GetTaggedEnergy(j);
			Double_t pb=GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(j));

			//due to kinematic not possible -> kick them out
			if(anglethetaproton_rek > 90){continue;}
			
 			if(beamphoton1E<200 || beamphoton1E>800){continue;}

			//energy dependent cuts
			Double_t oben_copl=180+15;
			Double_t unten_copl=180-15;
			Double_t oben_theta=12;
			Double_t unten_theta=-12;
			Double_t oben_mass=938+67;
			Double_t unten_mass=938-67;
			Double_t oben_mass_proton=938+67;
			Double_t unten_mass_proton=938-67;
			
			Double_t oben_inv=135+20;
			Double_t unten_inv=135-20;
			Double_t oben_inv_proton=135+20;
			Double_t unten_inv_proton=135-20;

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
			
			Double_t scalbackground_ppi0=-0.045*0.5;
			
			Double_t scalbackgroundpb = -0.045*pb;
			Double_t scalbackgroundpt = -0.045*pt;
			
			Double_t scalbackgroundpb_ppi0 = -0.5*0.0625*pb;
			Double_t scalbackgroundpt_ppi0 = -0.5*0.0625*pt;
			
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
					if(anglethetaproton_meas > 90){continue;}

					//...........................Histogramme vor jeglichen Cuts erstellen (3PED)..................................

					HistoManu::FillTH1_timeweighted(coplanarity_collerated,phi,time);
					HistoManu::FillTH1_timeweighted(thetaproton_collerated,thetadiff,time);
		
					HistoManu::FillTH1_timeweighted(missingmass_collerated_proton,missingmass,time);
		
					HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse34,time);
					HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse24,time);
					HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse23,time);
					HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse12,time);
					HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse13,time);
					HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse14,time);
		
					HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse12,invariantemasse34,time);
					HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse34,invariantemasse12,time);
					HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse24,invariantemasse13,time);
					HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse13,invariantemasse24,time);
					HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse14,invariantemasse23,time);
					HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse23,invariantemasse14,time);
	

					//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::single cuts:::::::::::::::::::::::::::::::::::::::::
					
					//__________________Polar angle______________________________________________________

					if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {
					
						HistoManu::FillTH1_timeweighted(coplanarity_theta_collerated,phi,time);
						HistoManu::FillTH1_timeweighted(missingmass_theta_collerated,missingmass,time);

						HistoManu::FillTH1_timeweighted(massesumme_theta_collerated,invariantemasse12,time);
						HistoManu::FillTH1_timeweighted(massesumme_theta_collerated,invariantemasse34,time);
						HistoManu::FillTH1_timeweighted(massesumme_theta_collerated,invariantemasse13,time);
						HistoManu::FillTH1_timeweighted(massesumme_theta_collerated,invariantemasse24,time);
						HistoManu::FillTH1_timeweighted(massesumme_theta_collerated,invariantemasse14,time);
						HistoManu::FillTH1_timeweighted(massesumme_theta_collerated,invariantemasse23,time);
					

						HistoManu::FillTH2_timeweighted(massegegenmasse_theta_collerated,invariantemasse12,invariantemasse34,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_theta_collerated,invariantemasse34,invariantemasse12,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_theta_collerated,invariantemasse13,invariantemasse24,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_theta_collerated,invariantemasse24,invariantemasse13,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_theta_collerated,invariantemasse23,invariantemasse14,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_theta_collerated,invariantemasse14,invariantemasse23,time);
	
					}

					//__________________Coplanarity______________________________________________________
					
					if((phi > ((unten_copl)) && phi < ((oben_copl)))){

						HistoManu::FillTH1_timeweighted(missingmass_copl_collerated,missingmass,time);

						HistoManu::FillTH1_timeweighted(massesumme_copl_collerated,invariantemasse12,time);
						HistoManu::FillTH1_timeweighted(massesumme_copl_collerated,invariantemasse34,time);
						HistoManu::FillTH1_timeweighted(massesumme_copl_collerated,invariantemasse13,time);
						HistoManu::FillTH1_timeweighted(massesumme_copl_collerated,invariantemasse24,time);
						HistoManu::FillTH1_timeweighted(massesumme_copl_collerated,invariantemasse14,time);
						HistoManu::FillTH1_timeweighted(massesumme_copl_collerated,invariantemasse23,time);

						HistoManu::FillTH2_timeweighted(massegegenmasse_copl_collerated,invariantemasse12,invariantemasse34,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_copl_collerated,invariantemasse34,invariantemasse12,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_copl_collerated,invariantemasse13,invariantemasse24,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_copl_collerated,invariantemasse24,invariantemasse13,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_copl_collerated,invariantemasse23,invariantemasse14,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_copl_collerated,invariantemasse14,invariantemasse23,time);
					
						
						if(GetRootinos()->HasTAPS(0)){
							HistoManu::FillTH1_timeweighted(thetaproton_copl_taps_collerated,thetadiff,time);
						}
						else
						{
							HistoManu::FillTH1_timeweighted(thetaproton_copl_collerated,thetadiff,time);
						}
					
					}
										
					//__________________invariant mass______________________________________________________
					
					if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){

						HistoManu::FillTH1_timeweighted(coplanarity_inv_collerated,phi,time);
						HistoManu::FillTH1_timeweighted(missingmass_inv_collerated_proton,missingmass,time);

					
						if(GetRootinos()->HasTAPS(0)){
							HistoManu::FillTH1_timeweighted(thetaproton_inv_taps_collerated,thetadiff,time);
						}
						else
						{
							HistoManu::FillTH1_timeweighted(thetaproton_inv_collerated,thetadiff,time);
						}

					}

		// 			::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::aufeinanderfolgende Cuts:::::::::::::::::::::::::::::::::::::::::
					

// 			--------------------------------------------------------CopLanaritY-Cuts----------------------
					
		// 			__________________Missing Mass______________________________________________________
					
					if(missingmass>(unten_mass) && missingmass<(oben_mass)){
					
						HistoManu::FillTH1_timeweighted(coplanarity_mass_collerated,phi,time);


		// 			__________________Polarwinkel______________________________________________________
					
						if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {

							HistoManu::FillTH1_timeweighted(coplanarity_mass_theta_collerated,phi,time);
					
						
		// 			__________________invariante Masse______________________________________________________
							if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
							//if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
								HistoManu::FillTH1_timeweighted(coplanarity_mass_theta_inv_collerated,phi,time);
								HistoManu::FillTH2_timeweighted(coplanarityverteilung_collerated,phi,beamphoton1E,time);
									
							}
						}
					}
										
// 			--------------------------------------------------------Missing Mass Cuts----------------------
															
		// 			__________________Polarwinkel______________________________________________________
					
					if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {
					
						HistoManu::FillTH1_timeweighted(missingmass_theta_collerated,missingmass,time);

		// 			__________________Coplanarity______________________________________________________
					
						if((phi > ( (unten_copl)) && phi < ((oben_copl)))){

							HistoManu::FillTH1_timeweighted(missingmass_theta_copl_collerated,missingmass,time);
					
					
		// 			__________________invariante Masse______________________________________________________
							if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
							//if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
								HistoManu::FillTH1_timeweighted(missingmass_theta_copl_inv_collerated,missingmass,time);
								HistoManu::FillTH2_timeweighted(missingmassverteilung_collerated_proton,missingmass,beamphoton1E,time);
							}
						}
					}
					
// 			--------------------------------------------------------PolarWinkel Cuts----------------------			
					
		// 			__________________Missing Mass______________________________________________________
					
					if(missingmass>(unten_mass) && missingmass<(oben_mass)){
					
						if(GetRootinos()->HasTAPS(0)){
							HistoManu::FillTH1_timeweighted(thetaproton_mass_taps_collerated,thetadiff,time);
						}
						else{
							HistoManu::FillTH1_timeweighted(thetaproton_mass_collerated,thetadiff,time);		
						}		
					
		// 			__________________Coplanarity______________________________________________________
					
						if((phi > ((unten_copl)) && phi < ((oben_copl)))){

							if(GetRootinos()->HasTAPS(0)){
								HistoManu::FillTH1_timeweighted(thetaproton_mass_copl_taps_collerated,thetadiff,time);
							}
							else{
								HistoManu::FillTH1_timeweighted(thetaproton_mass_copl_collerated,thetadiff,time);		
							}	
					
					
		// 			__________________invariante Masse______________________________________________________
							if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
// 							if(((invariantemasse1 > (pionmasse + unten_inv) && invariantemasse1 < (pionmasse + oben_inv)) && (invariantemasse2 > (pionmasse + unten_inv) && invariantemasse2 < (pionmasse + oben_inv)))){
					
								if(GetRootinos()->HasTAPS(0)){
									HistoManu::FillTH1_timeweighted(thetaproton_mass_copl_inv_taps_collerated,thetadiff,time);
									HistoManu::FillTH2_timeweighted(thetaverteilung_collerated_taps,thetadiff,beamphoton1E,time);
								}
								else{
									HistoManu::FillTH1_timeweighted(thetaproton_mass_copl_inv_collerated,thetadiff,time);	
									HistoManu::FillTH2_timeweighted(thetaverteilung_collerated,thetadiff,beamphoton1E,time);	
								}
							}
		
						}	
					} 
// 			................................invarianteMasse Cut...................................................................................
										
		// 			__________________Missing Mass______________________________________________________
					
					if(missingmass>(unten_mass) && missingmass<(oben_mass)){

						HistoManu::FillTH1_timeweighted(massesumme_mass_collerated_proton,invariantemasse12,time);
						HistoManu::FillTH1_timeweighted(massesumme_mass_collerated_proton,invariantemasse34,time);
						HistoManu::FillTH1_timeweighted(massesumme_mass_collerated_proton,invariantemasse13,time);
						HistoManu::FillTH1_timeweighted(massesumme_mass_collerated_proton,invariantemasse24,time);
						HistoManu::FillTH1_timeweighted(massesumme_mass_collerated_proton,invariantemasse14,time);
						HistoManu::FillTH1_timeweighted(massesumme_mass_collerated_proton,invariantemasse23,time);

						HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated_proton,invariantemasse12,invariantemasse34,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated_proton,invariantemasse34,invariantemasse12,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated_proton,invariantemasse13,invariantemasse24,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated_proton,invariantemasse24,invariantemasse13,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated_proton,invariantemasse23,invariantemasse14,time);
						HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated_proton,invariantemasse14,invariantemasse23,time);
									
					
		// 			__________________Polarwinkel______________________________________________________
					
						if((phi > ((unten_copl)) && phi < (oben_copl))){

							HistoManu::FillTH1_timeweighted(massesumme_mass_copl_collerated,invariantemasse12,time);
							HistoManu::FillTH1_timeweighted(massesumme_mass_copl_collerated,invariantemasse34,time);
							HistoManu::FillTH1_timeweighted(massesumme_mass_copl_collerated,invariantemasse13,time);
							HistoManu::FillTH1_timeweighted(massesumme_mass_copl_collerated,invariantemasse24,time);
							HistoManu::FillTH1_timeweighted(massesumme_mass_copl_collerated,invariantemasse14,time);
							HistoManu::FillTH1_timeweighted(massesumme_mass_copl_collerated,invariantemasse23,time);

							HistoManu::FillTH2_timeweighted(massegegenmasse_mass_copl_collerated,invariantemasse12,invariantemasse34,time);
							HistoManu::FillTH2_timeweighted(massegegenmasse_mass_copl_collerated,invariantemasse34,invariantemasse12,time);
							HistoManu::FillTH2_timeweighted(massegegenmasse_mass_copl_collerated,invariantemasse13,invariantemasse24,time);
							HistoManu::FillTH2_timeweighted(massegegenmasse_mass_copl_collerated,invariantemasse24,invariantemasse13,time);
							HistoManu::FillTH2_timeweighted(massegegenmasse_mass_copl_collerated,invariantemasse23,invariantemasse14,time);
							HistoManu::FillTH2_timeweighted(massegegenmasse_mass_copl_collerated,invariantemasse14,invariantemasse23,time);		
						
		// 			__________________Coplanarity______________________________________________________
																							
							if(thetadiff > (unten_theta) && thetadiff < (oben_theta)) {

								HistoManu::FillTH1_timeweighted(massesumme_mass_theta_copl_collerated,invariantemasse12,time);
								HistoManu::FillTH1_timeweighted(massesumme_mass_theta_copl_collerated,invariantemasse34,time);
								HistoManu::FillTH1_timeweighted(massesumme_mass_theta_copl_collerated,invariantemasse13,time);
								HistoManu::FillTH1_timeweighted(massesumme_mass_theta_copl_collerated,invariantemasse24,time);
								HistoManu::FillTH1_timeweighted(massesumme_mass_theta_copl_collerated,invariantemasse14,time);
								HistoManu::FillTH1_timeweighted(massesumme_mass_theta_copl_collerated,invariantemasse23,time);

								HistoManu::FillTH2_timeweighted(massegegenmasse_mass_theta_copl_collerated,invariantemasse12,invariantemasse34,time);
								HistoManu::FillTH2_timeweighted(massegegenmasse_mass_theta_copl_collerated,invariantemasse34,invariantemasse12,time);
								HistoManu::FillTH2_timeweighted(massegegenmasse_mass_theta_copl_collerated,invariantemasse13,invariantemasse24,time);
								HistoManu::FillTH2_timeweighted(massegegenmasse_mass_theta_copl_collerated,invariantemasse24,invariantemasse13,time);
								HistoManu::FillTH2_timeweighted(massegegenmasse_mass_theta_copl_collerated,invariantemasse23,invariantemasse14,time);
								HistoManu::FillTH2_timeweighted(massegegenmasse_mass_theta_copl_collerated,invariantemasse14,invariantemasse23,time);
	
								
							}							    				
						}		
					}
					//stuff after all cuts and creation of needed analysis histograms
	
					if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && (((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv)))){

						
						fitter.Set(photon1_4vektor, photon2_4vektor, photon3_4vektor, photon4_4vektor, proton_4vect_meas, beam_4vect);

						if(!fitter.Solve())	return;

						if(GetRootinos()->HasCB(0)){
							HistoManu::FillTH2_timeweighted(Check_CBdE_E,proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
						}

						if(GetRootinos()->HasTAPS(0)){
							HistoManu::FillTH2_timeweighted(Check_TAPSdE_E,proton_4vect_meas.E(),GetRootinos()->GetVetoEnergy(0),time);
						}

						if(cbtrigger=kTRUE || GetScalers()->GetNEntries()==0){
							HistoManu::FillTH1_timeweighted(h_energy_sum_2pi0_5ped,GetTrigger()->GetEnergySum(),time);
							HistoManu::FillTH1_timeandvalueweighted(h_energy_sum_2pi0_5ped_weightcos2pi0,GetTrigger()->GetEnergySum(),time,dwq(beamphoton1E, cospi0pi0, "2pi0",""));
							HistoManu::FillTH1_timeandvalueweighted(h_energy_sum_2pi0_5ped_weightcosppi0,GetTrigger()->GetEnergySum(),time,dwq(beamphoton1E,  0.5*(cosppi01+cosppi02), "2pi0","cosppi0"));
	
							HistoManu::FillTH3_timeweighted(cos2pi0_beamphoton_energysum_rek,cospi0pi0,beamphoton1E,GetTrigger()->GetEnergySum(),time);
							HistoManu::FillTH3_timeweighted(m2pi0_beamphoton_energysum_rek,invariantmass_pi0pi0,beamphoton1E,GetTrigger()->GetEnergySum(),time);
							HistoManu::FillTH3_timeandvalueweighted(cosppi0_beamphoton_energysum_rek,cosppi01,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);
							HistoManu::FillTH3_timeandvalueweighted(cosppi0_beamphoton_energysum_rek,cosppi02,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);
							HistoManu::FillTH3_timeandvalueweighted(mppi0_beamphoton_energysum_rek,invariantmass_ppi01,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);
							HistoManu::FillTH3_timeandvalueweighted(mppi0_beamphoton_energysum_rek,invariantmass_ppi02,beamphoton1E,GetTrigger()->GetEnergySum(),time,0.5);

							for(Int_t i=0; i<GetTracks()->GetNTracks() ;i++){//loop over all GetTracks()
								if(GetTracks()->HasTAPS(i)){//if track is from taps
									for(Int_t j=0; j<GetDetectorHits()->GetNBaF2Hits() ;j++){//loop over all hits
										if(GetDetectorHits()->GetBaF2Hits(j)!=GetTracks()->GetCentralCrystal(i)){continue;}//if central crystal index is unequal baf2 crystal index skip 	
									HistoManu::FillTH2_timeweighted(h_energy_taps_2pi0_5ped,GetDetectorHits()->GetBaF2Energy(j),GetTrigger()->GetEnergySum(),time);
									}
					
								}
							}
						}

						HistoManu::FillTH2_timeweighted(cosverteilung_collerated_proton,cospi0pi0,beamphoton1E,time);
						HistoManu::FillTH2_timeweighted(invverteilung_collerated_proton,invariantmass_pi0pi0,beamphoton1E,time);
						HistoManu::FillTH2_timeandvalueweighted(cosppi0verteilung_collerated_proton,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH2_timeandvalueweighted(cosppi0verteilung_collerated_proton,cosppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH2_timeandvalueweighted(invppi0verteilung_collerated_proton,invariantmass_ppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH2_timeandvalueweighted(invppi0verteilung_collerated_proton,invariantmass_ppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH2_timeweighted(cosverteilung_thetaproton_collerated_proton,cospi0pi0,anglethetaproton_meas, time);;					
						if(pt>0 || pt==5){
				
 							if(planesetting=="PARA"){

								HistoManu::FillTH3_timeweighted(kristallminus_targetplus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
								HistoManu::FillTH3_timeweighted(kristallminus_targetplus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

							}

							if(planesetting=="PERP"){

								HistoManu::FillTH3_timeweighted(kristallplus_targetplus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
								HistoManu::FillTH3_timeweighted(kristallplus_targetplus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);
							}
						}
		
						if(pt<0 || pt==5){

 							if(planesetting=="PARA"){

								HistoManu::FillTH3_timeweighted(kristallminus_targetminus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
								HistoManu::FillTH3_timeweighted(kristallminus_targetminus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

							}

 							if(planesetting=="PERP"){

								HistoManu::FillTH3_timeweighted(kristallplus_targetminus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
								HistoManu::FillTH3_timeweighted(kristallplus_targetminus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
								HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

							}
						}
					}

// 				
				}
//_________________________________________________ end of all events with proton identified!!!!!_________

	

//...........................Histogramme vor jeglichen Cuts erstellen(2+3ped)..................................
			HistoManu::FillTH1_timeweighted(missingmass_collerated_proton,missingmass,time);

			HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse34,time);
			HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse24,time);
			HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse23,time);
			HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse12,time);
			HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse13,time);
			HistoManu::FillTH1_timeweighted(massesumme_collerated_proton,invariantemasse14,time);

			HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse12,invariantemasse34,time);
			HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse34,invariantemasse12,time);
			HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse24,invariantemasse13,time);
			HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse13,invariantemasse24,time);
			HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse14,invariantemasse23,time);
			HistoManu::FillTH2_timeweighted(massegegenmasse_collerated_proton,invariantemasse23,invariantemasse14,time);

			if(GetRootinos()->GetNParticles()==1){

					if(!(thetadiff > (unten_theta) && thetadiff < (oben_theta)) ||  !(phi > (unten_copl) && phi < (oben_copl))){continue;} 

			}					
			
// 			::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::einzelne Cuts:::::::::::::::::::::::::::::::::::::::::
			
			
// 			__________________invariante Masse______________________________________________________
			if(((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv))){
			//if(((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)) && (invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){

				HistoManu::FillTH1_timeweighted(missingmass_inv_collerated,missingmass,time);
				HistoManu::FillTH2_timeweighted(missingmassverteilung_collerated,missingmass,beamphoton1E,time);

			}

			if(missingmass>(unten_mass) && missingmass<(oben_mass)){
									
				HistoManu::FillTH1_timeweighted(massesumme_mass_collerated,invariantemasse34,time);
				HistoManu::FillTH1_timeweighted(massesumme_mass_collerated,invariantemasse24,time);
				HistoManu::FillTH1_timeweighted(massesumme_mass_collerated,invariantemasse23,time);
				HistoManu::FillTH1_timeweighted(massesumme_mass_collerated,invariantemasse12,time);
				HistoManu::FillTH1_timeweighted(massesumme_mass_collerated,invariantemasse13,time);
				HistoManu::FillTH1_timeweighted(massesumme_mass_collerated,invariantemasse14,time);
	
				HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated,invariantemasse12,invariantemasse34,time);
				HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated,invariantemasse34,invariantemasse12,time);
				HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated,invariantemasse24,invariantemasse13,time);
				HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated,invariantemasse13,invariantemasse24,time);
				HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated,invariantemasse14,invariantemasse23,time);
				HistoManu::FillTH2_timeweighted(massegegenmasse_mass_collerated,invariantemasse23,invariantemasse14,time);
			}
					
						

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& (((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv)))&&GetTrigger()->GetNErrors()!=0){
				events_witherror->Fill(phimeson,cospi0pi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((invariantemasse1 > (unten_inv) && invariantemasse1 < (oben_inv)))&& ((invariantemasse2 > (unten_inv) && invariantemasse2 < (oben_inv)))){
				events_all->Fill(phimeson,cospi0pi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& (((invariantemasse12 > unten_inv && invariantemasse12 < oben_inv) && (invariantemasse34 > unten_inv && invariantemasse34 < oben_inv)) || ((invariantemasse13 > unten_inv && invariantemasse13 < oben_inv) && (invariantemasse24 > unten_inv && invariantemasse24 < oben_inv)) || ((invariantemasse14 > unten_inv && invariantemasse14 < 170) && (invariantemasse23 > unten_inv && invariantemasse23 < oben_inv)))){

				HistoManu::FillTH2_timeweighted(cosverteilung_collerated,cospi0pi0,beamphoton1E,time);
				HistoManu::FillTH2_timeweighted(invverteilung_collerated,invariantmass_pi0pi0,beamphoton1E,time);
				HistoManu::FillTH2_timeandvalueweighted(cosppi0verteilung_collerated,cosppi01,beamphoton1E,time,0.5);
				HistoManu::FillTH2_timeandvalueweighted(cosppi0verteilung_collerated,cosppi02,beamphoton1E,time,0.5);
				HistoManu::FillTH2_timeandvalueweighted(invppi0verteilung_collerated,invariantmass_ppi01,beamphoton1E,time,0.5);
				HistoManu::FillTH2_timeandvalueweighted(invppi0verteilung_collerated,invariantmass_ppi02,beamphoton1E,time,0.5);
			
 				if(pt>0 || pt==5){
				
 					if(planesetting=="PARA"){

						HistoManu::FillTH3_timeweighted(kristallminus_targetplus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						HistoManu::FillTH3_timeweighted(kristallminus_targetplus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetplus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

					}

					if(planesetting=="PERP"){

						HistoManu::FillTH3_timeweighted(kristallplus_targetplus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						HistoManu::FillTH3_timeweighted(kristallplus_targetplus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetplus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);
					}
				}
		
				if(pt<0 || pt==5){

					if(planesetting=="PARA"){

						HistoManu::FillTH3_timeweighted(kristallminus_targetminus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						HistoManu::FillTH3_timeweighted(kristallminus_targetminus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallminus_targetminus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

					}

 					if(planesetting=="PERP"){

						HistoManu::FillTH3_timeweighted(kristallplus_targetminus_collerated_proton,phimeson,cospi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pb_proton,phimeson,cospi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_collerated_pt_proton,phimeson,cospi0pi0,beamphoton1E,time,pt);
	
						HistoManu::FillTH3_timeweighted(kristallplus_targetminus_im_collerated_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_collerated_pb_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pb);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_collerated_pt_proton,phimeson,invariantmass_pi0pi0,beamphoton1E,time,pt);

						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_proton,phippi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pb_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pt_proton,phippi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_proton,phippi02,cosppi02,beamphoton1E,time,0.5);
							HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pb_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_ppi0_collerated_pt_proton,phippi02,cosppi02,beamphoton1E,time,scalsignalpt_ppi0);

						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pb_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pt_proton,invariantmass_ppi01,cosppi01,beamphoton1E,time,scalsignalpt_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,0.5);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pb_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpb_ppi0);
						HistoManu::FillTH3_timeandvalueweighted(kristallplus_targetminus_im_ppi0_collerated_pt_proton,phippi02,invariantmass_ppi02,beamphoton1E,time,scalsignalpt_ppi0);

					}
				}					
			}//end of analysis histos

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
cos2pi0_beamphoton_energysum_monte->Write();
cosppi0_beamphoton_energysum_monte->Write();

m2pi0_beamphoton_energysum_monte->Write();
mppi0_beamphoton_energysum_monte->Write();

//reconstructed
cos2pi0_beamphoton_energysum_rek->Write();
cosppi0_beamphoton_energysum_rek->Write();

m2pi0_beamphoton_energysum_rek->Write();
mppi0_beamphoton_energysum_rek->Write();

curDir1->cd();
poltable_energy->Write();
poltable_energy_weight->Write();
filenumber.Write();
TDirectory* curDir3  = curDir1->mkdir("Selektion");
curDir3->cd();
time_prompt->Write();
time1->Write();
time_side->Write();
//
missingmass_collerated->Write();
missingmass_inv_collerated->Write();
//
massesumme_collerated->Write();
massesumme_mass_collerated->Write();
massegegenmasse_collerated->Write();
massegegenmasse_mass_collerated->Write();


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
