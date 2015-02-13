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
Check_TAPS_TOF_photons_3ped		= new TH2F("Check_TAPS_TOF_photon_3ped", "TOF analysis (photon_3ped);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);
Check_TAPS_TOF_photons_2_3ped		= new TH2F("Check_TAPS_TOF_photon_2_3ped", "TOF analysis (photon_2_3ped);t_{tof} [ns];E_{dep} [Mev]", 	240, 0, 12, 600, 0, 600);

test = new TH1F("test","test",1100,-5.5,5.5);

//kinematic variables in dependence of energy
IM_energy_kplustplus = new TH3F("IM_energy_kplustplus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kplustminus = new TH3F("IM_energy_kplustminus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustplus = new TH3F("IM_energy_kminustplus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustminus = new TH3F("IM_energy_kminustminus", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);

MM_energy_kplustplus = new TH3F("MM_energy_kplustplus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kplustminus = new TH3F("MM_energy_kplustminus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustplus = new TH3F("MM_energy_kminustplus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustminus = new TH3F("MM_energy_kminustminus", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);

// //kinematic variables in dependence of energy
IM_energy_kplustplus_proton = new TH3F("IM_energy_kplustplus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kplustminus_proton = new TH3F("IM_energy_kplustminus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustplus_proton = new TH3F("IM_energy_kminustplus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);
IM_energy_kminustminus_proton = new TH3F("IM_energy_kminustminus_proton", "IM_energy; m_{#gamma,#gamma} [MeV]; cos(#theta_{#pi});E_{beam} [MeV]", 1000,0,1000,72,-1,1,200,200,800);

MM_energy_kplustplus_proton = new TH3F("MM_energy_kplustplus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kplustminus_proton = new TH3F("MM_energy_kplustminus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustplus_proton = new TH3F("MM_energy_kminustplus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);
MM_energy_kminustminus_proton = new TH3F("MM_energy_kminustminus_proton", "MM_energy; m_{mm} [MeV];cos(#theta_{#pi});E_{beam} [MeV]", 400,800,1200,72,-1,1,200,200,800);

coplanarityverteilung_collerated = new TH2F("coplanarityverteilung_collerated", "Coplanarity-Verteilung;#phi_{#pi}-#phi_{p}[deg]", 400,0,360,200,200,800);
coplanarityverteilung_collerated->Sumw2();

thetaverteilung_collerated = new TH2F("thetaverteilung_collerated", "Polar-Verteilung;#Delt #theta_{p}[deg]", 400,-200, 200,200,200,800);
thetaverteilung_collerated->Sumw2();

missingmassverteilung_collerated = new TH2F("missingmassverteilung_collerated", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,200,200,800);
missingmassverteilung_collerated->Sumw2();

missingmassverteilung_collerated_proton = new TH2F("missingmassverteilung_collerated_proton", "Missing Mass-Verteilung;m_{mm} [MeV]", 500,0,2000,200,200,800);
missingmassverteilung_collerated_proton->Sumw2();

invmassverteilung_collerated = new TH2F("invmassverteilung_collerated", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,200,200,800);
invmassverteilung_collerated->Sumw2();

invmassverteilung_collerated_proton = new TH2F("invmassverteilung_collerated_proton", "InvMass-Verteilung;m_{#gamma #gamma} [MeV]", 1000,0,1000,200,200,800);
invmassverteilung_collerated_proton->Sumw2();

//kinematic variables in dependence of energy
cosverteilung_collerated = new TH2F("cosverteilung_collerated", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosverteilung_collerated->Sumw2();

// //kinematic variables in dependence of energy
cosverteilung_collerated_proton = new TH2F("cosverteilung_collerated_proton", "Cos-Verteilung; cos(#theta_{#pi});E_{beam} [MeV]", 72,-1,1,200,200,800);
cosverteilung_collerated_proton->Sumw2();
// 
// 
// 

events_all = new TH3F("events_all","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
events_witherror = new TH3F("events_witherror","events with error/events all; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);

kristallminus_targetplus_collerated_proton = new TH3F("kristallminus_targetplus_collerated_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_proton->Sumw2();

kristallminus_targetminus_collerated_proton = new TH3F("kristallminus_targetminus_collerated_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_proton->Sumw2();

kristallplus_targetplus_collerated_proton = new TH3F("kristallplus_targetplus_collerated_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_proton->Sumw2();

kristallplus_targetminus_collerated_proton = new TH3F("kristallplus_targetminus_collerated_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_proton->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb_proton = new TH3F("kristallminus_targetplus_collerated_pb_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pb_proton->Sumw2();

kristallminus_targetminus_collerated_pb_proton = new TH3F("kristallminus_targetminus_collerated_pb_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pb_proton->Sumw2();

kristallplus_targetplus_collerated_pb_proton = new TH3F("kristallplus_targetplus_collerated_pb_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pb_proton->Sumw2();

kristallplus_targetminus_collerated_pb_proton = new TH3F("kristallplus_targetminus_collerated_pb_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pb_proton->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt_proton = new TH3F("kristallminus_targetplus_collerated_pt_proton","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pt_proton->Sumw2();

kristallminus_targetminus_collerated_pt_proton = new TH3F("kristallminus_targetminus_collerated_pt_proton","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pt_proton->Sumw2();

kristallplus_targetplus_collerated_pt_proton = new TH3F("kristallplus_targetplus_collerated_pt_proton","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pt_proton->Sumw2();

kristallplus_targetminus_collerated_pt_proton = new TH3F("kristallplus_targetminus_collerated_pt_proton","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pt_proton->Sumw2();

//Histogramme für 2pi0 System mit Cos-Verteilung___________________________________________________

kristallminus_targetplus_collerated = new TH3F("kristallminus_targetplus_collerated","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated->Sumw2();

kristallminus_targetminus_collerated = new TH3F("kristallminus_targetminus_collerated","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated->Sumw2();

kristallplus_targetplus_collerated = new TH3F("kristallplus_targetplus_collerated","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated->Sumw2();

kristallplus_targetminus_collerated = new TH3F("kristallplus_targetminus_collerated","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated->Sumw2();

//Histogramme gewichtet mit Beam Polarisation
kristallminus_targetplus_collerated_pb = new TH3F("kristallminus_targetplus_collerated_pb","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pb->Sumw2();

kristallminus_targetminus_collerated_pb = new TH3F("kristallminus_targetminus_collerated_pb","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pb->Sumw2();

kristallplus_targetplus_collerated_pb = new TH3F("kristallplus_targetplus_collerated_pb","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pb->Sumw2();

kristallplus_targetminus_collerated_pb = new TH3F("kristallplus_targetminus_collerated_pb","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pb->Sumw2();

//Histogramme gewichtet mit Target Polarisation
kristallminus_targetplus_collerated_pt = new TH3F("kristallminus_targetplus_collerated_pt","Kristall-45 Target positiv; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS};E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetplus_collerated_pt->Sumw2();

kristallminus_targetminus_collerated_pt = new TH3F("kristallminus_targetminus_collerated_pt","Kristall-45 Target negativ; #phi_{#pi}[deg];cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallminus_targetminus_collerated_pt->Sumw2();

kristallplus_targetplus_collerated_pt = new TH3F("kristallplus_targetplus_collerated_pt","Kristall+45 Target positiv; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]",96,-180,180,72,-1,1,200,200,800);
kristallplus_targetplus_collerated_pt->Sumw2();

kristallplus_targetminus_collerated_pt = new TH3F("kristallplus_targetminus_collerated_pt","Kristall+45 Target negativ; #phi_{#pi}[deg]; cos(#theta)_{#pi}^{CMS}; E_{beam} [MeV]", 96,-180,180,72,-1,1,200,200,800);
kristallplus_targetminus_collerated_pt->Sumw2();


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
pt=PPi0Analyse::targetpol(file_in);

// planesetting=PPi0Analyse::polplane(file_in);
if(linpol->GetPolPlane()==1){planesetting="PERP";}
if(linpol->GetPolPlane()==0){planesetting="PARA";}
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();
//FOR CARBON
// polsetting=PPi0Analyse::poledge(file_in);
// TString *data=new TString(polsetting);
// TDirectory* curDir1  = file_out->mkdir(Form("%s",data->Data()));

TDirectory* curDir1  = file_out->mkdir(Form("Butanol%i",(Int_t)linpol->GetEdgeSetting()));

TString* filename1 = new TString(file_in->GetPath());
TString* path11 = new TString(filename1->Tokenize("_")->At(filename1->Tokenize("_")->GetEntries()-1)->GetName());
path11->Resize(path11->Length()-7);
TNamed filenumber=TNamed(path11->Data(), "Filenumber");

curDir1->cd();
poltable_energy->Write();
poltable_energy_weight->Write();
filenumber.Write();
triggertest->Write();
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
Check_TAPS_TOF_photons_2_3ped->Write();

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
Check_TAPS_TOF_photons_3ped->Write();
Check_TAPS_TOF_proton->Write();

MM_energy_kplustplus_proton->Write();
MM_energy_kplustminus_proton->Write();
MM_energy_kminustplus_proton->Write();
MM_energy_kminustminus_proton->Write();


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
cosverteilung_collerated_proton->Write();
thetaverteilung_collerated->Write();
coplanarityverteilung_collerated->Write();
missingmassverteilung_collerated_proton->Write();
invmassverteilung_collerated_proton->Write();


file_out->Close();

	return kTRUE;
}

void PPi0Analyse::fOnBeforeEventProcessing() {

}


void	PPi0Analyse::ProcessEvent()	
{

// 			cout << pt << "\t" << planesetting << "\t" << polsetting << endl; 
for(Int_t i=0;i<trigger->GetNTriggerPattern();i++){	
triggertest->Fill(trigger->GetTriggerPattern(i));
}
	TLorentzVector pi0_4vect = photons->Particle(0)+photons->Particle(1);
	//invMass
	Double_t inv=pi0_4vect.M();

	//TEST the target polarization (no entry with pt=0)
	test->Fill(pt);

		for(Int_t j=0; j < tagger->GetNTagged();j++)
		{

 		if(pt==5){continue;}
//  		if(pt!=5){continue;}//only carbon runs

			Double_t time= tagger->GetTagged_t(j) - 0.5*(photons->GetTime(1)+photons->GetTime(0));
			time1->Fill(time);

			if(time > -20 && time < 10){
				time_prompt->Fill(time);
			}

				if((time > -200 && time < -100) ||(time > 100 && time < 300)){	
					time_side->Fill(time);
				}
			if(linpol->GetPolDegree(tagger->GetTagged_ch(j))==-1){continue;}//skip events with wrong polarisation!

			poltable_energy->Fill(tagger->GetPhotonBeam_E(j));
			poltable_energy_weight->Fill(tagger->GetPhotonBeam_E(j),linpol->GetPolDegree(tagger->GetTagged_ch(j)));

			//get target, beam and missng particle information
			TLorentzVector protonvektor_target(0.,0.,0.,938.272046);
			TLorentzVector  beam_4vect = TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(j),tagger->GetPhotonBeam_E(j));
			TLorentzVector  missingp_4vect = beam_4vect + protonvektor_target - pi0_4vect;
			Double_t anglethetaproton_rek = TMath::RadToDeg()*missingp_4vect.Vect().Theta();
			Double_t beamphoton1E=tagger->GetPhotonBeam_E(j);
			Double_t pb=linpol->GetPolDegree(tagger->GetTagged_ch(j));

			//due to kinematic not possible -> kick them out
			if(anglethetaproton_rek > 90){continue;}
			
			if(beamphoton1E<200 || beamphoton1E>800){continue;}

			//energy dependent cuts
			Double_t oben_copl=270.07535+(-0.37891)*beamphoton1E+(0.000637701)*beamphoton1E*beamphoton1E+(-3.62165e-07)*beamphoton1E*beamphoton1E*beamphoton1E;
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

			//scaling for timebackground
			Double_t scalbackgroundpb = -0.075*pb;
			Double_t scalbackgroundpt = -0.075*pt;

			//Missing Mass
			Double_t missingmass=missingp_4vect.M();

			//Boost-System
			TLorentzVector pi0_4vect_boost = CMVector(pi0_4vect, beam_4vect, protonvektor_target);
			TVector3 pi0_3vektor_boost = pi0_4vect_boost.Vect();
			Double_t cospi0 = pi0_3vektor_boost.CosTheta();
			Double_t phimeson = TMath::RadToDeg()*(pi0_4vect.Vect().Phi());

			//without any cuts!
			if(time > -20 && time < 10){
				MM->Fill(missingmass);
				IM->Fill(inv);
			}
	
			if((time > -200 && time < -100) ||(time > 100 && time < 300)){	
				MM->Fill(missingmass,-0.075);
				IM->Fill(inv,-0.075);
			}


			//with cuts
			if(time > -20 && time < 10){
		
				if(((inv > (unten_inv) && inv < (oben_inv)))){

					if(planesetting=="PARA"){

						if(pt>0){
							MM_energy_kminustplus->Fill(missingmass,cospi0,beamphoton1E);
						}
	
						if(pt<0){
							MM_energy_kminustminus->Fill(missingmass,cospi0,beamphoton1E);
						}
					}

					if(planesetting=="PERP"){

						if(pt>0 ){
							MM_energy_kplustplus->Fill(missingmass,cospi0,beamphoton1E);
						}
				
						if(pt<0 ){
							MM_energy_kplustminus->Fill(missingmass,cospi0,beamphoton1E);
						}
					}
					MM_all->Fill(missingmass);
					missingmassverteilung_collerated->Fill(missingmass,beamphoton1E);
				}
	
				if(((missingmass > (unten_mass) && missingmass < (oben_mass)))){

					if(planesetting=="PARA"){

						if(pt>0 ){
							IM_energy_kminustplus->Fill(inv,cospi0,beamphoton1E);
						}
						if(pt<0 ){
							IM_energy_kminustminus->Fill(inv,cospi0,beamphoton1E);
						}
					}

					if(planesetting=="PERP"){

						if(pt>0 ){
							IM_energy_kplustplus->Fill(inv,cospi0,beamphoton1E);
						}
						if(pt<0 ){
							IM_energy_kplustminus->Fill(inv,cospi0,beamphoton1E);
						}
					}
					IM_all->Fill(inv);
					invmassverteilung_collerated->Fill(inv,beamphoton1E);
				}
			}//prompt end

			if((time > -200 && time < -100) ||(time > 100 && time < 300)){	
	
				if(((inv > (unten_inv) && inv < (oben_inv)))){

					if(planesetting=="PARA"){

						if(pt>0 ){
							MM_energy_kminustplus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
						}
						if(pt<0 ){
							MM_energy_kminustminus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
						}
					}

					if(planesetting=="PERP"){

						if(pt>0 ){
							MM_energy_kplustplus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
						}		
						if(pt<0 ){
							MM_energy_kplustminus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
						}
					}						
					MM_all->Fill(missingmass,-0.075);
					missingmassverteilung_collerated->Fill(missingmass,beamphoton1E,-0.075);
				}
	
				if(((missingmass > (unten_mass) && missingmass < (oben_mass)))){

					if(planesetting=="PARA"){

						if(pt>0 ){
							IM_energy_kminustplus->Fill(inv,cospi0,beamphoton1E,-0.075);
						}
						if(pt<0 ){
							IM_energy_kminustminus->Fill(inv,cospi0,beamphoton1E,-0.075);
						}
					}

					if(planesetting=="PERP"){

						if(pt>0 ){
							IM_energy_kplustplus->Fill(inv,cospi0,beamphoton1E,-0.075);
						}
						if(pt<0 ){
							IM_energy_kplustminus->Fill(inv,cospi0,beamphoton1E,-0.075);
						}
					}
					IM_all->Fill(inv,-0.075);
					invmassverteilung_collerated->Fill(inv,beamphoton1E,-0.075);
				}
			}//sideband end

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))&&trigger->GetNError()!=0){
				events_witherror->Fill(phimeson,cospi0, beamphoton1E);
			}

			if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){
				events_all->Fill(phimeson,cospi0, beamphoton1E);
			}

			if(time > -20 && time < 10){

				if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){

					if(photons->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){
	
						Check_TAPS_TOF_photons_2_3ped->Fill(((photons->GetTime(0)-tagger->GetTagged_t(j))/photons->Particle(0).Vect().Mag())+1/0.299792458,photons->Particle(0).E());
					}
	
					if(photons->GetApparatus(1) == GTreeRawEvent::APPARATUS_TAPS){
	
						Check_TAPS_TOF_photons_2_3ped->Fill(((photons->GetTime(1)-tagger->GetTagged_t(j))/photons->Particle(1).Vect().Mag())+1/0.299792458,photons->Particle(1).E());
					}

					cosverteilung_collerated->Fill(cospi0,beamphoton1E);

					if(planesetting=="PARA"){
		
						if(pt>0 ){
							kristallminus_targetplus_collerated->Fill(phimeson,cospi0, beamphoton1E);
							kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0, beamphoton1E,pb);
							kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0, beamphoton1E,pt);
						}

						if(pt<0 ){
							kristallminus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E);
							kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,pb);
							kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,pt);
						}
					}

					if(planesetting=="PERP"){
		
						if(pt>0 ){
							kristallplus_targetplus_collerated->Fill(phimeson,cospi0,beamphoton1E);
							kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,pb);
							kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,pt);;
						}
			
						if(pt<0 ){
							kristallplus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E);
 							kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,pb);
 							kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,pt);
						}
					}						
				}
			}//end prompt

			if((time > -200 && time < -100) ||(time > 100 && time < 300)){	

				if((missingmass>(unten_mass) && missingmass<(oben_mass))&& ((inv > (unten_inv) && inv < (oben_inv)))){

					if(photons->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){
	
						Check_TAPS_TOF_photons_2_3ped->Fill(((photons->GetTime(0)-tagger->GetTagged_t(j))/photons->Particle(0).Vect().Mag())+1/0.299792458,photons->Particle(0).E(),-0.075);
					}
	
					if(photons->GetApparatus(1) == GTreeRawEvent::APPARATUS_TAPS){
	
						Check_TAPS_TOF_photons_2_3ped->Fill(((photons->GetTime(1)-tagger->GetTagged_t(j))/photons->Particle(1).Vect().Mag())+1/0.299792458,photons->Particle(1).E(),-0.075);
					}

					cosverteilung_collerated->Fill(cospi0,beamphoton1E,-0.075);

					if(planesetting=="PARA"){
		
						if(pt>0  ){
							kristallminus_targetplus_collerated->Fill(phimeson,cospi0, beamphoton1E,-0.075);
							kristallminus_targetplus_collerated_pb->Fill(phimeson,cospi0, beamphoton1E,scalbackgroundpb);
							kristallminus_targetplus_collerated_pt->Fill(phimeson,cospi0, beamphoton1E,scalbackgroundpt);
						}
			
						if(pt<0  ){
							kristallminus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E,-0.075);
							kristallminus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpb);
							kristallminus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpt);
						}
					}

					if(planesetting=="PERP"){
		
						if(pt>0  ){
							kristallplus_targetplus_collerated->Fill(phimeson,cospi0,beamphoton1E,-0.075);
							kristallplus_targetplus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpb);
							kristallplus_targetplus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpt);;
						}

						if(pt<0  ){
							kristallplus_targetminus_collerated->Fill(phimeson,cospi0,beamphoton1E,-0.075);
							kristallplus_targetminus_collerated_pb->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpb);
							kristallplus_targetminus_collerated_pt->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpt);
						}
					}						
				}
			}//end sideband


//_________________________________only one charged particle__________________________________________________________________________

				if(protons->GetNParticles()==1 || electrons->GetNParticles()==1 || chargedPi->GetNParticles()==1 || rootinos->GetNParticles()==1){

					//get charged particle information
					TLorentzVector proton_4vect_meas;
					if(electrons->GetNParticles()==1){proton_4vect_meas=electrons->Particle(0);}
					if(protons->GetNParticles()==1){proton_4vect_meas=protons->Particle(0);}
					if(chargedPi->GetNParticles()==1){proton_4vect_meas=chargedPi->Particle(0);}
					if(rootinos->GetNParticles()==1){proton_4vect_meas=rootinos->Particle(0);}

					//apply TCUTGs
					if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_CB && cutcb->IsInside(proton_4vect_meas.E(),rootinos->Get_dE(0))==1){kCB=10;}else{kCB=10;}
					if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS && cuttaps->IsInside(proton_4vect_meas.E(),rootinos->Get_dE(0))==1){kTAPS=10;}else{kTAPS=100;}

					//Phi-Difference
					Double_t anglephimeson = pi0_4vect.Vect().Phi();
					Double_t anglephiproton = proton_4vect_meas.Vect().Phi();
					Double_t phi = 360*(anglephimeson - anglephiproton)/(2*TMath::Pi());
					if(phi < 0){phi = phi + 360;}
			
					//theta differenz
					Double_t anglethetaproton_meas = TMath::RadToDeg()*proton_4vect_meas.Vect().Theta();
					Double_t thetadiff=anglethetaproton_rek-anglethetaproton_meas;

					//due to kinematic not possible -> kick them out
					if(anglethetaproton_meas > 90){continue;}

					//without any cuts!
					if(time > -20 && time < 10){

						if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_CB){
							Check_CBdE_E_nocuts->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0));
						}

						if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){
							Check_TAPSdE_E_nocuts->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0));
						}
						MM_proton->Fill(missingmass);
						theta_proton->Fill(thetadiff);
						coplanarity_proton->Fill(phi);
						IM_proton->Fill(inv);
					}
		
					if((time > -200 && time < -100) ||(time > 100 && time < 300)){	

						if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_CB){
							Check_CBdE_E_nocuts->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0),-0.075);
						}

						if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){
							Check_TAPSdE_E_nocuts->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0),-0.075);
						}

						MM_proton->Fill(missingmass,-0.075);
						IM_proton->Fill(inv,-0.075);
						theta_proton->Fill(thetadiff,-0.075);
						coplanarity_proton->Fill(phi,-0.075);
					}
				
					if(time > -20 && time < 10){

						if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){
							theta_all_proton->Fill(thetadiff);
							thetaverteilung_collerated->Fill(thetadiff,beamphoton1E);
						}
	
						if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (thetadiff > (unten_theta) && thetadiff < (oben_theta)) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){
							coplanarity_all_proton->Fill(phi);
							coplanarityverteilung_collerated->Fill(phi, beamphoton1E);
						}

						if((thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

							if(planesetting=="PARA"){
								if(pt>0 ){
									MM_energy_kminustplus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
		
								if(pt<0 ){
									MM_energy_kminustminus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
							}	

							if(planesetting=="PERP"){
								if(pt>0 ){
									MM_energy_kplustplus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
					
								if(pt<0 ){
									MM_energy_kplustminus_proton->Fill(missingmass,cospi0,beamphoton1E);
								}
							}	
							MM_all_proton->Fill(missingmass);
							missingmassverteilung_collerated_proton->Fill(missingmass,beamphoton1E);
						}
	
						if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl)))&&(kTAPS==10 || kCB==10)){

							if(planesetting=="PARA"){
								if(pt>0 ){
									IM_energy_kminustplus_proton->Fill(inv,cospi0,beamphoton1E);
								}
								if(pt<0 ){
									IM_energy_kminustminus_proton->Fill(inv,cospi0,beamphoton1E);
								}
							}

							if(planesetting=="PERP"){
								if(pt>0 ){
									IM_energy_kplustplus_proton->Fill(inv,cospi0,beamphoton1E);
								}
								if(pt<0 ){
									IM_energy_kplustminus_proton->Fill(inv,cospi0,beamphoton1E);
								}
							}							
							IM_all_proton->Fill(inv);
							invmassverteilung_collerated_proton->Fill(inv,beamphoton1E);
						}

					}//prompt end

					if((time > -200 && time < -100) ||(time > 100 && time < 300)){

						if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){
							theta_all_proton->Fill(thetadiff,-0.075);
							thetaverteilung_collerated->Fill(thetadiff,beamphoton1E,-0.075);

						}
	
						if((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton)) && (thetadiff > (unten_theta) && thetadiff < (oben_theta)) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){
							coplanarity_all_proton->Fill(phi,-0.075);
							coplanarityverteilung_collerated->Fill(phi, beamphoton1E,-0.075);
						}

						if((thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

							if(planesetting=="PARA"){
								if(pt>0 ){
									MM_energy_kminustplus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
								}
								if(pt<0 ){
									MM_energy_kminustminus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
								}
							}

							if(planesetting=="PERP"){				
								if(pt>0 ){
									MM_energy_kplustplus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
								}		
								if(pt<0 ){
									MM_energy_kplustminus->Fill(missingmass,cospi0,beamphoton1E,-0.075);
								}
							}								
							MM_all_proton->Fill(missingmass,-0.075);
							missingmassverteilung_collerated_proton->Fill(missingmass,beamphoton1E,-0.075);
						}
	
						if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl)))&&(kTAPS==10 || kCB==10)){

							if(planesetting=="PARA"){
// 							if(linpol->GetPolPlane()==0){
								if(pt>0 ){
									IM_energy_kminustplus->Fill(inv,cospi0,beamphoton1E,-0.075);
								}
								if(pt<0 ){
									IM_energy_kminustminus->Fill(inv,cospi0,beamphoton1E,-0.075);
								}
							}

							if(planesetting=="PERP"){
// 							if(linpol->GetPolPlane()==1){
								if(pt>0 ){
									IM_energy_kplustplus->Fill(inv,cospi0,beamphoton1E,-0.075);
								}
								if(pt<0 ){
									IM_energy_kplustminus->Fill(inv,cospi0,beamphoton1E,-0.075);
								}
							}							
							IM_all_proton->Fill(inv,-0.075);
							invmassverteilung_collerated_proton->Fill(inv,beamphoton1E,-0.075);
						}
					}//sideband end

					//stuff after all cuts and creation of needed analysis histograms_______________________________
					if(time > -20 && time < 10){
	
						if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

						if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_CB){
							Check_CBdE_E->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0));
						}

						if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){

							Check_TAPSdE_E->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0));
							Check_TAPS_TOF_proton->Fill(((rootinos->GetTime(0)-tagger->GetTagged_t(j))/proton_4vect_meas.Vect().Mag())+1/0.299792458,proton_4vect_meas.E());
						}

						if(photons->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){

							Check_TAPS_TOF_photons_3ped->Fill(((photons->GetTime(0)-tagger->GetTagged_t(j))/photons->Particle(0).Vect().Mag())+1/0.299792458,photons->Particle(0).E());
						}

						if(photons->GetApparatus(1) == GTreeRawEvent::APPARATUS_TAPS){

							Check_TAPS_TOF_photons_3ped->Fill(((photons->GetTime(1)-tagger->GetTagged_t(j))/photons->Particle(1).Vect().Mag())+1/0.299792458,photons->Particle(1).E());
						}
		
						cosverteilung_collerated_proton->Fill(cospi0,beamphoton1E);
							
							if(planesetting=="PARA"){
// 							if(linpol->GetPolPlane()==0){
			
								if(pt>0 ){
									kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0, beamphoton1E);
									kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0, beamphoton1E,pb);
									kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0, beamphoton1E,pt);
								}
	
				
								if(pt<0 ){
									kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E);
									kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,pb);
									kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,pt);
								}
	
							}

							if(planesetting=="PERP"){
// 							if(linpol->GetPolPlane()==1){
			
								if(pt>0 ){
									kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E);
									kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,pb);
									kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,pt);;
								}
	
				
								if(pt<0 ){
									kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E);
									kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,pb);
									kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,pt);
								}
	
							}						
	
						}
	
					}//end prompt
	
					if((time > -200 && time < -100) ||(time > 100 && time < 300)){	
	
						if(((missingmass>(unten_mass_proton) && missingmass<(oben_mass_proton))  && thetadiff > (unten_theta) && thetadiff < (oben_theta)) && (phi > ((unten_copl)) && phi < ((oben_copl))) && ((inv > (unten_inv_proton) && inv < (oben_inv_proton)))&&(kTAPS==10 || kCB==10)){

							if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_CB){
								Check_CBdE_E->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0),-0.075);
							}
	
							if(rootinos->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){
								Check_TAPSdE_E->Fill(proton_4vect_meas.E(),rootinos->Get_dE(0),-0.075);

								Check_TAPS_TOF_proton->Fill(((rootinos->GetTime(0)-tagger->GetTagged_t(j))/proton_4vect_meas.Vect().Mag())+1/0.299792458,proton_4vect_meas.E(),-0.075);

							}

							if(photons->GetApparatus(0) == GTreeRawEvent::APPARATUS_TAPS){
	
								Check_TAPS_TOF_photons_3ped->Fill(((photons->GetTime(0)-tagger->GetTagged_t(j))/photons->Particle(0).Vect().Mag())+1/0.299792458,photons->Particle(0).E(),-0.075);
							}
	
							if(photons->GetApparatus(1) == GTreeRawEvent::APPARATUS_TAPS){
	
								Check_TAPS_TOF_photons_3ped->Fill(((photons->GetTime(1)-tagger->GetTagged_t(j))/photons->Particle(1).Vect().Mag())+1/0.299792458,photons->Particle(1).E(),-0.075);
							}
		
							cosverteilung_collerated_proton->Fill(cospi0,beamphoton1E,-0.075);

							if(planesetting=="PARA"){
// 							if(linpol->GetPolPlane()==0){
			
								if(pt>0 ){
									kristallminus_targetplus_collerated_proton->Fill(phimeson,cospi0, beamphoton1E,-0.075);
									kristallminus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpb);
									kristallminus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpt);
								}
			
								if(pt<0 ){
									kristallminus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,-0.075);
									kristallminus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpb);
									kristallminus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpt);
								}
							}

							if(planesetting=="PERP"){
// 							if(linpol->GetPolPlane()==1){
			
								if(pt>0 ){
									kristallplus_targetplus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,-0.075);
									kristallplus_targetplus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpb);
									kristallplus_targetplus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpt);;
								}
				
								if(pt<0 ){
									kristallplus_targetminus_collerated_proton->Fill(phimeson,cospi0,beamphoton1E,-0.075);
									kristallplus_targetminus_collerated_pb_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpb);
									kristallplus_targetminus_collerated_pt_proton->Fill(phimeson,cospi0,beamphoton1E,scalbackgroundpt);
								}
							}						
						}
	
					}//end sideband
// 				
				}
// 				//all events with proton identified

		}//tagger
}


void PPi0Analyse::fOnEndProcessing() {

}


void	PPi0Analyse::ProcessScalerRead()
{
    //time.ScalerReadCorrection(5);
}
