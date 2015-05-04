


void	ShowHistogram(const TH1D* data, const TH1D* mcSignal, const TH1D* mcBG)
{
	if(!data)
	{
		std::cout << "Can not open Histogram." << std::endl;
		return;
	}
	Double_t	max	= data->GetMaximum();
	
	TH1D*	both	= mcSignal->Clone();
	both->Add(mcBG);
	Double_t	scale	= max/both->GetMaximum();
	both->Scale(scale);
	both->SetLineColor(kYellow);
	both->Draw();
	
	mcSignal->Scale(scale);
	mcSignal->SetLineColor(kMagenta);
	mcSignal->Draw("SAME");
	
	mcBG->Scale(scale);
	mcBG->SetLineColor(kGreen);
	mcBG->Draw("SAME");
	
	data->SetLineColor(kBlack);
	data->Draw("SAME");
}


void	OpenHistogram(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TString& str, const double signalScale)
{
	TH1D*	data		= (TH1D*)dataFile->Get(str);
	TH1D*	mcSignal	= (TH1D*)mcSignalFile->Get(str);
	mcSignal->Scale(signalScale);
	TH1D*	mcBG		= (TH1D*)mcBGFile->Get(str);
	ShowHistogram(data, mcSignal, mcBG);
}

void	ShowPhysicsCheckProton(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	TCanvas*		can	= new TCanvas("CanCheckProton", "CheckProton", 1500, 800);
	can->Divide(4,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_prAng", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_copl", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_prAng", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_copl", signalScale);
	can->cd(3);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_prAng", signalScale);
	can->cd(7);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_copl", signalScale);
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_prAng", signalScale);
	can->cd(8);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_copl", signalScale);
	
	out->cd();
	can->Write();
}
void	ShowPhysicsRaw(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	can	= new TCanvas("CanRaw", "Raw", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_IM", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_MM", signalScale);
	can->cd(3)->SetLogz();
	TH2D*	TOF	= (TH2D*)dataFile->Get("WithProton/Raw/_Raw_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub0IM", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub1IM", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub2IM", signalScale);
	
	out->cd();
	can->Write();
}
void	ShowPhysicsSubIM(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	can	= new TCanvas("CanSubIM", "SubIM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_IM", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_MM", signalScale);
	can->cd(3)->SetLogz();
	TOF	= (TH2D*)dataFile->Get("WithProton/SubIM_Cut/_subIMCut_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub0IM", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub1IM", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub2IM", signalScale);
	
	out->cd();
	can->Write();
}
void	ShowPhysicsMM(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	can	= new TCanvas("CanMM", "MM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_IM", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_MM", signalScale);
	can->cd(3)->SetLogz();
	TOF	= (TH2D*)dataFile->Get("WithProton/MM_Cut/_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub0IM", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub1IM", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub2IM", signalScale);
	
	out->cd();
	can->Write();
}

void	ShowPhysics(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	ShowPhysicsCheckProton(dataFile, mcSignalFile, mcBGFile, out, signalScale);
	ShowPhysicsRaw(dataFile, mcSignalFile, mcBGFile, out, signalScale);
	ShowPhysicsSubIM(dataFile, mcSignalFile, mcBGFile, out, signalScale);
	ShowPhysicsMM(dataFile, mcSignalFile, mcBGFile, out, signalScale);
}

void	ShowPhysics(const char* dataFileName, const char* mcSignalFileName, const char* mcBGFileName, const double signalScale)
{
	TFile*	dataFile		= TFile::Open(dataFileName);
	if(!dataFile)
	{
		std::cout << "Can not open dataFile " << dataFileName << std::endl;
		return;
	}
	TFile*	mcSignalFile	= TFile::Open(mcSignalFileName);
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << mcSignalFileName << std::endl;
		return;
	}
	TFile*	mcBGFile		= TFile::Open(mcBGFileName);
	if(!mcBGFile)
	{
		std::cout << "Can not open mcBGFile " << mcBGFileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("result.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	ShowPhysics(dataFile, mcSignalFile, mcBGFile, out, signalScale);
}

double	SignalScale(const TFile* mcSignalFile, const TFile* mcBGFile)
{
	TTree*	treeSignal	= (TTree*)mcSignalFile->Get("tagger");
	TTree*	treeBG		= (TTree*)mcBGFile->Get("tagger");
	
	treeSignal->Draw("taggedChannel>>hSignal");
	TH1D*	hSignal = (TH1D*)gDirectory->Get("hSignal");
	cout << hSignal->GetBinContent(1) << endl;
	
	treeBG->Draw("taggedChannel>>hBG");
	TH1D*	hBG = (TH1D*)gDirectory->Get("hBG");
	cout << hBG->GetBinContent(1) << endl;
	
	cout << (0.082*hBG->GetBinContent(1))/(3*hSignal->GetBinContent(1)) << endl;
	return (0.082*hBG->GetBinContent(1))/(3*hSignal->GetBinContent(1));
}

void	ShowPhysics(const char* dir = ".")
{
	std::strstream	acquSignalFileName;
	std::strstream	acquBGFileName;
	
	acquSignalFileName << dir << "/Acqu_g4_sim_etap_pi0pi0eta_00.root";
	acquBGFileName << dir << "/Acqu_g4_sim_pi0pi0pi0_6g_00.root";
	TFile*	mcSignalFile	= TFile::Open(acquSignalFileName.str().c_str());
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << acquSignalFileName << std::endl;
		return;
	}
	TFile*	mcBGFile		= TFile::Open(acquBGFileName.str().c_str());
	if(!mcBGFile)
	{
		std::cout << "Can not open mcBGFile " << acquBGFileName << std::endl;
		return;
	}
	double	signalScale	= SignalScale(mcSignalFile, mcBGFile);
	cout << "signalScale:   " << signalScale << endl;
	
	
	std::strstream	dataFileName;
	std::strstream	mcSignalFileName;
	std::strstream	mcBGFileName;
	
	dataFileName << dir << "/Result_CB.root";
	mcSignalFileName << dir << "/Phys_g4_sim_etap_pi0pi0eta_00.root";
	mcBGFileName << dir << "/Phys_g4_sim_pi0pi0pi0_6g_00.root";
	
	cout << "FitBins:" << endl;
	ShowPhysics(dataFileName.str().c_str(), mcSignalFileName.str().c_str(), mcBGFileName.str().c_str(), signalScale);
}
