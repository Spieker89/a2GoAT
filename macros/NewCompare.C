#include "macros/NewBestFit.C"


void	CompareSimDataShowHistogram(const TH1D* data, const TH1D* mcSignal, const TH1D* mcBG, const double minimum, const double maximum, const char* title, const char* axisTitle)
{
	if(!data)
	{
		std::cout << "Can not open Histogram." << std::endl;
		return;
	}
	Double_t	max	= data->GetMaximum();
	
	Int_t	rebin = 3;
	
	TH1D*	both	= mcSignal->Clone();
	both->Add(mcBG);
	Double_t	scale	= max/both->GetMaximum();
	both->Scale(scale);
	both->SetLineColor(kYellow);
	cout << "here " << minimum << "   " << maximum << endl;
	both->Rebin(rebin);
	both->Scale(1.0/rebin);
	both->GetXaxis()->SetRangeUser(minimum, maximum);
	both->SetTitle(title);
	both->GetXaxis()->SetTitle(axisTitle);
	both->SetStats(0);
	both->Draw();
	
	mcSignal->Scale(scale);
	mcSignal->SetLineColor(kMagenta);
	mcSignal->Rebin(rebin);
	mcSignal->Scale(1.0/rebin);
	mcSignal->GetXaxis()->SetRangeUser(minimum, maximum);
	mcSignal->Draw("SAME");
	
	mcBG->Scale(scale);
	mcBG->SetLineColor(kGreen);
	mcBG->Rebin(rebin);
	mcBG->Scale(1.0/rebin);
	mcBG->GetXaxis()->SetRangeUser(minimum, maximum);
	mcBG->Draw("SAME");
	
	data->SetLineColor(kBlack);
	data->Rebin(rebin);
	data->Scale(1.0/rebin);
	data->GetXaxis()->SetRangeUser(minimum, maximum);
	data->Draw("SAME");
}


void	CompareSimDataOpenHistogram(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TString& str, const double signalScale, const double min, const double max, const char* title, const char* axisTitle)
{
	TH1D*	data		= (TH1D*)dataFile->Get(str);
	TH1D*	mcSignal	= (TH1D*)mcSignalFile->Get(str);
	mcSignal->Scale(signalScale);
	TH1D*	mcBG		= (TH1D*)mcBGFile->Get(str);
	CompareSimDataShowHistogram(data, mcSignal, mcBG, min, max, title, axisTitle);
}

void	CompareSimData(const TFile* dataFile, const TFile* mcSinalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	can	= new TCanvas("canCompareSimData", "CompareSimData", 1500, 800);
    can->Divide(1,3);
    
    can->cd(1);
	CompareSimDataOpenHistogram(dataFile, mcSinalFile, mcBGFile, "WithProton/fitProton6/fitProton6_IM", signalScale, 800.0, 1050.0, "inv mass 6#gamma raw", "IM 6#gamma [MeV]");
    can->cd(2);
	CompareSimDataOpenHistogram(dataFile, mcSinalFile, mcBGFile, "WithProton/fitProton6/IM", signalScale, 800.0, 1050.0, "inv mass 6#gamma fitted", "IM 6#gamma [MeV]");
	can->cd(3);
	CompareSimDataOpenHistogram(dataFile, mcSinalFile, mcBGFile, "WithProton/fitProton6/SubAll", signalScale, 0.0, 600.0, "inv mass all 2#gamma combinations raw", "IM 2#gamma [MeV]");
	
	out->cd();
	can->Write();
}


void	CompareSimData(const char* dataFileName, const char* mcSignalFileName, const char* mcBGFileName, const double	signalScale)
{
	TFile*	dataFile		= TFile::Open(dataFileName);
	if(!dataFile)
	{
		std::cout << "Can not open dataFile " << dataFileName << std::endl;
		return;
	}
	TFile*	mcSignalFile		= TFile::Open(mcSignalFileName);
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << mcSignalFileName << std::endl;
		return;
	}
	TFile*	mcBGFile			= TFile::Open(mcBGFileName);
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
	
	CompareSimData(dataFile, mcSignalFile, mcBGFile, out, signalScale);
}

void	CompareSimData(const char* dir = ".")
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
	double	signalScale	= BestFitSignalScale(mcSignalFile, mcBGFile);
	cout << "signalScale:   " << signalScale << endl;
	
	
	std::strstream	dataFileName;
	std::strstream	mcSignalFileName;
	std::strstream	mcBGFileName;
	
	dataFileName << dir << "/Result_CB.root";
	mcSignalFileName << dir << "/Phys_g4_sim_etap_pi0pi0eta_00.root";
	mcBGFileName << dir << "/Phys_g4_sim_pi0pi0pi0_6g_00.root";
	
	cout << "FitBins:" << endl;
	CompareSimData(dataFileName.str().c_str(), mcSignalFileName.str().c_str(), mcBGFileName.str().c_str(), signalScale);
}
