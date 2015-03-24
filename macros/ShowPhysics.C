


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


void	OpenHistogram(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TString& str)
{
	TH1D*	data		= (TH1D*)dataFile->Get(str);
	TH1D*	mcSignal	= (TH1D*)mcSignalFile->Get(str);
	TH1D*	mcBG		= (TH1D*)mcBGFile->Get(str);
	ShowHistogram(data, mcSignal, mcBG);
}

void	ShowPhysicsCheckProton(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{
	TCanvas*		can	= new TCanvas("CanCheckProton", "CheckProton", 1500, 800);
	can->Divide(4,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_prAng");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_copl");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_prAng");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_copl");
	can->cd(3);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_prAng");
	can->cd(7);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_copl");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_prAng");
	can->cd(8);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_copl");
	
	out->cd();
	can->Write();
}
void	ShowPhysicsRaw(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{
	can	= new TCanvas("CanRaw", "Raw", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_MM");
	can->cd(3)->SetLogz();
	TH2D*	TOF	= (TH2D*)dataFile->Get("WithProton/Raw/_Raw_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub2IM");
	
	out->cd();
	can->Write();
}
void	ShowPhysicsSubIM(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{
	can	= new TCanvas("CanSubIM", "SubIM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_MM");
	can->cd(3)->SetLogz();
	TOF	= (TH2D*)dataFile->Get("WithProton/SubIM_Cut/_subIMCut_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub2IM");
	
	out->cd();
	can->Write();
}
void	ShowPhysicsMM(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{
	can	= new TCanvas("CanMM", "MM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_MM");
	can->cd(3)->SetLogz();
	TOF	= (TH2D*)dataFile->Get("WithProton/MM_Cut/_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub2IM");
	
	out->cd();
	can->Write();
}

void	ShowPhysics(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{
	ShowPhysicsCheckProton(dataFile, mcSignalFile, mcBGFile, out);
	ShowPhysicsRaw(dataFile, mcSignalFile, mcBGFile, out);
	ShowPhysicsSubIM(dataFile, mcSignalFile, mcBGFile, out);
	ShowPhysicsMM(dataFile, mcSignalFile, mcBGFile, out);
}

void	ShowPhysics(const char* dataFileName, const char* mcSignalFileName, const char* mcBGFileName)
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
	
	ShowPhysics(dataFile, mcSignalFile, mcBGFile, out);
}
