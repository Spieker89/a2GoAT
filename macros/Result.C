


void	ShowHistogram(const TH1D* data, const TH1D* mcSignal, const TH1D* mcBG)
{
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


void	FitProfile(const TH2D* data, const int i)
{
	TH1D*	slice		= (TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
	slice->SetAxisRange(850, 1020);
	slice->Draw();
	
	Double_t	max	= slice->GetMaximum();
	
	TF1*	fit = new TF1("fit","gaus(0) + gaus(3)",850,1020);
	fit->SetParameters(max, 958, 10, max/10, 930, 35);
	fit->SetParLimits(0, max/10, max *2);
	fit->SetParLimits(1, 950, 970);
	fit->SetParLimits(2, 5, 25);
	fit->SetParLimits(3, 0, max/2);
	fit->SetParLimits(4, 900, 980);
	fit->SetParLimits(5, 25, 75);
	
	slice->Fit(fit, "R");
	return;
}

void	ReconstructionEff(const TH2D* mcSignal)
{
	TH1D*	RecEff	= new TH1D("RecEff", "RecEff", 48, 0, 48);
	TH1D*	RecEffHelp	= new TH1D("help", "help", 48, 0, 48);
	TH1D*	slice;
	for(int i=0; i<48; i++)
	{
		slice		= (TH1D*)mcSignal->ProjectionX(TString("ReconstructionEff_Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
		RecEff->SetBinContent(i+1, slice->GetEntries());
		RecEffHelp->SetBinContent(i+1, 84910/48);
	}
	RecEff->Divide(RecEffHelp);
	RecEff->Draw();
}

void	Result(const char* dataFileName, const char* mcSignalFileName, const char* mcBGFileName)
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
	OpenHistogram(dataFile, mcSignalFile, mcBGFile,  "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_copl");
	
	out->cd();
	can->Write();
	
	
	
	can	= new TCanvas("CanRaw", "Raw", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_MM");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub2IM");
	
	out->cd();
	can->Write();
	
	
	
	can	= new TCanvas("CanSubIM", "SubIM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_MM");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub2IM");
	
	out->cd();
	can->Write();



	can	= new TCanvas("CanMM", "MM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_MM");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub2IM");
	
	out->cd();
	can->Write();
	
	
	
	can	= new TCanvas("CanFit", "Fit", 1500, 800);
	can->Divide(2,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/_fit4_IM");
	can->cd(2)->SetLogz();
	TH2D*	data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/TaggerBinning/_fit4_IM_Bins");
	data->Draw("COL");
	can->cd(3);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/_fit4_ChiSq");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/_fit4_ConfLev");
	
	out->cd();
	can->Write();
	
	can	= new TCanvas("CanFitBins", "FitBins", 1500, 800);
	can->Divide(8,5);
	for(int i=0; i<48; i++)
	{
		can->cd(i+1);
		FitProfile(data, i);
	}
	
	out->cd();
	can->Write();
	
	data		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/TaggerBinning/_fit4_IM_Bins");
	TH1D*	RecEff	= new TH1D("RecEff", "RecEff", 48, 0, 48);
	
	can	= new TCanvas("CanEff", "Eff", 1500, 800);
	can->cd();
	ReconstructionEff(data);
}
