


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

void	FitProfile(const TH2D* data, const int i, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& etapWidth, Double_t& BGWidth)
{
	TH1D*	slice		= (TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
	slice->SetAxisRange(850, 1020);
	slice->Draw();
	
	nEtapCut	= 0;
	for(int i=1; i<=slice->GetNbinsX(); i++)
	{
		if(slice->GetBinCenter(i)>(957-8.83) && slice->GetBinCenter(i)<(957+8.83))
			nEtapCut	+= slice->GetBinContent(i);
	}
	
	Double_t	max	= slice->GetMaximum();
	
	TF1*	fit = new TF1("fit","gaus(0) + gaus(3)",850,1020);
	fit->SetParameters(max, 958, etapWidth, max/10, 930, BGWidth);
	fit->SetParLimits(0, 0, max *2);
	fit->SetParLimits(1, 950, 970);
	fit->SetParLimits(2, etapWidth, etapWidth);
	fit->SetParLimits(3, 0, max/2);
	fit->SetParLimits(4, 900, 980);
	fit->SetParLimits(5, BGWidth, BGWidth);
	
	slice->Fit(fit, "R");
	
	nEtapFit	= fit->GetParameter(0) * fit->GetParameter(2) * 1.7724;
	dNEtapFit	= fit->GetParError(0) * fit->GetParameter(2) * 1.7724;
}

void	FitProfileSim(const TH2D* data, const int i, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& etapWidth, Double_t& BGWidth)
{
	TH1D*	slice		= (TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
	slice->SetAxisRange(850, 1020);
	slice->Draw();
	
	nEtapCut	= 0;
	for(int i=1; i<=slice->GetNbinsX(); i++)
	{
		if(slice->GetBinCenter(i)>(957-8.83) && slice->GetBinCenter(i)<(957+8.83))
			nEtapCut	+= slice->GetBinContent(i);
	}
	
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
	
	nEtapFit	= fit->GetParameter(0) * fit->GetParameter(2) * 1.7724;
	dNEtapFit	= fit->GetParError(0) * fit->GetParameter(2) * 1.7724;
	
	etapWidth	= fit->GetParameter(2);
	BGWidth		= fit->GetParameter(5);
}

void	ReconstructionEff(const TH2D* mcSignal, TH1D* RecEff)
{
	RecEff->Reset();
	TH1D*	RecEffHelp	= new TH1D("help", "help", 48, 0, 48);
	TH1D*	slice;
	for(int i=0; i<48; i++)
	{
		slice		= (TH1D*)mcSignal->ProjectionX(TString("ReconstructionEff_Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
		RecEff->SetBinContent(i+1, slice->GetEntries());
		if(i<40)
			RecEffHelp->SetBinContent(i+1, 84910.0/40.0);
		else
		{
			RecEff->SetBinContent(i+1, 1);
			RecEffHelp->SetBinContent(i+1, 1);
		}	
	}
	RecEff->Divide(RecEffHelp);
	RecEffHelp->Scale(40.0/84910.0);
	RecEffHelp->SetLineColor(kMagenta);
	RecEffHelp->Draw();
	RecEff->Draw("Same");
}

void TaggEff(Double_t* array, Double_t* darray)
{
	array[0]	= 48.0/100.0;
	array[1]	= 47.0/100.0;
	array[2]	= 47.0/100.0;
	array[3]	= 49.0/100.0;
	array[4]	= 51.0/100.0;
	array[5]	= 53.0/100.0;
	array[6]	= 54.0/100.0;
	array[7]	= 55.0/100.0;
	array[8]	= 56.0/100.0;
	array[9]	= 56.0/100.0;
	array[10]	= 58.0/100.0;
	array[11]	= 57.0/100.0;
	array[12]	= 59.0/100.0;
	array[13]	= 59.0/100.0;
	array[14]	= 60.0/100.0;
	array[15]	= 58.0/100.0;
	array[16]	= 61.0/100.0;
	array[17]	= 60.0/100.0;
	array[18]	= 60.0/100.0;
	array[19]	= 61.0/100.0;
	array[20]	= 63.0/100.0;
	array[21]	= 62.0/100.0;
	array[22]	= 64.0/100.0;
	array[23]	= 64.0/100.0;
	array[24]	= 66.0/100.0;
	array[25]	= 62.0/100.0;
	array[26]	= 65.0/100.0;
	array[27]	= 62.0/100.0;
	array[28]	= 66.0/100.0;
	array[29]	= 63.0/100.0;
	array[30]	= 65.0/100.0;
	array[31]	= 63.0/100.0;
	array[32]	= 65.0/100.0;
	array[33]	= 65.0/100.0;
	array[34]	= 65.0/100.0;
	array[35]	= 63.0/100.0;
	array[36]	= 66.0/100.0;
	array[37]	= 62.0/100.0;
	array[38]	= 68.0/100.0;
	array[39]	= 65.0/100.0;
	array[40]	= 66.0/100.0;
	array[41]	= 65.0/100.0;
	array[42]	= 65.0/100.0;
	array[43]	= 64.0/100.0;
	array[44]	= 66.0/100.0;
	array[45]	= 62.0/100.0;
	array[46]	= 66.0/100.0;
	array[47]	= 0;
	for(int i=0; i<48; i++)
		darray[i]	= 2.0/100.0;
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
	
	can	= new TCanvas("CanFitBins_Sim", "FitBins_Sim", 1500, 800);
	can->Divide(8,5);
	Double_t	nEtapFit[48];
	Double_t	dNEtapFit[48];
	Double_t	nEtapCut[48];
	Double_t	dNEtapCut[48];
	Double_t	x[48];
	Double_t	dx[48];
	Double_t	etapWidth[48];
	Double_t	BGWidth[48];
	for(int i=0; i<48; i++)
	{
		can->cd(i+1);
		FitProfileSim(data, i, nEtapFit[i], dNEtapFit[i], nEtapCut[i], etapWidth[i], BGWidth[i]);
		dNEtapCut[i]	= sqrt(nEtapCut[i]);
		x[i]			= i;
		dx[i]			= 0;
		//std::cout << nEtapFit[i] << "   " << dNEtapFit[i] << "   " << nEtapCut[i] << "   " << dNEtapCut[i] << std::endl;
	}
	out->cd();
	can->Write();
	
	
	can	= new TCanvas("CanFitBins", "FitBins", 1500, 800);
	can->Divide(8,5);
	for(int i=0; i<48; i++)
	{
		can->cd(i+1);
		FitProfile(data, i, nEtapFit[i], dNEtapFit[i], nEtapCut[i], etapWidth[i], BGWidth[i]);
		dNEtapCut[i]	= sqrt(nEtapCut[i]);
		x[i]			= i;
		dx[i]			= 0;
		//std::cout << nEtapFit[i] << "   " << dNEtapFit[i] << "   " << nEtapCut[i] << "   " << dNEtapCut[i] << std::endl;
	}
	out->cd();
	can->Write();
	
	
	can	= new TCanvas("CanResult", "Result", 1500, 800);
	can->Divide(2,2);
	TGraphErrors*	graph = new TGraphErrors(48, x, nEtapFit, dx, dNEtapFit);
	can->cd(1);
	graph->Draw();
	graph = new TGraphErrors(48, x, nEtapCut, dx, dNEtapCut);
	can->cd(2);
	graph->Draw();
	data		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/TaggerBinning/_fit4_IM_Bins");
	TH1D*	RecEff	= new TH1D("RecEff", "RecEff", 48, 0, 48);
	can->cd(3);
	ReconstructionEff(data, RecEff);
	TH1D*	help	= (TH1D*)dataFile->Get("EPT_ScalerCorT");
	can->cd(4);
	help->Draw();
	
	can	= new TCanvas("CanTaggEff", "TaggEff", 1500, 800);
	can->Divide(2,0);
	Double_t	taggEff[48];
	Double_t	dTaggEff[48];
	TaggEff(taggEff, dTaggEff);
	graph = new TGraphErrors(48, x, taggEff, dx, dTaggEff);
	can->cd(1);
	graph->Draw();
	Double_t	sc[48];
	Double_t	dsc[48];
	for(int i=0; i<48; i++)
	{
		sc[i]			= help->GetBinContent(i+1) * taggEff[i];
		dsc[i]			= help->GetBinContent(i+1) * dTaggEff[i];
	}
	graph = new TGraphErrors(48, x, sc, dx, dsc);
	can->cd(2);
	graph->Draw();
	
	can	= new TCanvas("CanEndResult", "EndResult", 1500, 800);
	can->Divide(2,2);
	
	Double_t	y[48];
	Double_t	dy[48];
	for(int i=0; i<48; i++)
	{
		y[i]			= nEtapFit[i] / RecEff->GetBinContent(i+1);
		dy[i]			= dNEtapFit[i] / RecEff->GetBinContent(i+1);
	}
	graph = new TGraphErrors(48, x, y, dx, dy);
	can->cd(1);
	graph->Draw();
	for(int i=0; i<48; i++)
	{
		y[i]			= nEtapCut[i] / RecEff->GetBinContent(i+1);
		dy[i]			= dNEtapCut[i] / RecEff->GetBinContent(i+1);
	}
	graph = new TGraphErrors(48, x, y, dx, dy);
	can->cd(2);
	graph->Draw();
	
	
	//5.5e-29
	Double_t	y[47];
	Double_t	dy[47];
	for(int i=0; i<47; i++)
	{
		y[i]			= nEtapFit[i] / (RecEff->GetBinContent(i+1) * sc[i]);
		dy[i]			= dNEtapFit[i] / (RecEff->GetBinContent(i+1) * sc[i]);
	}
	graph = new TGraphErrors(47, x, y, dx, dy);
	can->cd(3);
	graph->Draw();
	for(int i=0; i<47; i++)
	{
		y[i]			= nEtapCut[i] * 3.8 / (RecEff->GetBinContent(i+1) * sc[i]);
		dy[i]			= dNEtapCut[i] * 3.8 / (RecEff->GetBinContent(i+1) * sc[i]);
	}
	graph = new TGraphErrors(47, x, y, dx, dy);
	can->cd(4);
	graph->Draw();
}
