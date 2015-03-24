
void	FitBinsFitAll(const TH1D* h, Double_t* par)
{
	TF1*	fit = new TF1("fitfkt", "gaus(0)+gaus(3)", 750, 1100);
	fit->SetParameters(h->GetMaximum(), 957, 10, h->GetMaximum()/20, 950, 50);
	fit->SetParLimits(0, h->GetMaximum()/10, h->GetMaximum());
	fit->SetParLimits(1, 940, 990);
	fit->SetParLimits(2, 5, 25);
	fit->SetParLimits(3, 20, h->GetMaximum()/2);
	fit->SetParLimits(4, 700, 1200);
	fit->SetParLimits(5, 25, 100);

	h->Fit(fit, "R");
    h->SetAxisRange(750, 1100);
	h->Draw();
	
	par[0] = fit->GetParameter(0);
	par[1] = fit->GetParameter(1);
	par[2] = fit->GetParameter(2);
	par[3] = fit->GetParameter(3);
	par[4] = fit->GetParameter(4);
	par[5] = fit->GetParameter(5);
}

void	FitBinsFit(const TH1D* h, const Double_t* par)
{
	TF1*	fit = new TF1("fitfkt", "gaus(0)+gaus(3)", 750, 1100);
	fit->SetParameters(h->GetMaximum(), par[1], par[2], h->GetMaximum()/20, par[4], par[5]);
	fit->SetParLimits(0, 0, h->GetMaximum());
	fit->SetParLimits(1, par[1]-10, par[1]+10);
	fit->SetParLimits(2, par[2]/2, par[2]*2);
	fit->SetParLimits(3, 0, h->GetMaximum()/2);
	fit->SetParLimits(4, par[4]-25, par[4]+50);
	fit->SetParLimits(5, par[5], par[5]*2);

	h->Fit(fit, "R");
    h->SetAxisRange(750, 1100);
	h->Draw();
}

void	FitBins(const TFile* dataFile, const TFile* out, const Int_t addedChannels)
{
	std::vector<int>	ch;
	std::vector<int>	size;
	while(addedChannels*ch.size()<48)
	{
		if(48-(addedChannels*ch.size())>=addedChannels)
			size.push_back(addedChannels);
		else
			size.push_back(48-(addedChannels*ch.size()));
		ch.push_back(addedChannels*ch.size());
		std::cout << "starting number: " << ch.back() << "     channel count: " << size.back() << std::endl;
	}
	
	Double_t	par[6];
	TH2D*	data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/Final/TaggerBinning/BackgroundSubstraction/Final_IM_Bins_Prompt");
	FitBinsFitAll((TH1D*)data->ProjectionX(TString("BinUpTo35").Data(), 0, 35), par);
	
	can	= new TCanvas("FitBins", "FitBins", 1500, 800);
	can->Divide(4,TMath::Ceil(double(ch.size())/4));
	
	for(int i=0; i<ch.size(); i++)
	{
		can->cd(i+1);
		FitBinsFit((TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1), par);
	}
	out->cd();
	can->Write();
}

void	FitBins(const char* dataFileName, const Int_t addedChannels)
{
	TFile*	dataFile		= TFile::Open(dataFileName);
	if(!dataFile)
	{
		std::cout << "Can not open dataFile " << dataFileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("result.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	FitBins(dataFile, out, addedChannels);
}
