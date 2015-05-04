
void	FitScaleMC(const char* fileName)
{
	TFile*	file		= TFile::Open(fileName);
	if(!file)
	{
		std::cout << "Can not open file " << fileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("resultFitScaleMC.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file resultFitScaleMC.root" << std::endl;
		return;
	}	    
    
    TCanvas*	can	= new TCanvas("FitBinsMain", "FitBinsMain", 1500, 800);
	can->Divide(2, 1);
	
	TH2D*	hist2	= (TH2D*)file->Get("CalibCBCorr");
	TH1D*	hist	= hist2->ProjectionX()->Clone();
	TH1D*	histEta	= hist2->ProjectionX()->Clone();
	TH2D*	hist2uc		= (TH2D*)file->Get("CalibCB");
	TH1D*	histuc		= hist2uc->ProjectionX()->Clone();
	TH1D*	histEtauc	= hist2uc->ProjectionX()->Clone();
	if(!hist)
	{
		std::cout << "Can not open hist CalibCBCorr." << std::endl;
		return;
	}
	
	TF1*	fit = new TF1("fitfkt", "gaus(0)+pol2(3)", 100, 170);
	fit->SetParameters(hist->GetMaximum(), 135, 10, hist->GetMaximum()/3, 0, 0);
	fit->SetParLimits(0, hist->GetMaximum()/10, hist->GetMaximum()*2);
	fit->SetParLimits(1, 130, 140);
	fit->SetParLimits(2, 5, 25);
	fit->SetParLimits(3, 0, hist->GetMaximum());
	fit->SetParLimits(4, -100, 100);
	fit->SetParLimits(5, -10, 10);
    
	hist->Fit(fit, "R0");
    hist->SetAxisRange(100, 170);
    can->cd(1);
	hist->GetXaxis()->SetTitle("IM #gamma#gamma [MeV]");
	hist->Draw();
	fit->Draw("SAME");
    histuc->SetAxisRange(100, 170);
    histuc->SetLineColor(kMagenta); 
	histuc->Draw("SAME");
    histEta->SetAxisRange(500, 600);
    
	TF1*	fitEta = new TF1("fitfktEta", "gaus(0)+pol1(3)", 500, 600);
	fitEta->SetParameters(histEta->GetMaximum(), 547, 10, histEta->GetMaximum()/3, 0);
	fitEta->SetParLimits(0, histEta->GetMaximum()/10, histEta->GetMaximum()*2);
	fitEta->SetParLimits(1, 530, 570);
	fitEta->SetParLimits(2, 5, 50);
	fitEta->SetParLimits(3, 0, histEta->GetMaximum());
	fitEta->SetParLimits(4, -100, 100);

    histEta->SetAxisRange(500, 600);
	histEta->Fit(fitEta, "R0");
    can->cd(2);
	histEta->GetXaxis()->SetTitle("IM #gamma#gamma [MeV]");
	histEta->Draw();
	fitEta->Draw("SAME");
    histEtauc->SetLineColor(kMagenta);
	histEtauc->Draw("SAME");
}

void	FitScaleMCEta(const char* fileName)
{
	TFile*	file		= TFile::Open(fileName);
	if(!file)
	{
		std::cout << "Can not open file " << fileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("resultFitScaleMC.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file resultFitScaleMC.root" << std::endl;
		return;
	}	
    
	TH2D*	hist2	= (TH2D*)file->Get("CalibCBCorr");
	TH1D*	hist	= hist2->ProjectionX();
	if(!hist)
	{
		std::cout << "Can not open hist CalibCBCorr." << std::endl;
		return;
	}
	
	TF1*	fitEta = new TF1("fitfktEta", "gaus(0)+pol1(3)", 500, 600);
	fitEta->SetParameters(hist->GetMaximum(), 547, 10, hist->GetMaximum()/3, 0);
	fitEta->SetParLimits(0, hist->GetMaximum()/10, hist->GetMaximum()*2);
	fitEta->SetParLimits(1, 530, 570);
	fitEta->SetParLimits(2, 5, 50);
	fitEta->SetParLimits(3, 0, hist->GetMaximum());
	fitEta->SetParLimits(4, -100, 100);

    hist->SetAxisRange(500, 600);
	hist->Fit(fitEta, "R");
	hist->Draw();
}
