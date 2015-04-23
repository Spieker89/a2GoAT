
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
    
	TH2D*	hist2	= (TH2D*)file->Get("CalibCB");
	TH1D*	hist	= hist2->ProjectionX();
	if(!hist)
	{
		std::cout << "Can not open hist CalibCB." << std::endl;
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
    
	hist->Fit(fit, "R");
    hist->SetAxisRange(100, 170);
	hist->Draw();
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
    
	TH2D*	hist2	= (TH2D*)file->Get("CalibCB");
	TH1D*	hist	= hist2->ProjectionX();
	if(!hist)
	{
		std::cout << "Can not open hist CalibCB." << std::endl;
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
