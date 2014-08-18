
void	FitGauss(TH1* hist)
{
	TF1* ff = new TF1(TString(hist->GetName()).Append("f1"), "gaus");
	ff->SetParameters(hist->GetMaximum(), hist->GetMaximumBin(), 0);
    hist->Fit(TString(hist->GetName()).Append("f1"), "QN");
    cout << "FitGaus " << hist->GetName() << " µ: " << ff->GetParameter(1) << "   sigma: " << ff->GetParameter(2) << endl;
}

void	FitGaussGauss(TH1* hist)
{
	TF1* ff = new TF1(TString(hist->GetName()).Append("f1"), "gaus(0)+gaus(3)");
	ff->SetParameters(hist->GetMaximum()/2, 960, 10, hist->GetMaximum()/10, 950, 70);
	ff->SetParLimits(1, 930, 990);
	ff->SetParLimits(2, 1, 40);
	ff->SetParLimits(4, 900, 1100);
	ff->SetParLimits(5, 30, 200);
    hist->Fit(TString(hist->GetName()).Append("f1"), "Q");
    cout << "FitGausGaus   " << hist->GetName() << " µ: " << ff->GetParameter(1) << "   sigma: " << ff->GetParameter(2) << endl;
    cout << "   Background " << hist->GetName() << " µ: " << ff->GetParameter(4) << "   sigma: " << ff->GetParameter(5) << endl;
}


void	FindCutsInMC(const char* filename)
{
	TFile*	file	= TFile::Open(filename);
	if(!file)
	{
		cout << "could not open file " << filename << endl;
		return;
	}
	
	TFile*	out	= TFile::Open("resultFindCutsInMC.root", "RECREATE");
	if(!out)
	{
		cout << "could not create output file resultFindCutsInMC" << endl;
		return;
	}
	
	// eta
	
	TDirectory*	eta	= file->GetDirectory("eta");
	if(eta)
	{
		// Without proton
		TDirectory*	WithoutProton	= eta->GetDirectory("WithoutProton");
		if(WithoutProton)
		{
			TDirectory*	Raw	= WithoutProton->GetDirectory("Raw");
			if(Raw)
			{
				TH1* hist = (TH1*)Raw->Get("eta_Raw_sub0im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("eta_Raw_sub1im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("eta_Raw_sub2im");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory Raw" << endl;
				
			TDirectory*	IM_Cut	= WithoutProton->GetDirectory("IM_Cut");
			if(IM_Cut)
			{
				TH1* hist = (TH1*)IM_Cut->Get("eta_SubImCut_mm");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory IM_Cut" << endl;
				
			TDirectory*	Fit	= WithoutProton->GetDirectory("Fit");
			if(Fit)
			{
				TDirectory*	fit3	= Fit->GetDirectory("fit3");
				if(fit3)
				{
					TDirectory*	CutConfidenceLevel	= fit3->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("eta_fit_fit3CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit3" << endl;
				
				TDirectory*	fit4	= Fit->GetDirectory("fit4");
				if(fit4)
				{
					TDirectory*	CutConfidenceLevel	= fit4->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("eta_fit_fit4CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit4" << endl;
			}
			else
				cout << "could not open directory Fit" << endl;
		}
		else
			cout << "could not open directory WithoutProton" << endl;
			
		// With proton
		TDirectory*	WithProton	= eta->GetDirectory("WithProton");
		if(WithProton)
		{
			TDirectory*	Raw	= WithProton->GetDirectory("Raw");
			if(Raw)
			{
				TH1* hist = (TH1*)Raw->Get("eta_proton_Raw_sub0im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("eta_proton_Raw_sub1im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("eta_proton_Raw_sub2im");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory Raw" << endl;
				
			TDirectory*	IM_Cut	= WithProton->GetDirectory("IM_Cut");
			if(IM_Cut)
			{
				TH1* hist = (TH1*)IM_Cut->Get("eta_proton_SubImCut_mm");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory IM_Cut" << endl;
				
			TDirectory*	Fit	= WithProton->GetDirectory("Fit");
			if(Fit)
			{
				TDirectory*	fit3	= Fit->GetDirectory("fit3");
				if(fit3)
				{
					TDirectory*	CutConfidenceLevel	= fit3->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("eta_proton_fit_fit3CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit3" << endl;
				
				TDirectory*	fit4	= Fit->GetDirectory("fit4");
				if(fit4)
				{
					TDirectory*	CutConfidenceLevel	= fit4->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("eta_proton_fit_fit4CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit4" << endl;
			}
			else
				cout << "could not open directory Fit" << endl;
		}
		else
			cout << "could not open directory WithProton" << endl;
	}
	else
		cout << "could not open directory eta" << endl;
		
		
		
		
		
	
	// etap
	
	TDirectory*	etap	= file->GetDirectory("etap");
	if(etap)
	{
		// Without proton
		TDirectory*	WithoutProton	= etap->GetDirectory("WithoutProton");
		if(WithoutProton)
		{
			TDirectory*	Raw	= WithoutProton->GetDirectory("Raw");
			if(Raw)
			{
				TH1* hist = (TH1*)Raw->Get("etap_Raw_sub0im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("etap_Raw_sub1im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("etap_Raw_sub2im");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory Raw" << endl;
				
			TDirectory*	IM_Cut	= WithoutProton->GetDirectory("IM_Cut");
			if(IM_Cut)
			{
				TH1* hist = (TH1*)IM_Cut->Get("etap_SubImCut_mm");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory IM_Cut" << endl;
				
			TDirectory*	Fit	= WithoutProton->GetDirectory("Fit");
			if(Fit)
			{
				TDirectory*	fit3	= Fit->GetDirectory("fit3");
				if(fit3)
				{
					TDirectory*	CutConfidenceLevel	= fit3->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("etap_fit_fit3CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit3" << endl;
				
				TDirectory*	fit4	= Fit->GetDirectory("fit4");
				if(fit4)
				{
					TDirectory*	CutConfidenceLevel	= fit4->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("etap_fit_fit4CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit4" << endl;
			}
			else
				cout << "could not open directory Fit" << endl;
		}
		else
			cout << "could not open directory WithoutProton" << endl;
			
		// With proton
		TDirectory*	WithProton	= etap->GetDirectory("WithProton");
		if(WithProton)
		{
			TDirectory*	Raw	= WithProton->GetDirectory("Raw");
			if(Raw)
			{
				TH1* hist = (TH1*)Raw->Get("etap_proton_Raw_sub0im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("etap_proton_Raw_sub1im");
				FitGauss(hist);
				out->cd();
				hist->Write();
				
				hist = (TH1*)Raw->Get("etap_proton_Raw_sub2im");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory Raw" << endl;
				
			TDirectory*	IM_Cut	= WithProton->GetDirectory("IM_Cut");
			if(IM_Cut)
			{
				TH1* hist = (TH1*)IM_Cut->Get("etap_proton_SubImCut_mm");
				FitGauss(hist);
				out->cd();
				hist->Write();
			}
			else
				cout << "could not open directory IM_Cut" << endl;
				
			TDirectory*	Fit	= WithProton->GetDirectory("Fit");
			if(Fit)
			{
				TDirectory*	fit3	= Fit->GetDirectory("fit3");
				if(fit3)
				{
					TDirectory*	CutConfidenceLevel	= fit3->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("etap_proton_fit_fit3CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit3" << endl;
				
				TDirectory*	fit4	= Fit->GetDirectory("fit4");
				if(fit4)
				{
					TDirectory*	CutConfidenceLevel	= fit4->GetDirectory("CutConfidenceLevel");
					if(CutConfidenceLevel)
					{
						TH1* hist = (TH1*)CutConfidenceLevel->Get("etap_proton_fit_fit4CutCL");
						FitGaussGauss(hist);
						out->cd();
						hist->Write();
					}
					else
						cout << "could not open directory CutConfidenceLevel" << endl;
				}
				else
					cout << "could not open directory fit4" << endl;
			}
			else
				cout << "could not open directory Fit" << endl;
		}
		else
			cout << "could not open directory WithProton" << endl;
	}
	else
		cout << "could not open directory eta" << endl;
}
