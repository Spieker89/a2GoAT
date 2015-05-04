
void	FindCuts(const TH1D* data, const double sigmaCount, double& min, double& max, double mean = 135)
{
	TF1*	fit = new TF1("fitfkt", "gaus(0)", mean -100, mean+100);
	fit->SetParameters(data->GetMaximum(), mean, 10);
	
	data->GetXaxis()->SetTitle("IM #gamma#gamma [MeV]");	
	data->Fit(fit, "QR0");
	data->SetStats(0);	
    data->SetAxisRange(mean -100, mean+100);	
	data->Draw("SAME");
	fit->Draw("SAME");
	
	min	= fit->GetParameter(1) - (sigmaCount*fit->GetParameter(2));
	max	= fit->GetParameter(1) + (sigmaCount*fit->GetParameter(2));
	
	TLine*	minLine	= new TLine(min, 0, min, data->GetMaximum());
	TLine*	maxLine	= new TLine(max, 0, max, data->GetMaximum());
	
	minLine->Draw("SAME");
	maxLine->Draw("SAME");
}

void	FindCuts(const TFile* mcSignalFile, const FILE* out, const double sigmaCount)
{        
	can	= new TCanvas("canFindCuts", "FindCuts", 1500, 800);
    can->Divide(4,2);
    
    double	min[4], max[4];
    
    can->cd(1);
    FindCuts((TH1D*)mcSignalFile->Get("WithoutProton/Raw/_Raw_sub0IM"), sigmaCount, min[0], max[0], 547);
    cout << "SubIm Eta: " << min << " - " << max << endl;
    
    can->cd(2);
    FindCuts((TH1D*)mcSignalFile->Get("WithoutProton/Raw/_Raw_sub1IM"), sigmaCount, min[1], max[1]);
    cout << "SubIm Pi0a: " << min << " - " << max << endl;
    
    can->cd(3);
    FindCuts((TH1D*)mcSignalFile->Get("WithoutProton/Raw/_Raw_sub2IM"), sigmaCount, min[2], max[2]);
    cout << "SubIm Pi0b: " << min << " - " << max << endl;
    
    
    can->cd(4);
    FindCuts((TH1D*)mcSignalFile->Get("WithoutProton/SubIM_Cut/_subIMCut_MM"), sigmaCount, min[3], max[3], 938);
    cout << "MM: " << min << " - " << max << endl;
    
    fprintf(out, "Cut-Etap-SubIM:                 %lf %lf %lf %lf %lf %lf\n", min[0], max[0], min[1], max[1], min[2], max[2]);
    fprintf(out, "Cut-Etap-MM:                    %lf %lf\n", min[3], max[3]);
    fflush(out);
    
    
    can->cd(5);
    FindCuts((TH1D*)mcSignalFile->Get("WithProton/Raw/_Raw_sub0IM"), sigmaCount, min[0], max[0], 547);
    cout << "SubIm Eta: " << min << " - " << max << endl;
    
    can->cd(6);
    FindCuts((TH1D*)mcSignalFile->Get("WithProton/Raw/_Raw_sub1IM"), sigmaCount, min[1], max[1]);
    cout << "SubIm Pi0a: " << min << " - " << max << endl;
    
    can->cd(7);
    FindCuts((TH1D*)mcSignalFile->Get("WithProton/Raw/_Raw_sub2IM"), sigmaCount, min[2], max[2]);
    cout << "SubIm Pi0b: " << min << " - " << max << endl;
    
    
    can->cd(8);
    FindCuts((TH1D*)mcSignalFile->Get("WithProton/SubIM_Cut/_subIMCut_MM"), sigmaCount, min[3], max[3], 938);
    cout << "MM: " << min << " - " << max << endl;
    
    fprintf(out, "Cut-Etap-Proton-SubIM:          %lf %lf %lf %lf %lf %lf\n", min[0], max[0], min[1], max[1], min[2], max[2]);
    fprintf(out, "Cut-Etap-Proton-MM:             %lf %lf\n", min[3], max[3]);
    fflush(out);
}	

void	FindCuts(const char* mcSignalFileName, const double sigmaCount)
{
	TFile*	mcSignalFile	= TFile::Open(mcSignalFileName);
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << mcSignalFileName << std::endl;
		return;
	}
	FILE*	out				= fopen("configfiles/MyPhysics.dat", "w");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	FindCuts(mcSignalFile, out, sigmaCount);
}
