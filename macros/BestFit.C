

#include "macros/ShowPhysics.C"


void	BestFitCL(const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{
	can	= new TCanvas("BestFitCL", "BestFitCL", 1500, 800);
    can->Divide(4,2);
    can->cd(1);
    TH2D*	help		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit1/_ConfLev");
    help->Draw("COLZ");
    can->cd(2);
    help		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit3/_ConfLev");
    help->Draw("COLZ");
    can->cd(3);
    help		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/_ConfLev");
    help->Draw("COLZ");
    can->cd(4);
    help		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit7Proton/_ConfLev");
    help->Draw("COLZ");
    can->cd(5);
    help		= (TH2D*)mcBGFile->Get("WithProton/MM_Cut/fit1/_ConfLev");
    help->Draw("COLZ");
    can->cd(6);
    help		= (TH2D*)mcBGFile->Get("WithProton/MM_Cut/fit3/_ConfLev");
    help->Draw("COLZ");
    can->cd(7);
    help		= (TH2D*)mcBGFile->Get("WithProton/MM_Cut/fit4/_ConfLev");
    help->Draw("COLZ");
    can->cd(8);
    help		= (TH2D*)mcBGFile->Get("WithProton/MM_Cut/fit7Proton/_ConfLev");
    help->Draw("COLZ");
    
	out->cd();
	can->Write();
}

void	BestFitFit(const TH1D* hist1, const TH1D* hist2, Double_t* width)
{
	TH1D*	h = (TH1D*)hist1->Clone();
	h.Add(hist2);
	
	TF1*	fit = new TF1("fitfkt", "gaus(0)+gaus(3)", 750, 1100);
	fit->SetParameters(h->GetMaximum(), 957, 10, h->GetMaximum()/20, 950, 50);
	fit->SetParLimits(0, h->GetMaximum()/10, h->GetMaximum()*2);
	fit->SetParLimits(1, 940, 990);
	fit->SetParLimits(2, 5, 25);
	fit->SetParLimits(3, 0, h->GetMaximum()/2);
	fit->SetParLimits(4, 700, 1200);
	fit->SetParLimits(5, 25, 100);

	h->Fit(fit, "R");
    h->SetAxisRange(750, 1100);
	h->Draw();
	
	width[0] = fit->GetParameter(2);
	width[1] = fit->GetParError(2);
	width[2] = fit->GetParameter(5);
	width[3] = fit->GetParError(5);
}

void	BestFit(const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{	
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/Final/Final_IM");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit4/Final/Final_IM");
        help1->SetLineColor(kMagenta); 
		help1->SetAxisRange(850, 1050);       
		help1->SetStats(0);
		help1->GetXaxis()->SetTitle("inv. Mass (#eta#pi^{0}#pi^{0}) [MeV/c]");
		help1->SetTitle("");
        help1->Draw();    
		help2->SetStats(0);
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
        
	BestFitCL(mcSignalFile, mcBGFile, out);
	can	= new TCanvas("BestFit", "BestFit", 1500, 800);
    can->Divide(4,3);
	Double_t*	BestFitWidth[7];
	Double_t*	BestCounts[7];
    can->cd(1);
    {
        BestFitWidth[0]	= new Double_t[4];
        BestCounts[0]	= new Double_t[2];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit1/Final/Final_IM");
        BestCounts[0][0]	= help1->GetEntries();
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit1/Final/Final_IM");
        BestCounts[0][1]	= help2->GetEntries();
        BestFitFit(help1, help2, BestFitWidth[0]);
        help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
    }
    can->cd(2);
	{
        BestFitWidth[1]	= new Double_t[4];
        BestCounts[1]	= new Double_t[2];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit3/Final/Final_IM");
        BestCounts[1][0]	= help1->GetEntries();
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit3/Final/Final_IM");
        BestCounts[1][1]	= help2->GetEntries();
        BestFitFit(help1, help2, BestFitWidth[1]);
        help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
	}
    can->cd(3);
	{
        BestFitWidth[2]	= new Double_t[4];
        BestCounts[2]	= new Double_t[2];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/Final/Final_IM");
        BestCounts[2][0]	= help1->GetEntries();
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit4/Final/Final_IM");
        BestCounts[2][1]	= help2->GetEntries();
        BestFitFit(help1, help2, BestFitWidth[2]);
        help1->SetLineColor(kMagenta);
        help1->Draw();
        //help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
	}
    can->cd(4);
	{
        BestFitWidth[3]	= new Double_t[4];
        BestCounts[3]	= new Double_t[2];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit7Proton/Final/Final_IM");
        BestCounts[3][0]	= help1->GetEntries();
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit7Proton/Final/Final_IM");
        BestCounts[3][1]	= help2->GetEntries();
        BestFitFit(help1, help2, BestFitWidth[3]);
        help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
	}
	
    can->cd(5);
    {
        BestFitWidth[4]	= new Double_t[4];
        BestCounts[4]	= new Double_t[2];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithoutProton/MM_Cut/fit1/Final/Final_IM");
        BestCounts[4][0]	= help1->GetEntries();
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithoutProton/MM_Cut/fit1/Final/Final_IM");
        BestCounts[4][1]	= help2->GetEntries();
        BestFitFit(help1, help2, BestFitWidth[4]);
        help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
    }
    can->cd(6);
	{
        BestFitWidth[5]	= new Double_t[4];
        BestCounts[5]	= new Double_t[2];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithoutProton/MM_Cut/fit3/Final/Final_IM");
        BestCounts[5][0]	= help1->GetEntries();
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithoutProton/MM_Cut/fit3/Final/Final_IM");
        BestCounts[5][1]	= help2->GetEntries();
        BestFitFit(help1, help2, BestFitWidth[5]);
        help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
	}
    can->cd(7);
	{
        BestFitWidth[6]	= new Double_t[4];
        BestCounts[6]	= new Double_t[2];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithoutProton/MM_Cut/fit4/Final/Final_IM");
        BestCounts[6][0]	= help1->GetEntries();
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithoutProton/MM_Cut/fit4/Final/Final_IM");
        BestCounts[6][1]	= help2->GetEntries();
        BestFitFit(help1, help2, BestFitWidth[6]);
        help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");
	}
	Double_t	x[4];
	Double_t	dx[4];
	Double_t	y[4];
	Double_t	dy[4];
    can->cd(9);
	{
        for(int i=0; i<4; i++)
		{
			x[i]	= i+1;
			dx[i]	= 0;
			y[i]	= BestFitWidth[i][2] / BestFitWidth[i][0];
			dy[i]	= TMath::Sqrt(((BestFitWidth[i][3] / BestFitWidth[i][0])*(BestFitWidth[i][3] / BestFitWidth[i][0])) + ((BestFitWidth[i][2]*BestFitWidth[i][1] / BestFitWidth[i][0])*(BestFitWidth[i][2]*BestFitWidth[i][1] / BestFitWidth[i][0])));
		}
        TGraphErrors*	bestFit = new TGraphErrors(4, x, y, dx, dy);
		bestFit->Draw();
	}
	can->cd(10);
	{
        for(int i=0; i<4; i++)
		{
			x[i]	= i+1;
			dx[i]	= 0;
			y[i]	= BestCounts[i][1] / BestCounts[i][0];
			dy[i]	= TMath::Sqrt(((BestFitWidth[i][3] / BestFitWidth[i][0])*(BestFitWidth[i][3] / BestFitWidth[i][0])) + ((BestFitWidth[i][2]*BestFitWidth[i][1] / BestFitWidth[i][0])*(BestFitWidth[i][2]*BestFitWidth[i][1] / BestFitWidth[i][0])));
		}
        TGraphErrors*	bestFit = new TGraphErrors(4, x, y, dx, dy);
		bestFit->Draw();
	}
	
	out->cd();
	can->Write();
}

void	BestFit(const char* mcSignalFileName, const char* mcBGFileName)
{
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
	
	BestFit(mcSignalFile, mcBGFile, out);
}
