

#include "macros/ShowPhysics.C"


void	BestFitFit(const TH1D* hist1, const TH1D* hist2, Double_t* width, Double_t* count)
{
	//TH1D*	h = (TH1D*)hist1->Clone();
	TH1D*	h1 = (TH1D*)hist1->Clone();
	TH1D*	h2 = (TH1D*)hist2->Clone();
	//h->Add(hist2);
	
	TF1*	fit1 = new TF1("fitfkt1", "gaus", 750, 1100);
	fit1->SetParameters(h1->GetMaximum(), 957, 10);
	fit1->SetParLimits(0, h1->GetMaximum()/10, h1->GetMaximum()*2);
	fit1->SetParLimits(1, 940, 990);
	fit1->SetParLimits(2, 5, 25);
    fit1->SetLineColor(kBlack);
	TF1*	fit2 = new TF1("fitfkt2", "gaus", 750, 1100);
	fit2->SetParameters(h2->GetMaximum(), 957, 50);
	fit2->SetParLimits(0, 0, h2->GetMaximum()*2);
	fit2->SetParLimits(1, 700, 1200);
	fit2->SetParLimits(2, 25, 100);
    fit2->SetLineColor(kBlue);

	h1->Fit(fit1, "R0");
    h1->SetAxisRange(750, 1100);
    h1->SetLineColor(kMagenta);
	h1->Draw();
	fit1->Draw("SAME");
	h2->Fit(fit2, "R0");
    h2->SetAxisRange(750, 1100);
    h2->SetLineColor(kGreen);
	h2->Draw("SAME");
	fit2->Draw("SAME");
    //h->SetAxisRange(750, 1100);
    //h->SetLineColor(kBlack);
	//h->Draw("SAME");
	
	width[0] = fit1->GetParameter(2);
	width[1] = fit1->GetParError(2);
	width[2] = fit2->GetParameter(2);
	width[3] = fit2->GetParError(2);
	
	count[0] = fit1->GetParameter(0)*fit1->GetParameter(2)*TMath::Sqrt(TMath::Pi());
	count[1] = TMath::Sqrt( (fit1->GetParError(0)*fit1->GetParError(0)*fit1->GetParameter(2)*fit1->GetParameter(2)*TMath::Pi())+
							(fit1->GetParameter(0)*fit1->GetParameter(0)*fit1->GetParError(2)*fit1->GetParError(2)*TMath::Pi()));
	count[2] = fit2->GetParameter(0)*fit2->GetParameter(2)*TMath::Sqrt(TMath::Pi());
	count[3] = TMath::Sqrt( (fit2->GetParError(0)*fit2->GetParError(0)*fit2->GetParameter(2)*fit2->GetParameter(2)*TMath::Pi())+
							(fit2->GetParameter(0)*fit2->GetParameter(0)*fit2->GetParError(2)*fit2->GetParError(2)*TMath::Pi()));
}

void	BestFit(const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out)
{	
        /*TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fit4/IM");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fit4/IM");
        help1->SetLineColor(kMagenta); 
		help1->SetAxisRange(850, 1050);       
		help1->SetStats(0);
		help1->GetXaxis()->SetTitle("inv. Mass (#eta#pi^{0}#pi^{0}) [MeV/c]");
		help1->SetTitle("");
        help1->Draw();    
		help2->SetStats(0);
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
        
	can	= new TCanvas("BestFit", "BestFit", 1500, 800);
    can->Divide(4,3);
	Double_t*	BestFitWidth[8];
	Double_t*	BestCounts[8];
    can->cd(1);
    {
        BestFitWidth[0]	= new Double_t[4];
        BestCounts[0]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fit4/IM");
        BestCounts[0][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fit4/IM");
        BestCounts[0][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[0], BestCounts[0]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
    }
    can->cd(2);
	{
        BestFitWidth[1]	= new Double_t[4];
        BestCounts[1]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fitBeam4/IM");
        BestCounts[1][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fitBeam4/IM");
        BestCounts[1][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[1], BestCounts[1]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
	}
    can->cd(3);
	{
        BestFitWidth[2]	= new Double_t[4];
        BestCounts[2]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fitProton6/IM");
        BestCounts[2][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fitProton6/IM");
        BestCounts[2][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[2], BestCounts[2]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw();
        //help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
	}
    can->cd(4);
	{
        BestFitWidth[3]	= new Double_t[4];
        BestCounts[3]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fitBeamProton6/IM");
        BestCounts[3][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fitBeamProton6/IM");
        BestCounts[3][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[3], BestCounts[3]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
	}
	
    can->cd(5);
    {
        BestFitWidth[4]	= new Double_t[4];
        BestCounts[4]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fit4/Vertex/IM");
        BestCounts[4][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fit4/Vertex/IM");
        BestCounts[4][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[4], BestCounts[4]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
    }
    can->cd(6);
	{
        BestFitWidth[5]	= new Double_t[4];
        BestCounts[5]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fitBeam4/Vertex/IM");
        BestCounts[5][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fitBeam4/Vertex/IM");
        BestCounts[5][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[5], BestCounts[5]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
	}
    can->cd(7);
	{
        BestFitWidth[6]	= new Double_t[4];
        BestCounts[6]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fitProton6/Vertex/IM");
        BestCounts[6][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fitProton6/Vertex/IM");
        BestCounts[6][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[6], BestCounts[6]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
	}
    can->cd(8);
	{
        BestFitWidth[7]	= new Double_t[4];
        BestCounts[7]	= new Double_t[4];
        TH1D*	help1		= (TH1D*)mcSignalFile->Get("WithProton/fitBeamProton6/Vertex/IM");
        BestCounts[7][0]	= help1->Integral("width");
        TH1D*	help2		= (TH1D*)mcBGFile->Get("WithProton/fitBeamProton6/Vertex/IM");
        BestCounts[7][1]	= help2->Integral("width");
        BestFitFit(help1, help2, BestFitWidth[7], BestCounts[7]);
        /*help1->SetLineColor(kMagenta);
        help1->Draw("SAME");
        help2->SetLineColor(kGreen);
        help2->Draw("SAME");*/
	}
	Double_t	x[8];
	Double_t	dx[8];
	Double_t	y[8];
	Double_t	dy[8];
	can->cd(9);
	{
		for(int i=0; i<8; i++)
		{
			x[i]	= i+1;
			dx[i]	= 0;
			y[i]	= BestCounts[i][0];
			dy[i]	= BestCounts[i][1];
		}
        TGraphErrors*	bestFit = new TGraphErrors(8, x, y, dx, dy);
		bestFit->Draw();
	}
	can->cd(10);
	{
		for(int i=0; i<8; i++)
		{
			x[i]	= i+1;
			dx[i]	= 0;
			y[i]	= BestCounts[i][2];
			dy[i]	= BestCounts[i][3];
		}
        TGraphErrors*	bestFit = new TGraphErrors(8, x, y, dx, dy);
		bestFit->Draw();
	}
    can->cd(11);
	{
        for(int i=0; i<8; i++)
		{
			x[i]	= i+1;
			dx[i]	= 0;
			if(BestFitWidth[i][2]==0)
			{
				y[i]	= 0;
				dy[i]	= 1;
			}
			else
			{
				y[i]	= BestFitWidth[i][0] / BestFitWidth[i][2];
				dy[i]	= TMath::Sqrt(((BestFitWidth[i][1] / BestFitWidth[i][2])*(BestFitWidth[i][1] / BestFitWidth[i][2])) + ((BestFitWidth[i][0]*BestFitWidth[i][3] / (BestFitWidth[i][2]*BestFitWidth[i][2]))*(BestFitWidth[i][0]*BestFitWidth[i][3] / (BestFitWidth[i][2]*BestFitWidth[i][2]))));
			}
		}
        TGraphErrors*	bestFit = new TGraphErrors(8, x, y, dx, dy);
		bestFit->Draw();
	}
	can->cd(12);
	{
        for(int i=0; i<8; i++)
		{
			x[i]	= i+1;
			dx[i]	= 0;
			if(BestCounts[i][1]==0)
			{
				y[i]	= 0;
				dy[i]	= 1;
			}
			else
			{
				y[i]	= BestCounts[i][0] / BestCounts[i][2];
				dy[i]	= TMath::Sqrt((BestCounts[i][1]*BestCounts[i][1]/ (BestCounts[i][2]*BestCounts[i][2])) + (BestCounts[i][0]*BestCounts[i][0]*BestCounts[i][3]*BestCounts[i][3] / (BestCounts[i][2]*BestCounts[i][2]*BestCounts[i][2]*BestCounts[i][2])));
				//std::cout << y[i] << "     " << dy[i] << std::endl;
			}
		}
        TGraphErrors*	bestFit = new TGraphErrors(8, x, y, dx, dy);
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
