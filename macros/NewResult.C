

#include "macros/ShowPhysics.C"
#include "macros/NewBestFit.C"
#include "macros/NewFitBins.C"



void	Fit(const TH1D* data, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& etapWidth, Double_t& BGWidth)
{
	nEtapCut	= 0;
	for(int i=1; i<=data->GetNbinsX(); i++)
	{
		if(data->GetBinCenter(i)>(957-8.83) && data->GetBinCenter(i)<(957+8.83))
			nEtapCut	+= data->GetBinContent(i);
	}
	
	Double_t	max	= data->GetMaximum();
	
	TF1*	fit = new TF1("fit","gaus(0) + gaus(3)",850,1020);
	fit->SetParameters(max, 958, etapWidth, max/10, 930, BGWidth);
	fit->SetParLimits(0, 0, max *2);
	fit->SetParLimits(1, 950, 970);
	fit->SetParLimits(2, etapWidth, etapWidth);
	fit->SetParLimits(3, 0, max/2);
	fit->SetParLimits(4, 900, 980);
	fit->SetParLimits(5, BGWidth, BGWidth);
	
	data->Fit(fit, "R");
	
	nEtapFit	= fit->GetParameter(0) * fit->GetParameter(2) * 1.7724;
	dNEtapFit	= fit->GetParError(0) * fit->GetParameter(2) * 1.7724;
}
void	FitProfile(const TH2D* data, const int i, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& etapWidth, Double_t& BGWidth)
{
	TH1D*	slice		= (TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
	slice->SetAxisRange(850, 1020);
	slice->Draw();
	
	Fit(slice, nEtapFit, dNEtapFit, nEtapCut, etapWidth, BGWidth);
}

void	FitSim(const TH1D* data, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& width)
{
	nEtapCut	= 0;
	for(int i=1; i<=data->GetNbinsX(); i++)
	{
		if(data->GetBinCenter(i)>(957-8.83) && data->GetBinCenter(i)<(957+8.83))
			nEtapCut	+= data->GetBinContent(i);
	}
	
	Double_t	max	= data->GetMaximum();
	
	TF1*	fit = new TF1("fit","gaus(0)",850,1020);
	fit->SetParameters(max, 958, 10, max/10, 930, 35);
	fit->SetParLimits(0, max/10, max *2);
	fit->SetParLimits(1, 900, 1000);
	fit->SetParLimits(2, 5, 75);
	
	data->Fit(fit, "R");
	
	nEtapFit	= fit->GetParameter(0) * fit->GetParameter(2) * 1.7724;
	dNEtapFit	= fit->GetParError(0) * fit->GetParameter(2) * 1.7724;
	
	width	= fit->GetParameter(2);
}
void	FitProfileSim(const TH2D* data, const int i, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& width)
{
	TH1D*	slice		= (TH1D*)(data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1))->Clone();
	slice->SetAxisRange(850, 1020);
	slice->Draw();
		
	FitSim(slice, nEtapFit, dNEtapFit, nEtapCut, width);
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
	
	
	ShowPhysics(dataFile, mcSignalFile, mcBGFile, out);
	BestFit(mcSignalFile, mcBGFile, out);
	
	/*can	= new TCanvas("CanFit", "Fit", 1500, 800);
	can->Divide(3, 2);
	TH2D*	data;
	{
		can->cd(1)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/fit4/IM");
		data->Draw("COL");
		can->cd(2)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/fit4/_MM");
		data->Draw("COL");
		can->cd(3);
		OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/Final/Final_IM");
		can->cd(4)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/_Sub0IM");
		data->Draw("COL");
		can->cd(5)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/_Sub1IM");
		data->Draw("COL");
		can->cd(6)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/_Sub2IM");
		data->Draw("COL");
		//can->cd();
		//OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/_fit4_ChiSq");
		//can->cd(4);
		//OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/_fit4_ConfLev");
	}
	out->cd();
	can->Write();*/
	
	//TFile*	scalerFile		= TFile::Open("/home/ott/ScalerPhysics_CB.root");
	//if(!scalerFile)
	//{
	//	std::cout << "Can not open scalerFile /home/ott/ScalerPhysics_CB.root" << std::endl;
	//	return;
	//}
	FitBins(dataFile, dataFile, out, 3, "WithProton/fitBeam4/TaggerBinning/IM_Bins");
}
