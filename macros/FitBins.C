int	TaggE[48] = {
1577.31,
1573.83,
1570.43,
1567.06,
1563.71,
1560.39,
1557.11,
1553.8,
1550.53,
1547.24,
1543.98,
1540.7,
1537.44,
1534.17,
1530.91,
1527.66,
1524.4,
1521.13,
1517.89,
1514.63,
1511.36,
1508.11,
1504.86,
1501.62,
1498.39,
1495.16,
1491.94,
1488.71,
1485.49,
1482.26,
1479.03,
1475.8,
1472.58,
1469.36,
1466.13,
1462.89,
1459.67,
1456.4,
1453.17,
1449.92,
1446.66,
1443.42,
1440.17,
1436.94,
1433.7,
1430.49,
1427.32,
1424.42};

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
		array[i] += 0.02;
	for(int i=0; i<48; i++)
		darray[i]	= 2.0/100.0;
}
void ReconstructStart(Double_t* array)
{
	for(int i=0; i<38; i++)
		array[i]	= (84910.0-550.0)/(2*38);
	array[38]	= 550.0/2.0;
	for(int i=39; i<48; i++)
		array[i]	= 0.0;
}
void ReconstructEff(const TH2D* h, Double_t* array, Double_t* darray)
{
	Double_t	help[48];
	ReconstructStart(help);
	TH1D*	slice = 0;
	for(int i=0; i<48; i++)
	{
		slice	= (TH1D*)h->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
		if(help[i]==0)
		{
			array[i]	= 0.015;
			darray[i]	= 0.015;
		}
		else
		{
			array[i]	= slice->GetEntries()/help[i];
			Double_t	buf;
			buf			= TMath::Sqrt(slice->GetEntries())/help[i];
			darray[i]	= slice->GetEntries()*TMath::Sqrt(help[i])/(help[i]*help[i]);
			darray[i]	*= darray[i];
			darray[i]	+= buf*buf;
			darray[i]	= TMath::Sqrt(darray[i]);
		}
	}
}
void	FitBinsFitAll(const TH1D* hist, Double_t* par)
{
	TH1D* h	= hist->Clone();
	TF1*	fit = new TF1("fitfkt", "gaus(0)+gaus(3)", 850, 1050);
	fit->SetParameters(h->GetMaximum(), 957, 10, h->GetMaximum()/20, 950, 50);
	fit->SetParLimits(0, h->GetMaximum()/10, h->GetMaximum());
	fit->SetParLimits(1, 940, 990);
	fit->SetParLimits(2, 5, 25);
	fit->SetParLimits(3, 0, h->GetMaximum()/2);
	fit->SetParLimits(4, 700, 1200);
	fit->SetParLimits(5, 25, 100);

	h->Fit(fit, "R");
    h->SetAxisRange(850, 1050);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("inv. Mass (#eta#pi^{0}#pi^{0}) [MeV/c]");
    h->SetTitle("");
	h->Draw();
	
	par[0] = fit->GetParameter(0);
	par[1] = fit->GetParameter(1);
	par[2] = fit->GetParameter(2);
	par[3] = fit->GetParameter(3);
	par[4] = fit->GetParameter(4);
	par[5] = fit->GetParameter(5);
	
	TF1*	fit2 = new TF1("fitfktfg", "gaus(0)", 850, 1050);
	fit2->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));
    fit2->SetLineColor(kMagenta);
	fit2->Draw("SAMES");
	TF1*	fit3 = new TF1("fitfktbg", "gaus(0)", 850, 1050);
	fit3->SetParameters(fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5));
    fit3->SetLineColor(kGreen);
	fit3->Draw("SAMES");
}

void	FitBinsFit(const TH1D* hist, const Double_t TaggEMin, const Double_t TaggEMax, const Double_t* par, Double_t* parResult, Double_t* parResultError)
{
	TH1D* h = hist->Clone();
	TF1*	fit = new TF1("fitfkt", "gaus(0)+gaus(3)", 850, 1050);
	fit->SetParameters(h->GetMaximum(), par[1], par[2], h->GetMaximum()/20, par[4], par[5]);
	fit->SetParLimits(0, 0, h->GetMaximum());
	fit->SetParLimits(1, par[1]-10, par[1]+10);
	fit->SetParLimits(2, par[2]/2, par[2]*2);
	fit->SetParLimits(3, 0, h->GetMaximum()/2);
	fit->SetParLimits(4, par[4]-25, par[4]+50);
	fit->SetParLimits(5, par[5], par[5]*2);

	h->Fit(fit, "R");
	h->SetStats(0);
	h->SetTitle(TString("E_{#gamma} (").Append(TString().Itoa(TaggEMin, 10)).Append(" MeV - ").Append(TString().Itoa(TaggEMax, 10)).Append(" MeV) ").Data());
    h->SetAxisRange(850, 1050);
	h->Draw();
	
	parResult[0] = fit->GetParameter(0);
	parResult[1] = fit->GetParameter(1);
	parResult[2] = fit->GetParameter(2);
	parResult[3] = fit->GetParameter(3);
	parResult[4] = fit->GetParameter(4);
	parResult[5] = fit->GetParameter(5);
	
	parResultError[0] = fit->GetParError(0);
	parResultError[1] = fit->GetParError(1);
	parResultError[2] = fit->GetParError(2);
	parResultError[3] = fit->GetParError(3);
	parResultError[4] = fit->GetParError(4);
	parResultError[5] = fit->GetParError(5);
	
	TF1*	fit2 = new TF1("fitfktfg", "gaus(0)", 850, 1050);
	fit2->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));
    fit2->SetLineColor(kMagenta);
	fit2->Draw("SAMES");
	TF1*	fit3 = new TF1("fitfktbg", "gaus(0)", 850, 1050);
	fit3->SetParameters(fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5));
    fit3->SetLineColor(kGreen);
	fit3->Draw("SAMES");
}

void	FitBins(const TFile* dataFile, const TFile* scalerFile, const TFile* out, const Int_t addedChannels, const bool prompt = false)
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
	Double_t**	parResult 		= new Double_t*[ch.size()];
	Double_t**	parResultError 	= new Double_t*[ch.size()];
	TH2D*		data = 0;
	if(prompt)
	{
		std::cout << "using prompt hist." << std::endl;
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/Final/TaggerBinning/BackgroundSubstraction/Final_IM_Bins_Prompt");
	}
	else
	{
		std::cout << "using prompt random subtracted hist." << std::endl;
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/Final/TaggerBinning/Final_IM_Bins");
	}
		
	maincan	= new TCanvas("FitBinsMain", "FitBinsMain", 1500, 800);
	maincan->Divide(3, 4);
	maincan->cd(1);
	FitBinsFitAll((TH1D*)data->ProjectionX(TString("BinUpTo35").Data(), 0, -1), par);
	scan	= new TCanvas("FitBinsSCan", "FitBinsSCan", 1500, 800);
	scan->Divide(1, 1);
	scan->cd(1);
	FitBinsFitAll((TH1D*)data->ProjectionX(TString("BinUpTo35").Data(), 0, -1), par);
	
	can	= new TCanvas("FitBins", "FitBins", 1500, 800);
	can->Divide(4,TMath::Ceil(double(ch.size())/4));
	
	
	Double_t	x9[48];
	Double_t	dx9[48];
	for(int i=0; i<ch.size(); i++)
	{
		x9[i]	= 0;
		for(int k=ch[i]; k<(ch[i]+size[i]); k++)
			x9[i]	+= TaggE[k];
		x9[i]	/= size[i];
		x9[i]	/= 1000;
		dx9[i]	= 2.5 * addedChannels/1000;
	}
	
	for(int i=0; i<ch.size(); i++)
	{
		can->cd(i+1);
		parResult[i] 		= new Double_t[6];
		parResultError[i] 	= new Double_t[6];
		FitBinsFit((TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1), TaggE[ch[i]], TaggE[ch[i]+size[i]], par, parResult[i], parResultError[i]);
	}
	
	Double_t	x1[48];
	Double_t	dx1[48];
	Double_t	y1[48];
	Double_t	dy1[48];
	for(int i=0; i<ch.size(); i++)
	{
		x1[i]	= i+1;
		dx1[i]	= 0;
		y1[i]	= parResult[i][0] * parResult[i][2] * TMath::Sqrt(TMath::Pi());
		dy1[i]	= TMath::Sqrt(((parResultError[i][0] * parResult[i][2])*(parResultError[i][0] * parResult[i][2])) + ((parResult[i][0] * parResultError[i][2])*(parResult[i][0] * parResultError[i][2])));
	}
	maincan->cd(2);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x1, y1, dx1, dy1);
	graph->Draw();
	
	
	Double_t	x2[48];
	Double_t	dx2[48];
	Double_t	y2[48];
	Double_t	dy2[48];
	Double_t	taggEff[48];
	Double_t	dTaggEff[48];
	TaggEff(taggEff, dTaggEff);
	for(int i=0; i<ch.size(); i++)
	{
		x2[i]	= i+1;
		dx2[i]	= 0;
		y2[i]	= 0;
		for(int k=ch[i]; k<(ch[i]+size[i]); k++)
		{
			y2[i]	+= taggEff[k];
		}
		y2[i]	/= size[i];
		dy2[i]	= 2.0/100.0;
	}
	maincan->cd(3);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x2, y2, dx2, dy2);
	graph->Draw();
	
	
	Double_t	x3[48];
	Double_t	dx3[48];
	Double_t	y3[48];
	Double_t	dy3[48];
	for(int i=0; i<ch.size(); i++)
	{
		x3[i]	= i+1;
		dx3[i]	= 0;
		if(y2[i]==0)
		{
			y3[i]	= 0;
			dy3[i]	= 0;
		}
		else
		{
			y3[i]		= y1[i]/y2[i];
			Double_t	buf;
			buf			= dy1[i]/y2[i];
			dy3[i]		= dy2[i]*y1[i]/(y2[i]*y2[i]);
			dy3[i]		*= dy3[i];
			dy3[i]		+= buf*buf;
			dy3[i]	 	 = TMath::Sqrt(dy3[i]);
		}
	}
	maincan->cd(4);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x3, y3, dx3, dy3);
	graph->Draw();
	
	
	Double_t	x4[48];
	Double_t	dx4[48];
	Double_t	y4[48];
	Double_t	dy4[48];
	Double_t	recStart[48];
	ReconstructStart(recStart);	
	for(int i=0; i<ch.size(); i++)
	{
		x4[i]	= i+1;
		dx4[i]	= 0;
		y4[i]	= 0;
		for(int k=ch[i]; k<(ch[i]+size[i]); k++)
		{
			y4[i]	+= recStart[k];
		}
		//y4[i]	/= size[i];
		dy4[i]	= TMath::Sqrt(y4[i]);
	}
	maincan->cd(5);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x4, y4, dx4, dy4);
	graph->Draw();
	
	Double_t	x5[48];
	Double_t	dx5[48];
	Double_t	y5[48];
	Double_t	dy5[48];
	Double_t	recEff[48];
	Double_t	dRecEff[48];
	ReconstructEff(data, recEff, dRecEff);
	for(int i=0; i<ch.size(); i++)
	{
		x5[i]	= i+1;
		dx5[i]	= 0;
		y5[i]	= 0;
		dy5[i]	= 0;
		for(int k=ch[i]; k<(ch[i]+size[i]); k++)
		{
			y5[i]	+= recEff[k];
			dy5[i]	+= dRecEff[k]*dRecEff[k];
		}
		y5[i]	/= size[i];
		dy5[i]	 = TMath::Sqrt(dy5[i]);
	}
	maincan->cd(6);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x5, y5, dx5, dy5);
	graph->Draw();
	
	
	
	
	Double_t	x6[48];
	Double_t	dx6[48];
	Double_t	y6[48];
	Double_t	dy6[48];
	for(int i=0; i<ch.size(); i++)
	{
		x6[i]	= i+1;
		dx6[i]	= 0;
		if(y5[i]==0)
		{
			y6[i]	= 0.01;
			dy6[i]	= 0.01;
		}
		else
		{
			y6[i]		= y3[i]/y5[i];
			Double_t	buf;
			buf			= dy3[i]/y5[i];
			dy6[i]		= dy5[i]*y3[i]/(y5[i]*y5[i]);
			dy6[i]		*= dy6[i];
			dy6[i]		+= buf*buf;
			dy6[i]	 	 = TMath::Sqrt(dy6[i]);
		}
	}
	maincan->cd(7);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x6, y6, dx6, dy6);
	graph->Draw();
	
	out->cd();
	maincan->Write();
	can->Write();
	
	
	TH1D*	scaler		= (TH1D*)scalerFile->Get("EPT_ScalerCorT");
	if(!scaler)
	{
		std::cout << "Can not open scaler hist " << std::endl;
		return;
	}
	Double_t	x7[48];
	Double_t	dx7[48];
	Double_t	sc[48];
	Double_t	dsc[48];
	for(int i=0; i<ch.size(); i++)
	{
		x7[i]	= i+1;
		dx7[i]	= 0;
		sc[i] 	= 0;
		for(int k=ch[i]; k<(ch[i]+size[i]); k++)
			sc[i] 	+= scaler->GetBinContent(k+1);
		dsc[i]	= TMath::Sqrt(sc[i]);
	}
	maincan->cd(8);
    scaler->Draw();
	maincan->cd(9);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x7, sc, dx7, dsc);
	graph->Draw();
	
	Double_t	x8[48];
	Double_t	dx8[48];
	Double_t	y8[48];
	Double_t	dy8[48];
	for(int i=0; i<ch.size(); i++)
	{
		x8[i]	= i+1;
		dx8[i]	= 0;
		y8[i]	= (4/0.08491)*y6[i]/sc[i];
		Double_t	buf;
		buf			= dy6[i]/sc[i];
		dy8[i]		= dsc[i]*y6[i]/(sc[i]*sc[i]);
		dy8[i]		*= dy8[i];
		dy8[i]		+= buf*buf;
		dy8[i]	 	 = (4/0.08491)*TMath::Sqrt(dy8[i]);
	}
	maincan->cd(10);
    TGraphErrors*	graph = new TGraphErrors(ch.size(), x8, y8, dx8, dy8);
	graph->Draw();
	
	Double_t	y9[48];
	Double_t	dy9[48];
	for(int i=0; i<ch.size(); i++)
	{
		y9[i]	= y8[i]*1000000;
		dy9[i]	= dy8[i]*1000000;
	}
	endResultCan	= new TCanvas("FitBinsEndResult", "FitBinsEndResult", 1500, 800);
	endResultCan->Divide(1, 1);
	endResultCan->cd(1);
    TGraphErrors*	graph = new TGraphErrors(ch.size()-3, x9, y9, dx9, dy9);
    graph->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
    graph->GetYaxis()->SetTitle("[#mub]");
    graph->SetTitle("");
	graph->Draw();
	
	out->cd();
	maincan->Write();
	can->Write();
}

void	FitBins(const char* dataFileName, const char* scalerFileName, const Int_t addedChannels)
{
	TFile*	dataFile		= TFile::Open(dataFileName);
	if(!dataFile)
	{
		std::cout << "Can not open dataFile " << dataFileName << std::endl;
		return;
	}
	TFile*	mcSignalFile		= TFile::Open(mcSignalFileName);
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << mcSignalFileName << std::endl;
		return;
	}
	TFile*	scalerFile			= TFile::Open(scalerFileName);
	if(!scalerFile)
	{
		std::cout << "Can not open scalerFile " << scalerFileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("result.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	FitBins(dataFile, mcSignalFile, scalerFile, out, addedChannels);
}
