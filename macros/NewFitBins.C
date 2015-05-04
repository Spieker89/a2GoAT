#include "macros/NewBestFit.C"


double	TaggE[48] = {
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

struct	Value
{
	Double_t	value;
	Double_t	error;
};
struct	FitValuesSignal
{
	Value	factor;
	Value	mean;
	Value	sigma;
	Double_t	count;
};
struct	FitValuesData
{
	FitValuesSignal	signal;
	FitValuesSignal	bg;
};
struct	FitValuesBG
{
	FitValuesSignal	gaus;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
	}pol2;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
		Value	a3;
	}pol3;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
		Value	a3;
		Value	a4;
	}pol4;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
		Value	a3;
		Value	a4;
		Value	a5;
	}pol5;
	Double_t	count;
};
struct	FitValues
{
	FitValuesSignal	signal;
	FitValuesBG		bg;
};

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

void	FitBinsSimGaus(const TH1D* data, int channel, int size, FitValues& fitValues, const bool signal)
{
	TF1*	fit = new TF1("fitfktSignal", "gaus(0)", 850, 1050);
	if(signal)
	{
		data->SetLineColor(kRed);
		fit->SetLineColor(kMagenta);
		fit->SetParameters(data->GetMaximum(), 957, 10);
		fit->SetParLimits(1, 950, 965);
		fit->SetParLimits(2, 5, 20);
	}
	else
	{
		data->SetLineColor(kBlue);
		fit->SetLineColor(kBlue);
		fit->SetParameters(data->GetMaximum(), 930, 50);
		fit->SetParLimits(1, 850, 1050);
		fit->SetParLimits(2, 25, 250);
	}
	fit->SetParLimits(0, 0, data->GetMaximum()*1.5);
	
	TString	title("E_{#gamma} (");
	if(channel==0)
	{
		title.Append(TString().Itoa(int(TaggE[size]-((TaggE[size]-TaggE[size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[0]+((TaggE[0]-TaggE[1])/2)), 10));
		title.Append(" MeV) ");
	}
	else if(channel+size>=47)
	{
		title.Append(TString().Itoa(int(TaggE[47]-((TaggE[46]-TaggE[47])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	else
	{
		title.Append(TString().Itoa(int(TaggE[channel+size]-((TaggE[channel+size]-TaggE[channel+size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	data->SetTitle(title.Data());
	data->GetXaxis()->SetTitle("IM #gamma#gamma [MeV]");	
	data->Fit(fit, "R0");
	data->SetStats(0);		
    data->SetAxisRange(850, 1050);
	data->Draw("SAME");
	fit->Draw("SAME");
	
	if(signal)
	{
		fitValues.signal.factor.value	= fit->GetParameter(0);
		fitValues.signal.factor.error	= fit->GetParError(0);
		fitValues.signal.mean.value		= fit->GetParameter(1);
		fitValues.signal.mean.error		= fit->GetParError(1);
		fitValues.signal.sigma.value	= fit->GetParameter(2);
		fitValues.signal.sigma.error	= fit->GetParError(2);
		fitValues.signal.count			= data->Integral();
	}
	else
	{
		fitValues.bg.gaus.factor.value	= fit->GetParameter(0);
		fitValues.bg.gaus.factor.error	= fit->GetParError(0);
		fitValues.bg.gaus.mean.value	= fit->GetParameter(1);
		fitValues.bg.gaus.mean.error	= fit->GetParError(1);
		fitValues.bg.gaus.sigma.value	= fit->GetParameter(2);
		fitValues.bg.gaus.sigma.error	= fit->GetParError(2);
		fitValues.bg.count				= data->Integral();
	}
}
void	FitBinsSimPol(const TH1D* data, int channel, int size, FitValues& fitValues)
{
	TF1*	fit2 = new TF1("fitfktBG2", "pol2(0)", 850, 1050);
	TF1*	fit3 = new TF1("fitfktBG3", "pol3(0)", 850, 1050);
	TF1*	fit4 = new TF1("fitfktBG4", "pol4(0)", 850, 1050);
	TF1*	fit5 = new TF1("fitfktBG5", "pol5(0)", 850, 1050);
	fit2->SetLineColor(kGreen);
	fit3->SetLineColor(kOrange);
	fit4->SetLineColor(kYellow);
	fit5->SetLineColor(kCyan);
	fit2->SetParameters(0, 0, 0);
	fit3->SetParameters(0, 0, 0, 0);
	fit4->SetParameters(0, 0, 0, 0, 0);
	fit5->SetParameters(0, 0, 0, 0, 0, 0);
	
	data->Fit(fit2, "R0");
	data->Fit(fit3, "R0");
	data->Fit(fit4, "R0");
	data->Fit(fit5, "R0");
	fit2->Draw("SAME");
	fit3->Draw("SAME");
	fit4->Draw("SAME");
	fit5->Draw("SAME");
	
	fitValues.bg.pol2.a0.value	= fit2->GetParameter(0);
	fitValues.bg.pol2.a0.error	= fit2->GetParError(0);
	fitValues.bg.pol2.a1.value	= fit2->GetParameter(1);
	fitValues.bg.pol2.a1.error	= fit2->GetParError(1);
	fitValues.bg.pol2.a2.value	= fit2->GetParameter(2);
	fitValues.bg.pol2.a2.error	= fit2->GetParError(2);
	
	fitValues.bg.pol3.a0.value	= fit3->GetParameter(0);
	fitValues.bg.pol3.a0.error	= fit3->GetParError(0);
	fitValues.bg.pol3.a1.value	= fit3->GetParameter(1);
	fitValues.bg.pol3.a1.error	= fit3->GetParError(1);
	fitValues.bg.pol3.a2.value	= fit3->GetParameter(2);
	fitValues.bg.pol3.a2.error	= fit3->GetParError(2);
	fitValues.bg.pol3.a3.value	= fit3->GetParameter(3);
	fitValues.bg.pol3.a3.error	= fit3->GetParError(3);
	
	fitValues.bg.pol4.a0.value	= fit4->GetParameter(0);
	fitValues.bg.pol4.a0.error	= fit4->GetParError(0);
	fitValues.bg.pol4.a1.value	= fit4->GetParameter(1);
	fitValues.bg.pol4.a1.error	= fit4->GetParError(1);
	fitValues.bg.pol4.a2.value	= fit4->GetParameter(2);
	fitValues.bg.pol4.a2.error	= fit4->GetParError(2);
	fitValues.bg.pol4.a3.value	= fit4->GetParameter(3);
	fitValues.bg.pol4.a3.error	= fit4->GetParError(3);
	fitValues.bg.pol4.a4.value	= fit4->GetParameter(4);
	fitValues.bg.pol4.a4.error	= fit4->GetParError(4);
	
	fitValues.bg.pol5.a0.value	= fit5->GetParameter(0);
	fitValues.bg.pol5.a0.error	= fit5->GetParError(0);
	fitValues.bg.pol5.a1.value	= fit5->GetParameter(1);
	fitValues.bg.pol5.a1.error	= fit5->GetParError(1);
	fitValues.bg.pol5.a2.value	= fit5->GetParameter(2);
	fitValues.bg.pol5.a2.error	= fit5->GetParError(2);
	fitValues.bg.pol5.a3.value	= fit5->GetParameter(3);
	fitValues.bg.pol5.a3.error	= fit5->GetParError(3);
	fitValues.bg.pol5.a4.value	= fit5->GetParameter(4);
	fitValues.bg.pol5.a4.error	= fit5->GetParError(4);
	fitValues.bg.pol5.a5.value	= fit5->GetParameter(5);
	fitValues.bg.pol5.a5.error	= fit5->GetParError(5);
}
void	FitBinsSim(const TH1D* dataSignal, const TH1D* dataBG, int channel, int size, FitValues& fitValues)
{
	cout << "channel: " << channel << endl;
	FitBinsSimGaus(dataBG, channel, size, fitValues, false);
	FitBinsSimPol(dataBG, channel, size, fitValues);
	FitBinsSimGaus(dataSignal, channel, size, fitValues, true);
}

void	FitBinsSim(const TH2D* dataSignal, const TH2D* dataBG, const TFile* out, std::vector<int>& ch, std::vector<int>& size, const double signalScale, FitValues* fitValues)
{
	TCanvas*	can = new TCanvas("canFitBinsSim", "FitBinsSim", 1500, 800);
	std::cout << "can size: " << TMath::Ceil(TMath::Sqrt(ch.size())) << std::endl;
    can->Divide(TMath::Ceil(TMath::Sqrt(ch.size())), TMath::Ceil(TMath::Sqrt(ch.size())));
    
	for(int i=0; i<ch.size(); i++)
	{
		can->cd(i+1);
		TH1D*	dataSignalProjX	= ((TH1D*)dataSignal->ProjectionX(TString("SigBin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1))->Clone();
		dataSignalProjX.Scale(signalScale);
		FitBinsSim( dataSignalProjX,
					(TH1D*)dataBG->ProjectionX(TString("BGBin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1),
					ch[i],
					size[i],
					fitValues[i]);
	}
	
	TCanvas*	canMain = new TCanvas("canFitBinsSimMain", "FitBinsSimMain", 1500, 800);
    canMain->Divide(3,3);
    
	canMain->cd(1);
	FitBinsSimGaus((TH1D*)dataSignal->ProjectionX(TString("SignalBinSimMain").Data(), 1, 47), 0, 46, fitValues[ch.size()], true);
	canMain->cd(2);
	FitBinsSimGaus((TH1D*)dataBG->ProjectionX(TString("BGBinGausSimMain").Data(), 1, 47), 0, 46, fitValues[ch.size()], false);
	canMain->cd(3);
	FitBinsSimGaus((TH1D*)dataBG->ProjectionX(TString("BGBinGausSimMain2").Data(), 1, 47), 0, 46, fitValues[ch.size()], false);
	FitBinsSimPol((TH1D*)dataBG->ProjectionX(TString("BGBinPolSimMain").Data(), 1, 47), 0, 46, fitValues[ch.size()]);
	
	double*	x = new double[ch.size()];
	double*	dx = new double[ch.size()];
	canMain->cd(4);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			x[i]	= 0;
			for(int k=0; k<size[i]; k++)
				x[i]	+=	TaggE[ch[i]+k];
			x[i]	/= size[i];
			dx[i]	= 0;
			
			y[i]	= fitValues[i].signal.count;
			dy[i]	= TMath::Sqrt(fitValues[i].signal.count);
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("Entries Signal");
		graph->Draw();
	}
	
	canMain->cd(5);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].bg.count;
			dy[i]	= TMath::Sqrt(fitValues[i].bg.count);
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("Entries BG");
		graph->Draw();
	}
	
	canMain->cd(6);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.count/fitValues[i].bg.count;
			dy[i]	= y[i]*TMath::Sqrt((1/fitValues[i].signal.count) + (1/fitValues[i].bg.count));
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("Entries relation (signal/bg)");
		graph->Draw();
	}
	
	canMain->cd(7);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.factor.value;
			dy[i]	= fitValues[i].signal.factor.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("signal factor");
		graph->Draw();
	}
	
	canMain->cd(8);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.mean.value;
			dy[i]	= fitValues[i].signal.mean.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("signal mean");
		graph->Draw();
	}
	
	canMain->cd(9);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.sigma.value;
			dy[i]	= fitValues[i].signal.sigma.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("signal sigma");
		graph->Draw();
	}
	
	out->cd();
	canMain->Write();
	can->Write();
}


void	FitBins(const TH1D* data, int channel, int size, const FitValues& fitValues, FitValuesData& fitResult)
{
	TF1*	fit = new TF1("fitfkt", "gaus(0) + gaus(3)", 850, 1050);
	fit->SetParameters(data->GetMaximum(), fitValues.signal.mean.value, fitValues.signal.sigma.value, 0, fitValues.bg.gaus.mean.value, fitValues.bg.gaus.sigma.value);
	fit->SetParLimits(0, 0, data->GetMaximum());
	fit->SetParLimits(1, 0.9*fitValues.signal.mean.value, 1.1*fitValues.signal.mean.value);
	fit->SetParLimits(2, 0.9*fitValues.signal.sigma.value, 1.1*fitValues.signal.sigma.value);
	fit->SetParLimits(3, 0, data->GetMaximum());
	fit->SetParLimits(4, 0.9*fitValues.bg.gaus.mean.value, 1.1*fitValues.bg.gaus.mean.value);
	fit->SetParLimits(5, 0.9*fitValues.bg.gaus.sigma.value, 1.1*fitValues.bg.gaus.sigma.value);
	
	TString	title("E_{#gamma} (");
	if(channel==0)
	{
		title.Append(TString().Itoa(int(TaggE[size]-((TaggE[size]-TaggE[size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[0]+((TaggE[0]-TaggE[1])/2)), 10));
		title.Append(" MeV) ");
	}
	else if(channel+size>=47)
	{
		title.Append(TString().Itoa(int(TaggE[47]-((TaggE[46]-TaggE[47])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	else
	{
		title.Append(TString().Itoa(int(TaggE[channel+size]-((TaggE[channel+size]-TaggE[channel+size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	data->SetTitle(title.Data());
	data->GetXaxis()->SetTitle("IM #gamma#gamma [MeV]");	
	data->Fit(fit, "R0");
	data->SetStats(0);		
    data->SetAxisRange(850, 1050);
	data->Draw();
	fit->Draw("SAME");
	
	TF1*	fitSignal = new TF1("fitfktSignal", "gaus(0)", 850, 1050);
	TF1*	fitBG = new TF1("fitfktBG", "gaus(0)", 850, 1050);
	fitSignal->SetLineColor(kRed);
	fitBG->SetLineColor(kGreen);
	fitSignal->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));
	fitBG->SetParameters(fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5));
	fitSignal->Draw("SAME");
	fitBG->Draw("SAME");
	
	fitResult.signal.factor.value	= fit->GetParameter(0);
	fitResult.signal.factor.error	= fit->GetParError(0);
	fitResult.signal.mean.value		= fit->GetParameter(1);
	fitResult.signal.mean.error		= fit->GetParError(1);
	fitResult.signal.sigma.value	= fit->GetParameter(2);
	fitResult.signal.sigma.error	= fit->GetParError(2);
	fitResult.bg.factor.value	= fit->GetParameter(3);
	fitResult.bg.factor.error	= fit->GetParError(3);
	fitResult.bg.mean.value		= fit->GetParameter(4);
	fitResult.bg.mean.error		= fit->GetParError(4);
	fitResult.bg.sigma.value	= fit->GetParameter(5);
	fitResult.bg.sigma.error	= fit->GetParError(5);
}
void	FitBins(const TH2D* data, const TFile* out, std::vector<int>& ch, std::vector<int>& size, const double signalScale, FitValues* fitValues)
{
	TCanvas*	can = new TCanvas("canFitBins", "FitBins", 1500, 800);
    can->Divide(TMath::Ceil(TMath::Sqrt(ch.size())), TMath::Ceil(TMath::Sqrt(ch.size())));
    
    FitValuesData*	fitResult	= new FitValuesData[ch.size()+1];
    
	for(int i=0; i<ch.size(); i++)
	{
		can->cd(i+1);
		FitBins((TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1),
				ch[i],
				size[i],
				fitValues[i],
				fitResult[i]);
	}
	
	
	TCanvas*	canMain = new TCanvas("canFitBinsMain", "FitBinsMain", 1500, 800);
    canMain->Divide(3,3);
    
	canMain->cd(1);
	FitBins((TH1D*)data->ProjectionX(TString("SignalBinMain").Data(), 1, 38), 0, 37, fitValues[ch.size()], fitResult[ch.size()]);
	
	double*	x = new double[ch.size()];
	double*	dx = new double[ch.size()];
	canMain->cd(4);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			x[i]	= 0;
			for(int k=0; k<size[i]; k++)
				x[i]	+=	TaggE[ch[i]+k];
			x[i]	/= size[i];
			dx[i]	= 0;
			
			y[i]	= fitResult[i].signal.factor.value;
			dy[i]	= fitResult[i].signal.factor.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("factor Signal");
		graph->Draw();
	}
	
	canMain->cd(5);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].signal.mean.value;
			dy[i]	= fitResult[i].signal.mean.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("mean Signal");
		graph->Draw();
	}
	
	canMain->cd(6);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].signal.sigma.value;
			dy[i]	= fitResult[i].signal.sigma.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("sigma Signal");
		graph->Draw();
	}
	
	canMain->cd(7);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].bg.factor.value;
			dy[i]	= fitResult[i].bg.factor.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("factor bg");
		graph->Draw();
	}
	
	canMain->cd(8);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].bg.mean.value;
			dy[i]	= fitResult[i].bg.mean.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("mean bg");
		graph->Draw();
	}
	
	canMain->cd(9);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].bg.sigma.value;
			dy[i]	= fitResult[i].bg.sigma.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("sigma bg");
		graph->Draw();
	}
	
	out->cd();
	can->Write();
	canMain->Write();
	
	if(fitResult)	delete fitResult;
}

void	FitBins(const TFile* dataFile, const TFile* mcSinalFile, const TFile* mcBGFile, const TFile* scalerFile, const TFile* out, const Int_t addedChannels, const char* histName, const double signalScale)
{
	std::vector<int>	ch;
	std::vector<int>	size;
	while(addedChannels*ch.size()<47)
	{
		if(47-(addedChannels*ch.size())>=addedChannels)
			size.push_back(addedChannels);
		else
			size.push_back(47-(addedChannels*ch.size()));
		ch.push_back(addedChannels*ch.size());
		/*if(ch.back()==0)
		{
			std::cout << "taggedEnergyMin: " << TaggE[0]<< std::endl;
			std::cout << "taggedEnergyMin: " << TaggE[1]<< std::endl;
			std::cout << "taggedEnergyMin: " << (TaggE[1]-TaggE[0])/2 << std::endl;
			std::cout << "taggedEnergyMin: " << TaggE[0]-((TaggE[1]-TaggE[0])/2) << std::endl;
			taggedEnergyMin.puch_back(TaggE[0]-((TaggE[1]-TaggE[0])/2));
		}
		else
			taggedEnergyMin.puch_back(TaggE[ch.back()]-((TaggE[ch.back()]-TaggE[ch.back()-1])/2));
			
		if(47-(addedChannels*ch.size())>=addedChannels)
			taggedEnergyMax.puch_back(TaggE[ch.back()]+((TaggE[ch.back()+1]-TaggE[ch.back()])/2));
		else
			taggedEnergyMax.puch_back(TaggE[46]+((TaggE[46]-TaggE[45])/2));*/
			
		std::cout << "starting number: " << ch.back() << "     channel count: " << size.back() << std::endl;
	}
    
    FitValues*	fitValues	= new FitValues[ch.size()+1];
	TH2D*		dataSignal	= (TH2D*)mcSinalFile->Get(histName);
	TH2D*		dataBG		= (TH2D*)mcBGFile->Get(histName);
    
	FitBinsSim(dataSignal, dataBG, out, ch, size, signalScale, fitValues);
	
	TH2D*		data		= (TH2D*)dataFile->Get(histName);
	FitBins(data, out, ch, size, signalScale, fitValues);
	
	if(fitValues)	delete fitValues;
}


void	FitBins(const char* dataFileName, const char* mcSignalFileName, const char* mcBGFileName, const char* scalerFileName, const Int_t addedChannels, const char* histName, const double	signalScale)
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
	TFile*	mcBGFile			= TFile::Open(mcBGFileName);
	if(!mcBGFile)
	{
		std::cout << "Can not open mcBGFile " << mcBGFileName << std::endl;
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
	
	FitBins(dataFile, mcSignalFile, mcBGFile, scalerFile, out, addedChannels, histName, signalScale);
}

void	FitBins(const char* dir = ".", const Int_t addedChannels = 3, const char* histName = "WithProton/fitProton6/TaggerBinning/IM_Bins")
{
	std::strstream	acquSignalFileName;
	std::strstream	acquBGFileName;
	
	acquSignalFileName << dir << "/Acqu_g4_sim_etap_pi0pi0eta_00.root";
	acquBGFileName << dir << "/Acqu_g4_sim_pi0pi0pi0_6g_00.root";
	TFile*	mcSignalFile	= TFile::Open(acquSignalFileName.str().c_str());
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << acquSignalFileName << std::endl;
		return;
	}
	TFile*	mcBGFile		= TFile::Open(acquBGFileName.str().c_str());
	if(!mcBGFile)
	{
		std::cout << "Can not open mcBGFile " << acquBGFileName << std::endl;
		return;
	}
	double	signalScale	= BestFitSignalScale(mcSignalFile, mcBGFile);
	cout << "signalScale:   " << signalScale << endl;
	
	
	std::strstream	dataFileName;
	std::strstream	mcSignalFileName;
	std::strstream	mcBGFileName;
	std::strstream	scalerFileName;
	
	dataFileName << dir << "/Result_CB.root";
	mcSignalFileName << dir << "/Phys_g4_sim_etap_pi0pi0eta_00.root";
	mcBGFileName << dir << "/Phys_g4_sim_pi0pi0pi0_6g_00.root";
	scalerFileName << dir << "/ScalerPhysics_CB.root";
	
	cout << "FitBins:" << endl;
	FitBins(dataFileName.str().c_str(), mcSignalFileName.str().c_str(), mcBGFileName.str().c_str(), scalerFileName.str().c_str(), addedChannels, histName, signalScale);
}
