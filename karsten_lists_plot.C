#include <iostream>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>


using namespace std;
#include <fstream>
int karsten_lists_plot(const char *name1, const char *name2) {

 	const double DEFAULT_VALUE = -1.0;
  //************************************
  

	ifstream in3;
 	in3.open("/disk/nobackup/001/spieker/Mainz1/a2GoAT/runNrliste_may.txt");
				Int_t runNumber_intern;

	if (in3.is_open()) 
	{
		while (1) 
		{
			in3 >> runNumber_intern;

			TFile fout(Form("test_%i.root",runNumber_intern),"RECREATE");
			TTree tout("tout","correation");
			int evN;           tout.Branch("evN", &evN);
			int isPID, isMWPC;      tout.Branch("isPID", &isPID);  tout.Branch("isMWPC", &isMWPC);
			double theta_PID, theta_MWPC; tout.Branch("theta_PID",&theta_PID); tout.Branch("theta_MWPC",&theta_MWPC);
			int detector_PID, detector_MWPC; tout.Branch("detector_PID", &detector_PID); tout.Branch("detector_MWPC", &detector_MWPC);
			int runNumber1, runNumber2; tout.Branch("runNumber1", &runNumber1); tout.Branch("runNumber2", &runNumber2);
			Float_t theta;
			Int_t evNumber; 
			Int_t runNumber;
			Int_t detector;

			Int_t bigg=0;
			Int_t smll=1000000000;

			ifstream in1;
			in1.open(Form("liste_PID_May14_prompt_%i.txt",runNumber_intern));
			ifstream in2;
			in2.open(Form("liste_PID_May14_old_prompt_%i.txt",runNumber_intern));
			Int_t N11=0;
			Int_t N22=0;

			if (in1.is_open()) {
				while (1) {
					
					in1 >> theta >> evNumber >> runNumber >> detector;
					if(evNumber==0)cout << evNumber << endl;
					if (!in1.good()) break;
					if(evNumber>bigg)bigg=evNumber;
					if(evNumber<smll)smll=evNumber;
						
					
					N11++;
					}
			}
			in1.close();

			if (in2.is_open()) {
				while (1) {
				
					in2 >> theta >> evNumber >> runNumber >> detector;
					if (!in2.good()) break;
					if(evNumber>bigg)bigg=evNumber;
					if(evNumber<smll)smll=evNumber;
					
				
					N22++;
					}
				}
			in2.close();
			
			cout << "runNr = " <<runNumber_intern << "\t" << "N1 = " << N11 << ", N2 = " << N22 << endl;

			const Int_t Np=bigg-smll;
			cout << bigg << "\t" << smll << "\t" << Np << endl;
			double *val1 = new double[Np+1];
			double *val2 = new double[Np+1];
			double *runNr1 = new double[Np+1];
			double *runNr2 = new double[Np+1];
			double *detector_PID_array = new double[Np+1];
			double *detector_MWPC_array = new double[Np+1];
			for(int i=0;i<Np+1;i++) {
				val1[i] = DEFAULT_VALUE; val2[i] = DEFAULT_VALUE;
			} 
			
				in1.open(Form("liste_PID_May14_prompt_%i.txt",runNumber_intern));
				if (in1.is_open()) {
					while (1) {
					
					in1 >> theta >> evNumber >> runNumber >> detector;
					if (!in1.good()) break;
					evN=evNumber;
					val1[evN-smll]=theta;
					detector_PID_array[evN-smll]=detector;
					runNr1[evN-smll]=runNumber;
			
			
					}
				}
				in1.close();
					
				in2.open(Form("liste_PID_May14_old_prompt_%i.txt",runNumber_intern));
				if (in2.is_open()) {
					while (1) {
					
					in2 >> theta >> evNumber >> runNumber >> detector;
					if (!in2.good()) break;
					evN=evNumber;
					val2[evN-smll]=theta;
					detector_MWPC_array[evN-smll]=detector;
					runNr2[evN-smll]=runNumber;
			
					}
				}
				in2.close();
			
			
			for(int i=0;i<Np+1;i++) {
			
				evN = smll+i;
				isPID=0; isMWPC=0;
				if(val1[i]<0 && val2[i]<0) continue;
				theta_PID=val1[i]; theta_MWPC=val2[i];
				detector_PID=detector_PID_array[i]; detector_MWPC=detector_MWPC_array[i]; 
				runNumber1=runNr1[i]; runNumber2=runNr2[i]; 
				if(runNumber1==0)runNumber1=runNr2[i];
				if(runNumber2==0)runNumber2=runNr1[i];
				if(val1[i]>0) isPID=1;
				if(val2[i]>0) isMWPC=1;
				tout.Fill();
			}

		if (!in3.good()) break;

		tout.Write();
		fout.Close();
		}	
	}

  return 0;
}
