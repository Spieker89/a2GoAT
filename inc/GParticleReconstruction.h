
#ifndef __GParticleReconstruction_h__
#define __GParticleReconstruction_h__

#include "GTreeManager.h"
#define DEFAULT_PI0_IM_WIDTH 20.0
#define DEFAULT_ETA_IM_WIDTH 44.0
#define DEFAULT_ETAP_IM_WIDTH 60.0

#define	pdg_rootino 0
#define	pdg_pi0 1
#define	pdg_eta 2
#define	pdg_etaP 3
#define	pdg_proton 4
#define	pdg_chpion 5
#define	pdg_electron 6

#define m_pi0 135.0
#define	m_eta 545.0
#define	m_etaP 958.0
#define	m_proton 938.3
#define	m_chpion 139.6
#define	m_electron 0.511

class	GParticleReconstruction : public GTreeManager
{
private:

	Bool_t 		ReconstructChargedParticles;
	Bool_t		ReconstructMesons;

	char 		cutfilename[256];
	char 		cutname[256];	
	
	TFile* 		CutFile;
	TCutG* 		Cut;
   	TCutG* 		Cut_CB_proton; 
   	TCutG* 		Cut_CB_pion; 
   	TCutG*		Cut_CB_electron;
   	TCutG* 		Cut_TAPS_proton; 
   	TCutG* 		Cut_TAPS_pion; 
   	TCutG*		Cut_TAPS_electron;

   	Int_t		Cut_CB_proton_active;   	
   	Int_t		Cut_TAPS_proton_active;  
   	Int_t 		Cut_proton_active;
   	
   	Int_t		Cut_CB_pion_active;
   	Int_t		Cut_TAPS_pion_active;
   	Int_t 		Cut_pion_active;
   	
   	Int_t		Cut_CB_electron_active; 
   	Int_t		Cut_TAPS_electron_active; 
   	Int_t 		Cut_electron_active;
   	   	   	
   	TCutG* 		Cut_proton;
   	TCutG*		Cut_pion;
   	TCutG*		Cut_electron;
	
	Double_t	width_pi0;
	Double_t	width_eta;
	Double_t	width_etaP;

	Int_t* 		PDG;
	Int_t* 		Identified;
	
	Int_t 		nParticles;
   	Int_t 		i;
   	
protected:
	    
public:

	GParticleReconstruction();
	~GParticleReconstruction();

	virtual Bool_t	Init(const char* intreefile, const char* outtreefile) {return kTRUE;}
	virtual void	Analyse() {;}
	virtual void	Reconstruct();

    
	Bool_t	PostInit();
	void 	ChargedReconstruction();
	void 	MesonReconstruction();	
	void	AddParticle(Int_t pdg_code, Int_t nindex, Int_t index_list[]);
	void	AddParticle(Int_t pdg_code, Int_t i) 
			{Int_t index_list[1]; index_list[0] = i; AddParticle(pdg_code, 1, index_list);}

	TCutG*	OpenCutFile(Char_t* filename, Char_t* cutname);

	// Make some variables available for sorting/printing
 

};


#endif

