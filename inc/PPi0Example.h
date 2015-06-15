#ifndef __PPi0Example_h__
#define __PPi0Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"

#include "HistoManu.h"
#include "HistoManu2.h"
#include "HistoManu3.h"
class	PPi0Example  : public PPhysics
{
private:

     
    HistoManu*	IM;
    HistoManu*	MM;
    HistoManu2* MM_energy;
    HistoManu3* MM_energy_inv;

	TH1F* MM1;
	TH2F* MM1_energy;
    	TH3F* MM1_energy_inv;
	TH1F *IM1;
	TH1F *IM1_side;


protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PPi0Example();
    virtual ~PPi0Example();
    virtual Bool_t  Init();

};
#endif
