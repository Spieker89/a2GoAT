#ifndef __GH1Example_h__
#define __GH1Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GHistEvent.h"
#include "GCheckProton.h"

class	GH1Example  : public GTreeManager
{
private:
    GHistEvent3Mesons   hist_eta;
    GCheckProton        check_eta_proton;
    GHistEvent3Mesons   hist_etap;
    GCheckProton        check_etap_proton;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    GH1Example();
    virtual ~GH1Example();

    //virtual Bool_t	Init(const char* configfile);

};
#endif
