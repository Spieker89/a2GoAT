#ifndef __GH1Example_h__
#define __GH1Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GHistEvent.h"

class	GH1Example  : public GTreeManager
{
private:
    GHistEvent     hist_eta;
    GHistEvent     hist_etap;

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
