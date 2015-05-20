#ifndef __PPi0Example_h__
#define __PPi0Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"

class	PPi0Example  : public PPhysics
{
private:

     
    GH1*	IM;
    GH1*	MM;
	TH1F* MM1;
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
