#ifndef __MyPhysics_h__
#define __MyPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GAnalysis3Mesons.h"

class	MyPhysics  : public GTreeManager
{
private:
    GAnalysis3Mesons        hist;
    GAnalysis3MesonsProton  hist_proton;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    MyPhysics();
    virtual ~MyPhysics();

    virtual Bool_t	Init(const char* configfile);

};
#endif
