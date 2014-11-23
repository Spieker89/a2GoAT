#ifndef __GHistEvent_h__
#define __GHistEvent_h__

#include <TLorentzVector.h>

#include "GH1.h"


class	GHistEvent  : public GHistLinked
{
private:

protected:
    GH1     im;
    GH1     mm;
    GH1     theta;
    GH1     phi;

public:
    GHistEvent(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistEvent();

    virtual void    CalcResult()                                                                                                                                {im.CalcResult(); mm.CalcResult();}
    virtual Int_t   Fill(Double_t x)                                                                                                                            {}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime)                                                                       {im.Fill(IM, taggerTime); mm.Fill(MM, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime, const Int_t taggerChannel)                                            {im.Fill(IM, taggerTime, taggerChannel); mm.Fill(MM, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Theta, const Double_t Phi, const Double_t taggerTime)                             {im.Fill(IM, taggerTime); mm.Fill(MM, taggerTime); theta.Fill(IM, taggerTime); phi.Fill(MM, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Theta, const Double_t Phi, const Double_t taggerTime, const Int_t taggerChannel)  {im.Fill(IM, taggerTime, taggerChannel); mm.Fill(MM, taggerTime); theta.Fill(IM, taggerTime); phi.Fill(MM, taggerTime);}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "")                                                                                                                {im.Reset(option); mm.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE)                           {im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); mm.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);}
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};




class	GHistEvent3Mesons   : public GHistEvent
{
private:
    GH1     sub0_im;
    GH1     sub1_im;
    GH1     sub2_im;
    GH1     sub_theta;
    GH1     sub_phi;

protected:

public:
    GHistEvent3Mesons(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistEvent3Mesons();

    virtual void    CalcResult()                                                                                                                                                                {GHistEvent::CalcResult(); sub0_im.CalcResult(); sub1_im.CalcResult(); sub2_im.CalcResult();}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime)                               {GHistEvent::Fill(IM, MM, taggerTime); sub0_im.Fill(SUB0_IM, taggerTime); sub1_im.Fill(SUB1_IM, taggerTime); sub2_im.Fill(SUB2_IM, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime, const Int_t taggerChannel)    {GHistEvent::Fill(IM, MM, taggerTime, taggerChannel); sub0_im.Fill(SUB0_IM, taggerTime); sub1_im.Fill(SUB1_IM, taggerTime); sub2_im.Fill(SUB2_IM, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t* SUB_THETA, const Double_t* SUB_PHI, const Double_t taggerTime);
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t* SUB_THETA, const Double_t* SUB_PHI, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option)                                                                                                                                                     {GHistEvent::Reset(option); sub0_im.Reset(option); sub1_im.Reset(option); sub2_im.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)                                                                    {GHistEvent::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); sub0_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); sub1_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); sub2_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);}
};

#endif
