#include "GHistEvent.h"
#include "GTreeTagger.h"
#include "GTreeMeson.h"




GHistEvent::GHistEvent(const char* name, const char* title, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    im(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 1500, 0, 1500, 48, kFALSE),
    mm(TString(name).Append("_mm"), TString(title).Append(" mis. Mass"), 2000, 0, 2000, 48, kFALSE)
{

}

GHistEvent::~GHistEvent()
{

}

void   GHistEvent::CalcResult()
{
    im.CalcResult();
    mm.CalcResult();
}

void    GHistEvent::Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime)
{
    im.Fill(IM, taggerTime);
    mm.Fill(MM, taggerTime);
}

void    GHistEvent::Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime, const Int_t taggerChannel)
{
    im.Fill(IM, taggerTime, taggerChannel);
    mm.Fill(MM, taggerTime);
}

void    GHistEvent::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
    mm.PrepareWriteList(arr, TString(name).Append("_MM").Data());
}

void    GHistEvent::Reset(Option_t* option)
{
    im.Reset(option);
    mm.Reset(option);
}

void    GHistEvent::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    mm.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}












GHistEvent3Mesons::GHistEvent3Mesons(const char* name, const char* title, Bool_t linkHistogram) :
    GHistEvent(name, title, linkHistogram),
    sub0_im(TString(name).Append("_sub0im"), TString(title).Append(" sub Part. 0 inv. Mass"), 800, 0, 800, 48, kFALSE),
    sub1_im(TString(name).Append("_sub1im"), TString(title).Append(" sub Part. 1 inv. Mass"), 400, 0, 400, 48, kFALSE),
    sub2_im(TString(name).Append("_sub2im"), TString(title).Append(" sub Part. 2 inv. Mass"), 400, 0, 400, 48, kFALSE)
{

}

GHistEvent3Mesons::~GHistEvent3Mesons()
{

}

void   GHistEvent3Mesons::CalcResult()
{
    GHistEvent::CalcResult();
    sub0_im.CalcResult();
    sub1_im.CalcResult();
    sub2_im.CalcResult();
}

void    GHistEvent3Mesons::Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime)
{
    sub0_im.Fill(SUB0_IM, taggerTime);
    sub1_im.Fill(SUB1_IM, taggerTime);
    sub2_im.Fill(SUB2_IM, taggerTime);
    GHistEvent::Fill(IM, MM, taggerTime);
}

void    GHistEvent3Mesons::Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime, const Int_t taggerChannel)
{
    sub0_im.Fill(SUB0_IM, taggerTime);
    sub1_im.Fill(SUB1_IM, taggerTime);
    sub2_im.Fill(SUB2_IM, taggerTime);
    GHistEvent::Fill(IM, MM, taggerTime, taggerChannel);
}

void    GHistEvent3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistEvent::PrepareWriteList(arr, name);

    sub0_im.PrepareWriteList(arr, TString(name).Append("_sub0IM").Data());
    sub1_im.PrepareWriteList(arr, TString(name).Append("_sub1IM").Data());
    sub2_im.PrepareWriteList(arr, TString(name).Append("_sub2IM").Data());
}

void    GHistEvent3Mesons::Reset(Option_t* option)
{
    GHistEvent::Reset(option);
    sub0_im.Reset(option);
    sub1_im.Reset(option);
    sub2_im.Reset(option);
}

void    GHistEvent3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    GHistEvent::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub0_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub1_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub2_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
