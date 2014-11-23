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


void    GHistEvent::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        mm.PrepareWriteList(arr, TString(name).Append("_MM").Data());
    }
    else
    {
        im.PrepareWriteList(arr);
        mm.PrepareWriteList(arr);
    }
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

void    GHistEvent3Mesons::Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t* SUB_THETA, const Double_t* SUB_PHI, const Double_t taggerTime)
{
    GHistEvent::Fill(IM, MM, taggerTime);
    sub0_im.Fill(SUB0_IM, taggerTime);
    sub1_im.Fill(SUB1_IM, taggerTime);
    sub2_im.Fill(SUB2_IM, taggerTime);
    for(int i=0; i<6; i++)
    {
        sub_theta.Fill(SUB_THETA[i], taggerTime);
        sub_theta.Fill(SUB_PHI[i], taggerTime);
    }
}
void    GHistEvent3Mesons::Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t* SUB_THETA, const Double_t* SUB_PHI, const Double_t taggerTime, const Int_t taggerChannel)
{
    GHistEvent::Fill(IM, MM, taggerTime, taggerChannel);
    sub0_im.Fill(SUB0_IM, taggerTime);
    sub1_im.Fill(SUB1_IM, taggerTime);
    sub2_im.Fill(SUB2_IM, taggerTime);
    for(int i=0; i<6; i++)
    {
        sub_theta.Fill(SUB_THETA[i], taggerTime);
        sub_theta.Fill(SUB_PHI[i], taggerTime);
    }
}

void    GHistEvent3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistEvent::PrepareWriteList(arr, name);

    if(name)
    {
        sub0_im.PrepareWriteList(arr, TString(name).Append("_sub0IM").Data());
        sub1_im.PrepareWriteList(arr, TString(name).Append("_sub1IM").Data());
        sub2_im.PrepareWriteList(arr, TString(name).Append("_sub2IM").Data());
    }
    else
    {
        sub0_im.PrepareWriteList(arr);
        sub1_im.PrepareWriteList(arr);
        sub2_im.PrepareWriteList(arr);
    }
}
