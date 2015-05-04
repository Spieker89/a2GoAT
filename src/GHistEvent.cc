#include "GHistEvent.h"
#include "GTreeTagger.h"
#include "GTreeMeson.h"


GHistEventAngles::GHistEventAngles(const char* name, const char* title, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    theta(TString(name).Append("_theta"), TString(title).Append(" Theta"), 1800, 0, 180, kFALSE),
    phi(TString(name).Append("_phi"), TString(title).Append(" Phi"), 3600, -180, 180, kFALSE),
    thetaCM(TString(name).Append("_thetaCM"), TString(title).Append(" Theta CM"), 1800, 0, 180, 48, kFALSE)
{

}

GHistEventAngles::~GHistEventAngles()
{

}


void    GHistEventAngles::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        theta.PrepareWriteList(arr, TString(name).Append("_Theta").Data());
        phi.PrepareWriteList(arr, TString(name).Append("_Phi").Data());
        thetaCM.PrepareWriteList(arr, TString(name).Append("_ThetaCM").Data());
    }
    else
    {
        theta.PrepareWriteList(arr);
        phi.PrepareWriteList(arr);
        thetaCM.PrepareWriteList(arr);
    }
}








GHistEvent::GHistEvent(const char* name, const char* title, Bool_t linkHistogram) :
    GHistEventAngles(name, title, linkHistogram),
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

    GHistEventAngles::PrepareWriteList(arr, name);

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
    sub0_im(TString(name).Append("_sub0im"), TString(title).Append(" sub Part. 0 inv. Mass"), 800, 0, 800, kFALSE),
    sub1_im(TString(name).Append("_sub1im"), TString(title).Append(" sub Part. 1 inv. Mass"), 400, 0, 400, kFALSE),
    sub2_im(TString(name).Append("_sub2im"), TString(title).Append(" sub Part. 2 inv. Mass"), 400, 0, 400, kFALSE)
{

}

GHistEvent3Mesons::~GHistEvent3Mesons()
{

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
















GHistEventProton::GHistEventProton(const char* name, const char* title, Bool_t linkHistogram) :
    GHistEvent(name, title, linkHistogram),
    protonE(TString(name).Append("_protonE"), TString(title).Append(" Proton E"), 800, 0, 800, kFALSE),
    protonAngles(TString(name).Append("_protonAn"), TString(title).Append(" ProtonAn"), kFALSE)
{

}

GHistEventProton::~GHistEventProton()
{

}


void    GHistEventProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistEvent::PrepareWriteList(arr, name);

    if(name)
    {
        protonE.PrepareWriteList(arr, TString(name).Append("_protonE").Data());
        protonAngles.PrepareWriteList(arr, TString(name).Append("_protonEcalc").Data());
    }
    else
    {
        protonE.PrepareWriteList(arr);
        protonAngles.PrepareWriteList(arr);
    }
}





GHistEvent3MesonsProton::GHistEvent3MesonsProton(const char* name, const char* title, Bool_t linkHistogram) :
    GHistEventProton(name, title, linkHistogram),
    sub0_im(TString(name).Append("_sub0im"), TString(title).Append(" sub Part. 0 inv. Mass"), 800, 0, 800, kFALSE),
    sub1_im(TString(name).Append("_sub1im"), TString(title).Append(" sub Part. 1 inv. Mass"), 400, 0, 400, kFALSE),
    sub2_im(TString(name).Append("_sub2im"), TString(title).Append(" sub Part. 2 inv. Mass"), 400, 0, 400, kFALSE)
{

}

GHistEvent3MesonsProton::~GHistEvent3MesonsProton()
{

}


void    GHistEvent3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistEventProton::PrepareWriteList(arr, name);

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
