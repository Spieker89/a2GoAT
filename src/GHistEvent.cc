#include "GHistEvent.h"
#include "GTreeTagger.h"
#include "GTreeMeson.h"




GHistEvent::GHistEvent(const char* name, const char* title, const char *dirName) :
    im(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 1500, 0, 1500, kTRUE, dirName),
    mm(TString(name).Append("_mm"), TString(title).Append(" mis. Mass"), 2000, 0, 2000, kTRUE, dirName)
{

}

GHistEvent::~GHistEvent()
{

}

void    GHistEvent::Fill(const Double_t IM, const Double_t MM, const Double_t taggerTime, const Double_t taggerChannel)
{
    im.Fill(IM, taggerTime, taggerChannel);
    mm.Fill(MM, taggerTime, taggerChannel);
}

void    GHistEvent::Fill(const TLorentzVector& part, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
        {
            im.Fill(part.M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            mm.Fill((tagger.GetVectorProtonTarget(i)-part).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
        else
        {
            im.Fill(part.M(), tagger.GetTagged_t(i), 0);
            mm.Fill((tagger.GetVectorProtonTarget(i)-part).M(), tagger.GetTagged_t(i), 0);
        }
    }
}

void    GHistEvent::Fill(const TLorentzVector& part, const TLorentzVector& rest, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
        {
            im.Fill(part.M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            mm.Fill((tagger.GetVectorProtonTarget(i)-part-rest).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
        else
        {
            im.Fill(part.M(), tagger.GetTagged_t(i), 0);
            mm.Fill((tagger.GetVectorProtonTarget(i)-part-rest).M(), tagger.GetTagged_t(i), 0);
        }
    }
}

void    GHistEvent::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    mm.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}












GHistEvent3Mesons::GHistEvent3Mesons(const char* name, const char* title, const char* dirName) :
    GHistEvent(name, title, dirName),
    sub0_im(TString(name).Append("_sub0im"), TString(title).Append(" sub Part. 0 inv. Mass"), 800, 0, 800, kTRUE, dirName),
    sub1_im(TString(name).Append("_sub1im"), TString(title).Append(" sub Part. 1 inv. Mass"), 400, 0, 400, kTRUE, dirName),
    sub2_im(TString(name).Append("_sub2im"), TString(title).Append(" sub Part. 2 inv. Mass"), 400, 0, 400, kTRUE, dirName)
{

}

GHistEvent3Mesons::~GHistEvent3Mesons()
{

}

void    GHistEvent3Mesons::Fill(const Double_t IM, const Double_t MM, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime, const Double_t taggerChannel)
{
    sub0_im.Fill(SUB0_IM, taggerTime, taggerChannel);
    sub1_im.Fill(SUB1_IM, taggerTime, taggerChannel);
    sub2_im.Fill(SUB2_IM, taggerTime, taggerChannel);
    GHistEvent::Fill(IM, MM, taggerTime, taggerChannel);
}

void    GHistEvent3Mesons::Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
        {
            im.Fill(meson.Particle(0).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            mm.Fill((tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            sub0_im.Fill((meson.SubPhotons(0, 0)+meson.SubPhotons(0, 1)).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            sub1_im.Fill((meson.SubPhotons(0, 2)+meson.SubPhotons(0, 3)).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            sub2_im.Fill((meson.SubPhotons(0, 4)+meson.SubPhotons(0, 5)).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
        else
        {
            im.Fill(meson.Particle(0).M(), tagger.GetTagged_t(i));
            mm.Fill((tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M(), tagger.GetTagged_t(i));
            sub0_im.Fill((meson.SubPhotons(0, 0)+meson.SubPhotons(0, 1)).M(), tagger.GetTagged_t(i));
            sub1_im.Fill((meson.SubPhotons(0, 2)+meson.SubPhotons(0, 3)).M(), tagger.GetTagged_t(i));
            sub2_im.Fill((meson.SubPhotons(0, 4)+meson.SubPhotons(0, 5)).M(), tagger.GetTagged_t(i));
        }
    }
}

void    GHistEvent3Mesons::Fill(const GTreeMeson& meson, const TLorentzVector& rest, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
        {
            im.Fill(meson.Particle(0).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            mm.Fill((tagger.GetVectorProtonTarget(i)-meson.Particle(0)-rest).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            sub0_im.Fill((meson.SubPhotons(0, 0)+meson.SubPhotons(0, 1)).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            sub1_im.Fill((meson.SubPhotons(0, 2)+meson.SubPhotons(0, 3)).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            sub2_im.Fill((meson.SubPhotons(0, 4)+meson.SubPhotons(0, 5)).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
        else
        {
            im.Fill(meson.Particle(0).M(), tagger.GetTagged_t(i));
            mm.Fill((tagger.GetVectorProtonTarget(i)-meson.Particle(0)-rest).M(), tagger.GetTagged_t(i));
            sub0_im.Fill((meson.SubPhotons(0, 0)+meson.SubPhotons(0, 1)).M(), tagger.GetTagged_t(i));
            sub1_im.Fill((meson.SubPhotons(0, 2)+meson.SubPhotons(0, 3)).M(), tagger.GetTagged_t(i));
            sub2_im.Fill((meson.SubPhotons(0, 4)+meson.SubPhotons(0, 5)).M(), tagger.GetTagged_t(i));
        }
    }
}

void    GHistEvent3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    sub0_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub1_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub2_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    GHistEvent::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
