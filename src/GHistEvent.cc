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
            mm.Fill((tagger.GetVectorProtonTarget(0)-part).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
        else
        {
            im.Fill(part.M(), tagger.GetTagged_t(i), 0);
            mm.Fill((tagger.GetVectorProtonTarget(0)-part).M(), tagger.GetTagged_t(i), 0);
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
            mm.Fill((tagger.GetVectorProtonTarget(0)-part-rest).M(), tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
        else
        {
            im.Fill(part.M(), tagger.GetTagged_t(i), 0);
            mm.Fill((tagger.GetVectorProtonTarget(0)-part-rest).M(), tagger.GetTagged_t(i), 0);
        }
    }
}

void    GHistEvent::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    mm.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}


