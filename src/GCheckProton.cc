#include "GCheckProton.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"



GCheckProtonHist::GCheckProtonHist(const char* name, const char* title, const char* dirName) :
    protonAngeDiff(TString(name).Append("_ProtonAngleDiff"), TString(title).Append(" Proton Angle Diff."), 1000, 0, 100, kTRUE, dirName),
    protonAngeDiffSmalest(TString(name).Append("_ProtonAngleDiffSmalest"), TString(title).Append(" Smalest Proton Angle Diff."), 1000, 0, 100, kTRUE, dirName),
    protonCoplanarity(TString(name).Append("_Coplanarity"), TString(title).Append(" Coplanarity"), 3600, 0, 360, kTRUE, dirName)
{

}

GCheckProtonHist::~GCheckProtonHist()
{

}







GCheckProton::GCheckProton(const char* name, const char* title, const char* dirName, const Double_t ProtonAngleDiff_max, const Double_t ProtonCoplanarity_min, const Double_t ProtonCoplanarity_max) :
    raw(TString(name).Append("_Raw"), TString(title).Append(" Raw"), TString(dirName).Append("/CheckProton/Raw")),
    cutProtonAngle(TString(name).Append("_cutProtonAngle"), TString(title).Append(" Cut Proton Angle"), TString(dirName).Append("/CheckProton/CutProtonAngle")),
    cutCoplanarity(TString(name).Append("_cutCoplanarity"), TString(title).Append(" Cut Coplanarity"), TString(dirName).Append("/CheckProton/CutCoplanarity")),
    cutBoth(TString(name).Append("_cutBoth"), TString(title).Append(" Cut Both"), TString(dirName).Append("/CheckProton/CutBoth")),
    CutProtonAngleDiff(ProtonAngleDiff_max)
{
    CutProtonCoplanarity[0] = ProtonCoplanarity_min;
    CutProtonCoplanarity[1] = ProtonCoplanarity_max;
}

GCheckProton::~GCheckProton()
{

}

Bool_t  GCheckProton::Check(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    if(proton.GetNParticles()>1)
        return kFALSE;
    if(tagger.GetNTagged() == 0)
        return kFALSE;

    Bool_t      passedAngleDiff   = kFALSE;
    Double_t    helpAngleDiff;
    Double_t    smalestAngleDiff    = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(0)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
    Double_t    smalestAngleDiffTaggerTime     = tagger.GetTagged_t(0);
    Double_t    smalestAngleDiffTaggerBin      = tagger.GetTagged_ch(0);
    Double_t    helpCoplanarity     = TMath::RadToDeg()*TMath::Abs(meson.Particle(0).Phi()-proton.Particle(0).Phi());

    if(helpCoplanarity>CutProtonCoplanarity[0] && helpCoplanarity<CutProtonCoplanarity[1])
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            if(CreateHistogramsForTaggerBinning)
            {
                helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
                if(helpAngleDiff < smalestAngleDiff)
                {
                    smalestAngleDiff            = helpAngleDiff;
                    smalestAngleDiffTaggerTime  = tagger.GetTagged_t(i);
                    smalestAngleDiffTaggerBin   = tagger.GetTagged_ch(i);
                }
                raw.Fill(helpAngleDiff, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                cutCoplanarity.Fill(helpAngleDiff, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                if(helpAngleDiff<CutProtonAngleDiff)
                {
                    passedAngleDiff = kTRUE;
                    cutProtonAngle.Fill(helpAngleDiff, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    cutBoth.Fill(helpAngleDiff, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                }
            }
            else
            {
                helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
                if(helpAngleDiff < smalestAngleDiff)
                {
                    smalestAngleDiff            = helpAngleDiff;
                    smalestAngleDiffTaggerTime  = tagger.GetTagged_t(i);
                    smalestAngleDiffTaggerBin   = tagger.GetTagged_ch(i);
                }
                raw.Fill(helpAngleDiff, tagger.GetTagged_t(i), 0);
                cutCoplanarity.Fill(helpAngleDiff, tagger.GetTagged_t(i), 0);
                if(helpAngleDiff<CutProtonAngleDiff)
                {
                    passedAngleDiff = kTRUE;
                    cutProtonAngle.Fill(helpAngleDiff, tagger.GetTagged_t(i), 0);
                    cutBoth.Fill(helpAngleDiff, tagger.GetTagged_t(i), 0);
                }
            }
        }
        raw.Fill(smalestAngleDiff, helpCoplanarity, smalestAngleDiffTaggerTime, smalestAngleDiffTaggerBin);
        cutCoplanarity.Fill(smalestAngleDiff, helpCoplanarity, smalestAngleDiffTaggerTime, smalestAngleDiffTaggerBin);
        if(passedAngleDiff == kTRUE)
        {
            cutProtonAngle.Fill(smalestAngleDiff, helpCoplanarity, smalestAngleDiffTaggerTime, smalestAngleDiffTaggerBin);
            cutBoth.Fill(smalestAngleDiff, helpCoplanarity, smalestAngleDiffTaggerTime, smalestAngleDiffTaggerBin);
            return kTRUE;
        }
    }
    else
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            if(CreateHistogramsForTaggerBinning)
            {
                helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
                if(helpAngleDiff < smalestAngleDiff)
                {
                    smalestAngleDiff            = helpAngleDiff;
                    smalestAngleDiffTaggerTime  = tagger.GetTagged_t(i);
                    smalestAngleDiffTaggerBin   = tagger.GetTagged_ch(i);
                }
                raw.Fill(helpAngleDiff, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                if(helpAngleDiff<CutProtonAngleDiff)
                {
                    passedAngleDiff = kTRUE;
                    cutProtonAngle.Fill(helpAngleDiff, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                }
            }
            else
            {
                helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
                if(helpAngleDiff < smalestAngleDiff)
                {
                    smalestAngleDiff            = helpAngleDiff;
                    smalestAngleDiffTaggerTime  = tagger.GetTagged_t(i);
                    smalestAngleDiffTaggerBin   = tagger.GetTagged_ch(i);
                }
                raw.Fill(helpAngleDiff, tagger.GetTagged_t(i), 0);
                if(helpAngleDiff<CutProtonAngleDiff)
                {
                    passedAngleDiff = kTRUE;
                    cutProtonAngle.Fill(helpAngleDiff, tagger.GetTagged_t(i), 0);
                }
            }
        }
        raw.Fill(smalestAngleDiff, helpCoplanarity, smalestAngleDiffTaggerTime, smalestAngleDiffTaggerBin);
        if(passedAngleDiff == kTRUE)
            cutProtonAngle.Fill(smalestAngleDiff, helpCoplanarity, smalestAngleDiffTaggerTime, smalestAngleDiffTaggerBin);
    }
    return kFALSE;
}

void    GCheckProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    cutProtonAngle.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    cutCoplanarity.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    cutBoth.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GCheckProton::SetCuts(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity)
{
    CutProtonAngleDiff      = maxProtonAngleDiff;
    CutProtonCoplanarity[0] = minCoplanarity;
    CutProtonCoplanarity[1] = maxCoplanarity;
}
