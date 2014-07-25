#include "GCheckProton.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"




GCheckProton::GCheckProton(const char* name, const char* title, const char* dirName, const Double_t ProtonAngleDiff_max, const Double_t ProtonCoplanarity_min, const Double_t ProtonCoplanarity_max) :
    protonAngeDiff(TString(name).Append("_ProtonAngleDiff"), TString(title).Append(" Proton Angle Diff."), 1000, 0, 100, kTRUE, dirName),
    protonAngeDiffSmalest(TString(name).Append("_ProtonAngleDiffSmalest"), TString(title).Append(" Smalest Proton Angle Diff."), 1000, 0, 100, kTRUE, dirName),
    protonCoplanarity(TString(name).Append("_Coplanarity"), TString(title).Append(" Coplanarity"), 3600, 0, 360, kTRUE, dirName),
    protonCoplanarityAfterAngleDiff(TString(name).Append("_CoplanarityAfterAngleDiff"), TString(title).Append(" Coplanarity After Angle Diff. Cut"), 3600, 0, 360, kTRUE, dirName),
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
    Bool_t      found   = kFALSE;
    Double_t    helpAngleDiff;
    Double_t    smalestAngleDiff    = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(0)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
    Double_t    smalestAngleDiffTaggerIndex    = 0;
    Double_t    smalestAngleDiffTaggerBin      = tagger.GetTagged_ch(0);
    Double_t    helpCoplanarity     = TMath::RadToDeg()*TMath::Abs(meson.Particle(0).Phi()-proton.Particle(0).Phi());
    protonCoplanarity.Fill(helpCoplanarity);

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
        {
            helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
            if(helpAngleDiff < smalestAngleDiff)
            {
                smalestAngleDiff            = helpAngleDiff;
                smalestAngleDiffTaggerIndex = i;
                smalestAngleDiffTaggerBin   = tagger.GetTagged_ch(i);
            }
            protonAngeDiff.Fill(helpAngleDiff, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            if(helpAngleDiff<CutProtonAngleDiff)
                found = kTRUE;
        }
        else
        {
            helpAngleDiff   = TMath::RadToDeg()*(tagger.GetVectorProtonTarget(i)-meson.Particle(0)).Angle(proton.Particle(0).Vect());
            if(helpAngleDiff < smalestAngleDiff)
            {
                smalestAngleDiff            = helpAngleDiff;
                smalestAngleDiffTaggerIndex = i;
                smalestAngleDiffTaggerBin   = tagger.GetTagged_ch(i);
            }
            protonAngeDiff.Fill(helpAngleDiff, tagger.GetTagged_t(i));
            if(helpAngleDiff<CutProtonAngleDiff)
                found = kTRUE;
        }
    }

    if(found == kTRUE)
        protonCoplanarityAfterAngleDiff.Fill(helpCoplanarity);

    if(helpCoplanarity<CutProtonCoplanarity[0] || helpCoplanarity>CutProtonCoplanarity[1])
        return kFALSE;

    if(CreateHistogramsForTaggerBinning)
        protonAngeDiffSmalest.Fill(smalestAngleDiff, smalestAngleDiffTaggerIndex, smalestAngleDiffTaggerBin);
    else
        protonAngeDiffSmalest.Fill(smalestAngleDiff, smalestAngleDiffTaggerIndex);

    return found;
}

void    GCheckProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    protonAngeDiff.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    protonCoplanarity.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
