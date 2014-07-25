#ifndef __GCheckProton_h__
#define __GCheckProton_h__


#include "GH1.h"


class GTreeTagger;
class GTreeParticle;
class GTreeMeson;


class   GCheckProton
{
private:
    GH1         protonAngeDiff;
    GH1         protonAngeDiffSmalest;
    GH1         protonCoplanarity;
    GH1         protonCoplanarityAfterAngleDiff;

    Double_t    CutProtonAngleDiff;
    Double_t    CutProtonCoplanarity[2];

protected:

public:
    GCheckProton(const char* name, const char* title, const char* dirName, const Double_t ProtonAngleDiff_max = 10, const Double_t ProtonCoplanarity_min = 160, const Double_t ProtonCoplanarity_max = 200);
    ~GCheckProton();

    Bool_t  Check(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
};





#endif
