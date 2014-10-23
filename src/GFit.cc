#include "GFit.h"
#include "GTreeMeson.h"
#include "GTreeTagger.h"


GFitPulls4Vector::GFitPulls4Vector(const char* name, const char* title) :
    Pull_Px(TString(name).Append("_Px"), TString(title).Append(" Px"), 2000, -10, 10, kFALSE),
    Pull_Py(TString(name).Append("_Py"), TString(title).Append(" Py"), 2000, -10, 10, kFALSE),
    Pull_Pz(TString(name).Append("_Pz"), TString(title).Append(" Pz"), 2000, -10, 10, kFALSE),
    Pull_E(TString(name).Append("_E"), TString(title).Append(" E"), 2000, -10, 10, kFALSE)
{

}

GFitPulls4Vector::~GFitPulls4Vector()
{

}


GFitPulls6Photons::GFitPulls6Photons(const char* name, const char* title) :
    g0(TString(name).Append("_g0"), TString(title).Append(" Photon 0")),
    g1(TString(name).Append("_g1"), TString(title).Append(" Photon 1")),
    g2(TString(name).Append("_g2"), TString(title).Append(" Photon 2")),
    g3(TString(name).Append("_g3"), TString(title).Append(" Photon 3")),
    g4(TString(name).Append("_g4"), TString(title).Append(" Photon 4")),
    g5(TString(name).Append("_g5"), TString(title).Append(" Photon 5"))
{

}

GFitPulls6Photons::~GFitPulls6Photons()
{

}








GFit::GFit(const char* name, const char* title, const Bool_t IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(IsEtap),
    fit3(6, 3, 0),
    fit4(6, 4, 0),
    fit3_ConfidenceLevel(TString(name).Append("_fit3_ConfidenceLevel"), TString(title).Append(" Fit 3 Con. ConfidenceLevel"), 1000, 0, 1, kFALSE),
    fit4_ConfidenceLevel(TString(name).Append("_fit4_ConfidenceLevel"), TString(title).Append(" Fit 4 Con. ConfidenceLevel"), 1000, 0, 1, kFALSE),
    fit3_ChiSq(TString(name).Append("_fit3_ChiSq"), TString(title).Append(" Fit 3 Con. ChiSq"), 1000, 0, 100, kFALSE),
    fit4_ChiSq(TString(name).Append("_fit4_ChiSq"), TString(title).Append(" Fit 4 Con. ChiSq"), 1000, 0, 100, kFALSE),
    fit3_Pulls(TString(name).Append("_fit3_Pulls"), TString(title).Append(" Fit 3 Con. Pull")),
    fit4_Pulls(TString(name).Append("_fit4_Pulls"), TString(title).Append(" Fit 4 Con. Pull")),
    im_fit3(TString(name).Append("_fit3"), TString(title).Append(" Fit 3 Con."), 1500, 0, 1500, kFALSE),
    im_fit4(TString(name).Append("_fit4"), TString(title).Append(" Fit 4 Con."), 1500, 0, 1500, kFALSE),
    cutConfidenceLevel3(0.1),
    cutConfidenceLevel4(0.1),
    im_fit3_cutCL(TString(name).Append("_fit3CutCL"), TString(title).Append(" Fit 3 Con. cut Con. Level"), 1500, 0, 1500, kFALSE),
    im_fit4_cutCL(TString(name).Append("_fit4CutCL"), TString(title).Append(" Fit 4 Con. cut Con. Level"), 1500, 0, 1500, kFALSE),
    fit3_Pulls_cutCL(TString(name).Append("_fit3CutCL_Pulls"), TString(title).Append(" Fit 3 Con. cut Con. Level Pulls")),
    fit4_Pulls_cutCL(TString(name).Append("_fit4CutCL_Pulls"), TString(title).Append(" Fit 4 Con. cut Con. Level Pulls"))
{
    GammaResFile   = new TFile("~/GammaRes.root");
    GammaEloss     = (TH2F*)GammaResFile->Get("Eloss");
    GammaERes      = (TH2F*)GammaResFile->Get("EResIter");
    GammaThetaRes  = (TH2F*)GammaResFile->Get("ThetaRes;1");
    GammaPhiRes    = (TH2F*)GammaResFile->Get("PhiRes;1");
}

GFit::~GFit()
{

}

void   GFit::CalcResult()
{
    fit3_ConfidenceLevel.CalcResult();
    fit4_ConfidenceLevel.CalcResult();
    fit3_ChiSq.CalcResult();
    fit4_ChiSq.CalcResult();
    fit3_Pulls.CalcResult();
    fit4_Pulls.CalcResult();
    im_fit3.CalcResult();
    im_fit4.CalcResult();
    im_fit3_cutCL.CalcResult();
    im_fit4_cutCL.CalcResult();
    fit3_Pulls_cutCL.CalcResult();
    fit4_Pulls_cutCL.CalcResult();
}

void    GFit::Fit(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    GKinFitterParticle  particle[6];

    Int_t Ebin  = 0;
    Int_t Thbin = 0;
    Float_t resth = 0;
    Float_t resph = 0;
    Float_t resE  = 0;

    fit3.Reset();
    fit4.Reset();

    for(int i=0; i<6; i++)
    {
        Ebin  = GammaEloss->GetXaxis()->FindFixBin(meson.SubPhotons(0, i).E());
        Thbin = GammaEloss->GetYaxis()->FindFixBin(meson.SubPhotons(0, i).Theta()*TMath::RadToDeg());
        // Get resolutions
        resth = GammaThetaRes->GetBinContent(Ebin, Thbin);
        resph = GammaPhiRes->GetBinContent(Ebin, Thbin);
        resE  = GammaERes->GetBinContent(Ebin, Thbin);
        if(resth==0 || resph==0 || resE==0 ) return; // If energy or angle is out of calibrated  range!
        // Now set particle parameters
        //                     LorentzVector
        particle[i].Set4Vector(meson.SubPhotons(0, i));
        //std::cout << "Res: " << resth << ", " << resph << ", " << photons->Particle(i).E()*resE << std::endl;
        particle[i].SetResolutions(resth, resph, 2 *meson.SubPhotons(0, i).E()*resE);
        fit3.AddPosKFParticle(particle[i]);
        fit4.AddPosKFParticle(particle[i]);
    }

    Fit3(meson, tagger, CreateHistogramsForTaggerBinning);
    Fit4(meson, tagger, CreateHistogramsForTaggerBinning);
}

void    GFit::Fit3(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Int_t	sub[6] = {0, 1, 2, 3, 4, 5};
    if(isEtap)
        fit3.AddSubInvMassConstraint(2, &sub[0], MASS_ETA);
    else
        fit3.AddSubInvMassConstraint(2, &sub[0], MASS_PI0);

    fit3.AddSubInvMassConstraint(2, &sub[2], MASS_PI0);
    fit3.AddSubInvMassConstraint(2, &sub[4], MASS_PI0);

    if(fit3.Solve()>=0)
    {
        fit3_ConfidenceLevel.Fill(fit3.ConfidenceLevel());
        fit3_ChiSq.Fill(fit3.GetChi2());
        fit3_Pulls.Fill(fit3);
        im_fit3.Fill(fit3.GetTotalFitParticle().Get4Vector().M(), tagger, CreateHistogramsForTaggerBinning);
        if(fit3.ConfidenceLevel()>cutConfidenceLevel3)
        {
            fit3_Pulls_cutCL.Fill(fit3);
            im_fit3_cutCL.Fill(fit3.GetTotalFitParticle().Get4Vector().M(), tagger, CreateHistogramsForTaggerBinning);
        }
    }
}

void    GFit::Fit4(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Int_t       sub[6] = {0, 1, 2, 3, 4, 5};
    Double_t    minChiSq    = 1000000000;
    Double_t    im          = 0;
    Int_t       channel     = 0;

    Bool_t  found   = kFALSE;

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        fit4.ResetConstraints();
        if(isEtap)
            fit4.AddSubInvMassConstraint(2, &sub[0], MASS_ETA);
        else
            fit4.AddSubInvMassConstraint(2, &sub[0], MASS_PI0);

        fit4.AddSubInvMassConstraint(2, &sub[2], MASS_PI0);
        fit4.AddSubInvMassConstraint(2, &sub[4], MASS_PI0);
        fit4.AddSubMissMassConstraint(tagger.GetVectorProtonTarget(i), 6, sub, MASS_PROTON);

        if(fit4.Solve()>=0)
        {
            if(fit4.GetChi2()<minChiSq)
            {
                minChiSq    = fit4.GetChi2();
                im          = fit4.GetTotalFitParticle().Get4Vector().M();
                channel     = tagger.GetTagged_ch(i);
                found       = kTRUE;
            }
        }
    }

    if(found == kTRUE)
    {
        fit4_ConfidenceLevel.Fill(fit4.ConfidenceLevel());
        fit4_ChiSq.Fill(fit4.GetChi2());
        fit4_Pulls.Fill(fit4);
        im_fit4.Fill(im, 0, channel);
        if(fit4.ConfidenceLevel()>cutConfidenceLevel4)
        {
            fit4_Pulls_cutCL.Fill(fit4);
            im_fit4_cutCL.Fill(fit4.GetTotalFitParticle().Get4Vector().M(), tagger, CreateHistogramsForTaggerBinning);
        }
    }
}

void    GFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory("Fit_3");
    fit3_ConfidenceLevel.PrepareWriteList(folder, TString(name).Append("3_ConfidenceLevel"));
    fit3_ChiSq.PrepareWriteList(folder, TString(name).Append("3_ChiSq"));
    im_fit3.PrepareWriteList(folder, TString(name).Append("3_IM"));
    GHistWriteList* subFolder  = folder->GetDirectory("Pulls");
    fit3_Pulls.PrepareWriteList(subFolder, TString(name).Append("3_Pulls"));
    subFolder  = folder->GetDirectory("Cut_ConfidenceLevel");
    im_fit3_cutCL.PrepareWriteList(subFolder, TString(name).Append("3_IM_CutConL"));
    fit3_Pulls_cutCL.PrepareWriteList(subFolder, TString(name).Append("3_Pulls_CutConL"));

    folder  = arr->GetDirectory("Fit_4");
    fit4_ConfidenceLevel.PrepareWriteList(folder, TString(name).Append("4_ConfidenceLevel"));
    fit4_ChiSq.PrepareWriteList(folder, "fit4_ChiSq");
    im_fit4.PrepareWriteList(folder, "fit4_IM");
    subFolder  = folder->GetDirectory("Pulls");
    fit4_Pulls.PrepareWriteList(subFolder, "fit4_Pulls");
    subFolder  = folder->GetDirectory("Cut_ConfidenceLevel");
    im_fit4_cutCL.PrepareWriteList(subFolder, "fit4_IM_CutConL");
    subFolder  = subFolder->GetDirectory("Pulls");
    fit4_Pulls_cutCL.PrepareWriteList(subFolder, "fit4_Pulls_CutConL");
}

void    GFit::Reset(Option_t* option)
{
    fit3_ConfidenceLevel.Reset(option);
    fit4_ConfidenceLevel.Reset(option);
    fit3_ChiSq.Reset(option);
    fit4_ChiSq.Reset(option);
    fit3_Pulls.Reset(option);
    fit4_Pulls.Reset(option);
    im_fit3.Reset(option);
    im_fit4.Reset(option);
    im_fit3_cutCL.Reset(option);
    im_fit4_cutCL.Reset(option);
    fit3_Pulls_cutCL.Reset(option);
    fit4_Pulls_cutCL.Reset(option);
}

void    GFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    fit3_ConfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_ConfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit3_ChiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fit4_ChiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_fit3_cutCL.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_fit4_cutCL.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
