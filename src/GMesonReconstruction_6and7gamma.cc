#include "GMesonReconstruction_6and7gamma.h"
#include "GParticleReconstruction.h"

using namespace std;


GMesonReconstruction_6and7gamma::GMesonReconstruction_6and7gamma()    :
    width_pi0(22),
    width_eta(40),
    width_etap(60),
    TOFAddPhoton6Hits("TOFAddPhoton6Hits", "TOFAddPhoton6Hits", 300, -15, 15, 800, 0, 800),
    TOFSubPhoton6Hits("TOFSubPhoton6Hits", "TOFSubPhoton6Hits", 300, -15, 15, 800, 0, 800),
    TOFAddPhoton7Hits("TOFAddPhoton7Hits", "TOFAddPhoton7Hits", 300, -15, 15, 800, 0, 800),
    TOFSubPhoton7Hits("TOFSubPhoton7Hits", "TOFSubPhoton7Hits", 300, -15, 15, 800, 0, 800),
    TOFAddProton("TOFAdd", "TOFAdd", 300, -15, 15, 800, 0, 800),
    TOFSubProton("TOFSub", "TOFSub", 300, -15, 15, 800, 0, 800),
    ChiSqDist("ChiSqDist", "ChiSqDist", 500, 0, 1, 500, 0, 1),
    ChiSqDist6Hits("ChiSqDist6Hits", "ChiSqDist6Hits", 500, 0, 1, 500, 0, 1),
    ChiSqDist7Hits("ChiSqDist7Hits", "ChiSqDist7Hits", 500, 0, 1, 500, 0, 1)
{
    GHistBGSub::InitCuts(-20, 20, -530, -30);
    GHistBGSub::AddRandCut(30, 530);
}

GMesonReconstruction_6and7gamma::~GMesonReconstruction_6and7gamma()
{

}

Bool_t GMesonReconstruction_6and7gamma::Start()
{
    GetNeutralPions()->CloseForInput();
    GetEtas()->CloseForInput();
    GetEtaPrimes()->CloseForInput();

    if(!TraverseEntries(0, GetPhotons()->GetNEntries()))		return kFALSE;

    return kTRUE;
}

Bool_t	GMesonReconstruction_6and7gamma::Init()
{

    string  config = ReadConfig("Do-Meson-Reconstruction");
    if (strcmp(config.c_str(), "nokey") == 0)
    {
        meson_theta_min = 0.0;
        meson_theta_max = 180.0;
    }
    else if(sscanf( config.c_str(), "%*d %lf %lf\n", &meson_theta_min, &meson_theta_max) == 2)
    {
        cout << "meson reconstruction is active over theta range [" <<
        meson_theta_min << "," << meson_theta_max <<"]" << endl;
    }
    else
    {
        meson_theta_min = 0.0;
        meson_theta_max = 180.0;
    }

    config = ReadConfig("Cut-IM-Width-Pi0");
    sscanf( config.c_str(), "%lf\n", &width_pi0);
    if(width_pi0) cout << "Pi0 IM width cut set to " << width_pi0 << " MeV" << endl;
    else
    {
        width_pi0 = DEFAULT_PI0_IM_WIDTH;
        cout << "Pi0 IM width cut set to default (" << width_pi0 << " MeV)" << endl;
    }

    config = ReadConfig("Cut-IM-Width-Eta");
    sscanf( config.c_str(), "%lf\n", &width_eta);
    if(width_pi0) cout << "Eta IM width cut set to " << width_eta << " MeV" << endl;
    else
    {
        width_eta = DEFAULT_ETA_IM_WIDTH;
        cout << "Pi0 IM width cut set to default (" << width_eta << " MeV)" << endl;
    }

    config = ReadConfig("Cut-IM-Width-Eta-Prime");
    sscanf( config.c_str(), "%lf\n", &width_etap);
    if(width_etap) cout << "Eta-Prime IM width cut set to " << width_etap << " MeV" << endl;
    else
    {
        width_etap = DEFAULT_ETAP_IM_WIDTH;
        cout << "Eta-Prime IM width cut set to default (" << width_etap << " MeV)" << endl;
    }
    cout << endl;

    /*TOFCut = new TCutG("TOFCut", 8);
    TOFCut->SetPoint(0, -10, 40);
    TOFCut->SetPoint(1, -4, 40);
    TOFCut->SetPoint(2, -1.5, 250);
    TOFCut->SetPoint(3, 0, 100);
    TOFCut->SetPoint(4, 0, 400);
    TOFCut->SetPoint(5, -1.5, 500);
    TOFCut->SetPoint(6, -5, 200);
    TOFCut->SetPoint(7, -10, 100);*/

    return kTRUE;
}

Bool_t  GMesonReconstruction_6and7gamma::ProcessEventWithoutFilling()
{
    GetNeutralPions()->Clear();
    GetEtas()->Clear();
    GetEtaPrimes()->Clear();

    ChiSq3Pi0 = 1e10;
    ChiSqEtap = 1e10;

    if(GetNReconstructed()==6)
    {
        Reconstruct6g();
        for(int t=0; t<GetTagger()->GetNTagged(); t++)
        {
            for(int i=0; i<GetTracks()->GetNTracks(); i++)
            {
                if(i==GetProtons()->GetTrackIndex(0)) continue;
                if(GetTracks()->GetTheta(i)>21) continue;
                TOFAddPhoton6Hits.Fill(GetTagger()->GetTaggedTime(t)+GetTracks()->GetTime(i), GetTracks()->GetClusterEnergy(i), GetTagger()->GetTaggedTime(t));
                TOFSubPhoton6Hits.Fill(GetTagger()->GetTaggedTime(t)-GetTracks()->GetTime(i), GetTracks()->GetClusterEnergy(i), GetTagger()->GetTaggedTime(t));
            }
        }
        //if(minDecayIndex!=3)
            ChiSqDist6Hits.Fill(TMath::Prob(ChiSqEtap, 3), TMath::Prob(ChiSq3Pi0, 3));
    }
    else if(GetNReconstructed()==7)
    {
        if(!Reconstruct7g()) return kFALSE;
        for(int t=0; t<GetTagger()->GetNTagged(); t++)
        {
            for(int i=0; i<GetTracks()->GetNTracks(); i++)
            {
                if(i==GetProtons()->GetTrackIndex(0)) continue;
                if(GetTracks()->GetTheta(i)>21) continue;
                TOFAddPhoton7Hits.Fill(GetTagger()->GetTaggedTime(t)+GetTracks()->GetTime(i), GetTracks()->GetClusterEnergy(i), GetTagger()->GetTaggedTime(t));
                TOFSubPhoton7Hits.Fill(GetTagger()->GetTaggedTime(t)-GetTracks()->GetTime(i), GetTracks()->GetClusterEnergy(i), GetTagger()->GetTaggedTime(t));
            }
            TOFAddProton.Fill(GetTagger()->GetTaggedTime(t)+GetProtons()->GetTime(0), GetProtons()->GetClusterEnergy(0), GetTagger()->GetTaggedTime(t));
            TOFSubProton.Fill(GetTagger()->GetTaggedTime(t)-GetProtons()->GetTime(0), GetProtons()->GetClusterEnergy(0), GetTagger()->GetTaggedTime(t));
        }
        //if(minDecayIndex!=3)
            ChiSqDist7Hits.Fill(TMath::Prob(ChiSqEtap, 3), TMath::Prob(ChiSq3Pi0, 3));
    }
    else
        return kFALSE;

    //std::cout << TMath::Prob(ChiSqEtap, 3) << "   " << TMath::Prob(ChiSq3Pi0, 3) << std::endl;
    //std::cout << ChiSqEtap << "   " << ChiSq3Pi0 << std::endl;
    ChiSqDist.Fill(TMath::Prob(ChiSqEtap, 3), TMath::Prob(ChiSq3Pi0, 3));

    return kTRUE;
}

void    GMesonReconstruction_6and7gamma::Reconstruct6g()
{
    TLorentzVector  meson[15][3];
    Double_t        help[2][3];
    Double_t        ChiSq[15][4];

    for(int i=0; i<15; i++)
    {
        meson[i][0] = GetPhotons()->Particle(perm6g[i][0]) + GetPhotons()->Particle(perm6g[i][1]);
        meson[i][1] = GetPhotons()->Particle(perm6g[i][2]) + GetPhotons()->Particle(perm6g[i][3]);
        meson[i][2] = GetPhotons()->Particle(perm6g[i][4]) + GetPhotons()->Particle(perm6g[i][5]);
        help[0][0]     = (MASS_ETA - meson[i][0].M())/width_eta;
        help[0][1]     = (MASS_ETA - meson[i][1].M())/width_eta;
        help[0][2]     = (MASS_ETA - meson[i][2].M())/width_eta;
        help[1][0]     = (MASS_PI0 - meson[i][0].M())/width_pi0;
        help[1][1]     = (MASS_PI0 - meson[i][1].M())/width_pi0;
        help[1][2]     = (MASS_PI0 - meson[i][2].M())/width_pi0;
        ChiSq[i][0] = (help[0][0]*help[0][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][1] = (help[1][0]*help[1][0]) + (help[0][1]*help[0][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][2] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[0][2]*help[0][2]);
        ChiSq[i][3] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
    }

    minChiSq        = ChiSq[0][0];
    minDecayIndex   = 0;
    minIndex        = 0;

    for(int i=0; i<15; i++)
    {
        for(int d=0; d<4; d++)
        {
            if(d==3)
            {
                if(ChiSq[i][d]<=ChiSq3Pi0)
                    ChiSq3Pi0   = ChiSq[i][d];
            }
            else
            {
                if(ChiSq[i][d]<=ChiSqEtap)
                    ChiSqEtap   = ChiSq[i][d];
            }
            if(ChiSq[i][d]<=minChiSq)
            {
                minChiSq        = ChiSq[i][d];
                minIndex        = i;
                minDecayIndex   = d;
            }
        }
    }

    if(minDecayIndex == 3)      //found 3Pi0
    {
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        //GetNeutralPions()->AddParticle(meson[minIndex][0], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        //GetNeutralPions()->AddParticle(meson[minIndex][1], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        //GetNeutralPions()->AddParticle(meson[minIndex][2], 2, daughter_pdg, daughter_index);*/

        reconstructedEta    = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = GetPhotons()->Particle(perm6g[minIndex][0]);
        daughter[1] = GetPhotons()->Particle(perm6g[minIndex][1]);
        daughter[2] = GetPhotons()->Particle(perm6g[minIndex][2]);
        daughter[3] = GetPhotons()->Particle(perm6g[minIndex][3]);
        daughter[4] = GetPhotons()->Particle(perm6g[minIndex][4]);
        daughter[5] = GetPhotons()->Particle(perm6g[minIndex][5]);
        GetPhotons()->Particle(0) = daughter[0];
        GetPhotons()->Particle(1) = daughter[1];
        GetPhotons()->Particle(2) = daughter[2];
        GetPhotons()->Particle(3) = daughter[3];
        GetPhotons()->Particle(4) = daughter[4];
        GetPhotons()->Particle(5) = daughter[5];
        GetEtas()->AddParticle(0, 6, 0, daughter_index, daughter);
        return;
    }

    //found Eta2Pi0

    if(minDecayIndex == 0)  //Eta is meson[i][0]
    {
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = GetPhotons()->Particle(perm6g[minIndex][0]);
        daughter[1] = GetPhotons()->Particle(perm6g[minIndex][1]);
        daughter[2] = GetPhotons()->Particle(perm6g[minIndex][2]);
        daughter[3] = GetPhotons()->Particle(perm6g[minIndex][3]);
        daughter[4] = GetPhotons()->Particle(perm6g[minIndex][4]);
        daughter[5] = GetPhotons()->Particle(perm6g[minIndex][5]);
        GetPhotons()->Particle(0) = daughter[0];
        GetPhotons()->Particle(1) = daughter[1];
        GetPhotons()->Particle(2) = daughter[2];
        GetPhotons()->Particle(3) = daughter[3];
        GetPhotons()->Particle(4) = daughter[4];
        GetPhotons()->Particle(5) = daughter[5];
    }
    else if(minDecayIndex == 1)  //Eta is meson[i][1]
    {
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = GetPhotons()->Particle(perm6g[minIndex][2]);
        daughter[1] = GetPhotons()->Particle(perm6g[minIndex][3]);
        daughter[2] = GetPhotons()->Particle(perm6g[minIndex][0]);
        daughter[3] = GetPhotons()->Particle(perm6g[minIndex][1]);
        daughter[4] = GetPhotons()->Particle(perm6g[minIndex][4]);
        daughter[5] = GetPhotons()->Particle(perm6g[minIndex][5]);
        GetPhotons()->Particle(0) = daughter[0];
        GetPhotons()->Particle(1) = daughter[1];
        GetPhotons()->Particle(2) = daughter[2];
        GetPhotons()->Particle(3) = daughter[3];
        GetPhotons()->Particle(4) = daughter[4];
        GetPhotons()->Particle(5) = daughter[5];
    }
    else if(minDecayIndex == 2)  //Eta is meson[i][2]
    {
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][2];
        daughter_index[5] = perm6g[minIndex][3];
        daughter[0] = GetPhotons()->Particle(perm6g[minIndex][4]);
        daughter[1] = GetPhotons()->Particle(perm6g[minIndex][5]);
        daughter[2] = GetPhotons()->Particle(perm6g[minIndex][0]);
        daughter[3] = GetPhotons()->Particle(perm6g[minIndex][1]);
        daughter[4] = GetPhotons()->Particle(perm6g[minIndex][2]);
        daughter[5] = GetPhotons()->Particle(perm6g[minIndex][3]);
        GetPhotons()->Particle(0) = daughter[0];
        GetPhotons()->Particle(1) = daughter[1];
        GetPhotons()->Particle(2) = daughter[2];
        GetPhotons()->Particle(3) = daughter[3];
        GetPhotons()->Particle(4) = daughter[4];
        GetPhotons()->Particle(5) = daughter[5];
    }

    reconstructedEtap   = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
    GetEtaPrimes()->AddParticle(0, 6, 0, daughter_index, daughter);
}

void    GMesonReconstruction_6and7gamma::Reconstruct6g(TLorentzVector** vec)
{
    TLorentzVector  meson[15][3];
    Double_t        help[2][3];
    Double_t        ChiSq[15][4];

    for(int i=0; i<15; i++)
    {
        meson[i][0] = *vec[perm6g[i][0]];
        meson[i][0] += *vec[perm6g[i][1]];
        meson[i][1] = *vec[perm6g[i][2]] + *vec[perm6g[i][3]];
        meson[i][2] = *vec[perm6g[i][4]] + *vec[perm6g[i][5]];
        help[0][0]     = (MASS_ETA - meson[i][0].M())/width_eta;
        help[0][1]     = (MASS_ETA - meson[i][1].M())/width_eta;
        help[0][2]     = (MASS_ETA - meson[i][2].M())/width_eta;
        help[1][0]     = (MASS_PI0 - meson[i][0].M())/width_pi0;
        help[1][1]     = (MASS_PI0 - meson[i][1].M())/width_pi0;
        help[1][2]     = (MASS_PI0 - meson[i][2].M())/width_pi0;
        ChiSq[i][0] = (help[0][0]*help[0][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][1] = (help[1][0]*help[1][0]) + (help[0][1]*help[0][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][2] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[0][2]*help[0][2]);
        ChiSq[i][3] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
    }

    minChiSq        = ChiSq[0][0];
    minDecayIndex   = 0;
    minIndex        = 0;

    for(int i=0; i<15; i++)
    {
        for(int d=0; d<4; d++)
        {
            if(d==3)
            {
                if(ChiSq[i][d]<=ChiSq3Pi0)
                    ChiSq3Pi0   = ChiSq[i][d];
            }
            else
            {
                if(ChiSq[i][d]<=ChiSqEtap)
                    ChiSqEtap   = ChiSq[i][d];
            }
            if(ChiSq[i][d]<=minChiSq)
            {
                minChiSq        = ChiSq[i][d];
                minIndex        = i;
                minDecayIndex   = d;
            }
        }
    }

    if(minDecayIndex == 3)      //found 3Pi0
    {
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        //GetNeutralPions()->AddParticle(meson[minIndex][0], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        //GetNeutralPions()->AddParticle(meson[minIndex][1], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        //GetNeutralPions()->AddParticle(meson[minIndex][2], 2, daughter_pdg, daughter_index);*/

        reconstructedEta    = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = *vec[perm6g[minIndex][0]];
        daughter[1] = *vec[perm6g[minIndex][1]];
        daughter[2] = *vec[perm6g[minIndex][2]];
        daughter[3] = *vec[perm6g[minIndex][3]];
        daughter[4] = *vec[perm6g[minIndex][4]];
        daughter[5] = *vec[perm6g[minIndex][5]];
        return;
    }

    //found Eta2Pi0
    if(minDecayIndex == 0)  //Eta is meson[i][0]
    {
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = *vec[perm6g[minIndex][0]];
        daughter[1] = *vec[perm6g[minIndex][1]];
        daughter[2] = *vec[perm6g[minIndex][2]];
        daughter[3] = *vec[perm6g[minIndex][3]];
        daughter[4] = *vec[perm6g[minIndex][4]];
        daughter[5] = *vec[perm6g[minIndex][5]];
    }
    else if(minDecayIndex == 1)  //Eta is meson[i][1]
    {
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = *vec[perm6g[minIndex][2]];
        daughter[1] = *vec[perm6g[minIndex][3]];
        daughter[2] = *vec[perm6g[minIndex][0]];
        daughter[3] = *vec[perm6g[minIndex][1]];
        daughter[4] = *vec[perm6g[minIndex][4]];
        daughter[5] = *vec[perm6g[minIndex][5]];
    }
    else if(minDecayIndex == 2)  //Eta is meson[i][2]
    {
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][2];
        daughter_index[5] = perm6g[minIndex][3];
        daughter[0] = *vec[perm6g[minIndex][4]];
        daughter[1] = *vec[perm6g[minIndex][5]];
        daughter[2] = *vec[perm6g[minIndex][0]];
        daughter[3] = *vec[perm6g[minIndex][1]];
        daughter[4] = *vec[perm6g[minIndex][2]];
        daughter[5] = *vec[perm6g[minIndex][3]];
    }

    reconstructedEtap   = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
}

bool    GMesonReconstruction_6and7gamma::Reconstruct7g()
{
    Double_t    bestChiSq;
    Double_t    bestIndex;
    TLorentzVector* vec[6];
    for(int l=0; l<6; l++)
        vec[l] = new TLorentzVector();
    bestChiSq   = 1e10;
    bestIndex   = 0;
    bool    found = false;
    for(int i=0; i<7; i++)
    {
        if(GetPhotons()->Particle(i).Theta()*TMath::RadToDeg() > 21) continue;
        found = true;
        int k   = 0;
        for(int l=0; l<7; l++)
        {
            if(l!=i)
            {
                *vec[k] = GetPhotons()->Particle(l);
                k++;
            }
        }
        Reconstruct6g(vec);
        if(minChiSq < bestChiSq)
        {
            bestChiSq   = minChiSq;
            bestIndex   = i;
        }
    }
    if(!found)  return false;

    int k           = 0;
    int trackIndex[6];
    for(int i=0; i<7; i++)
    {
        if(i!=bestIndex)
        {
            *vec[k] = GetPhotons()->Particle(i);
            trackIndex[k] = GetPhotons()->GetTrackIndex(i);
            k++;
        }
    }
    Reconstruct6g(vec);
    GetProtons()->AddParticle(GetTracks()->GetClusterEnergy(bestIndex),
                              GetTracks()->GetTheta(bestIndex),
                              GetTracks()->GetPhi(bestIndex),
                              MASS_PROTON,
                              GetTracks()->GetTime(bestIndex),
                              GetTracks()->GetClusterSize(bestIndex),
                              GetTracks()->GetCentralCrystal(bestIndex),
                              GetTracks()->GetCentralVeto(bestIndex),
                              GetTracks()->GetDetectors(bestIndex),
                              GetTracks()->GetVetoEnergy(bestIndex),
                              GetTracks()->GetMWPC0Energy(bestIndex),
                              GetTracks()->GetMWPC1Energy(bestIndex),
                              bestIndex);
    GetPhotons()->RemoveAllParticles();
    for(int i=0; i<6; i++)
    {
        GetPhotons()->AddParticle(GetTracks()->GetClusterEnergy(trackIndex[i]),
                                  GetTracks()->GetTheta(trackIndex[i]),
                                  GetTracks()->GetPhi(trackIndex[i]),
                                  MASS_PROTON,
                                  GetTracks()->GetTime(trackIndex[i]),
                                  GetTracks()->GetClusterSize(trackIndex[i]),
                                  GetTracks()->GetCentralCrystal(trackIndex[i]),
                                  GetTracks()->GetCentralVeto(trackIndex[i]),
                                  GetTracks()->GetDetectors(trackIndex[i]),
                                  GetTracks()->GetVetoEnergy(trackIndex[i]),
                                  GetTracks()->GetMWPC0Energy(trackIndex[i]),
                                  GetTracks()->GetMWPC1Energy(trackIndex[i]),
                                  trackIndex[i]);
    }

    if(minDecayIndex == 3)      //found 3Pi0
    {
        GetEtas()->AddParticle(0, 6, 0, daughter_index, daughter);
        return true;
    }

    GetEtaPrimes()->AddParticle(0, 6, 0, daughter_index, daughter);
    return true;
}
void  GMesonReconstruction_6and7gamma::ProcessEvent()
{
    if(!ProcessEventWithoutFilling())   return;

    GetEventParameters()->SetNReconstructed(GetNReconstructed());
    GetEventParameters()->Fill();

    //GetPhotons()->Fill();
    //GetProtons()->Fill();
    GetNeutralPions()->Fill();
    GetEtas()->Fill();
    GetEtaPrimes()->Fill();
    FillReadList();
}

TLorentzVector  GMesonReconstruction_6and7gamma::SetMass(const TLorentzVector& vec, const Double_t mass)
{
    Double_t P;
    Double_t E;
    Double_t T;

    T = vec.E() - vec.M();
    E = T + mass;
    P = TMath::Sqrt(E*E - mass*mass);

    return TLorentzVector(vec.Vect().Unit().x() * P, vec.Vect().Unit().y() * P, vec.Vect().Unit().z() * P, E);
}

TCutG*	GMesonReconstruction_6and7gamma::OpenCutFile(Char_t* file, Char_t* name)
{
    TCutG *cut;

    TFile cutFile(file, "READ");

    if( !cutFile.IsOpen() ) {
        cerr << "Can't open cut file: " << file << endl;
        throw false;
    }

    // Try to find a TCutG with the name we want
    // GetObject checks the type to be TCutG,
    // see http://root.cern.ch/root/html534/TDirectory.html#TDirectory:GetObject
    cutFile.GetObject(name, cut);

    if( !cut ) {
        cerr << "Could not find a TCutG with the name " << name << " in " << file << endl;
        throw false;
    }

    cutFile.Close();

    cout << "cut file " << file << 	" opened (Cut-name = " << name << ")"<< endl;
    return cut;
}

Int_t		GMesonReconstruction_6and7gamma::perm6g[15][6]=
{
    {0,1,2,3,4,5},
    {0,1,2,4,3,5},
    {0,1,2,5,4,3},

    {0,2,1,3,4,5},
    {0,2,1,4,3,5},
    {0,2,1,5,4,3},

    {0,3,2,1,4,5},
    {0,3,2,4,1,5},
    {0,3,2,5,4,1},

    {0,4,2,3,1,5},
    {0,4,2,1,3,5},
    {0,4,2,5,1,3},

    {0,5,2,3,4,1},
    {0,5,2,4,3,1},
    {0,5,2,1,4,3}
};
