#include "TestPhysics.h"

#include "Event.h"
#include "Particle.h"
#include "ParticleType.h"
#include "utils/combinatorics.h"
#include "TH1D.h"
#include <memory>
#include <iostream>
#include "plot/HistogramFactories.h"
#include "TCanvas.h"

#include "plot/root_draw.h"
#include "Detector.h"

using namespace std;
using namespace ant;

ParticleCombinatoricsTest::ParticleCombinatoricsTest()
{
    HistogramFactory::SetName("TestPhysics");
    const BinSettings im_binning(100,0,250);
    const BinSettings energy_binning(100,0,250);
    const BinSettings npart_binning(10,0,10);

    ggim     = SmartHist<double>::makeHist("2 #gamma IM", "M_{#gamma #gamma} [MeV]","#", im_binning);
    gggim    = SmartHist<double>::makeHist("3 #gamma im","M_{#gamma #gamma #gamma} [MeV]","#", im_binning);
    nphotons = SmartHist<int>::makeHist("Number of photons", "N", "", npart_binning);
    nprotons = SmartHist<int>::makeHist("Number of protons","N","",npart_binning);

    // Build a map of ParticleType -> Histogram, and fill it
    for( auto& type : ParticleTypeDatabase() ) {
        EHists[&type] = SmartHistFactory::KinEnergy(type.PrintName()+" Energy");
    }

}


void ParticleCombinatoricsTest::ProcessEvent(const Event &event)

{
    const ParticleList& photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);
    const ParticleList& protons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Proton);
    const ParticleList& all = event.Reconstructed().Particles().GetAll();

    for( auto& particle : all ) {

        // fill the histogram corresponding to the partice type of the current particle
        auto entry = EHists.find(&particle->Type());
        if( entry != EHists.end()) {
            entry->second.Fill(particle);
        }

    }

    nphotons.Fill(photons.size());
    nprotons.Fill(protons.size());

    auto combinations2 = makeCombination(photons,2);
    do {
        TLorentzVector v;
        for( auto& i: combinations2 ) {
            v += *i;
        }

        ggim.Fill(v.M());

    } while(combinations2.next());

    auto combinations3 = makeCombination(photons,3);
    do {
        TLorentzVector v;
        for( auto& i: combinations3 ) {
            v += *i;
        }

        gggim.Fill(v.M());

    } while(combinations3.next());
}

void ParticleCombinatoricsTest::Finish()
{

}

void ParticleCombinatoricsTest::ShowResult()
{
    canvas cc("ParticleCombinatoricsTest");
    cc << ggim << gggim << nphotons << nprotons << endc;

    canvas cc2("ParticleCombinatoricsTest2");
    for( auto& e : EHists) {
        cc2 << e.second;
    }
    cc2 << endc;

}
