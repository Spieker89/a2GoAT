#include "MyPhysics.h"



MyPhysics::MyPhysics()    :
    hist_raw("raw", "raw"),
    hist_SubImCut("SubImCut", "SubImCut"),
    hist_MMCut("MMCut", "MMCut"),
    fit3(kTRUE),
    fit4(kTRUE),
    fit3Beam(kTRUE),
    fit4Beam(kTRUE),
    fit4Proton(kTRUE),
    fit3BeamProton(kTRUE),
    fit4BeamProton(kTRUE),
    fit3_im("fit3_im", "fit3_im", 2000, 0, 2000, 48),
    fit3_cs("fit3_cs", "fit3_cs", 1000, 0, 1000, 48),
    fit3_cl("fit3_cl", "fit3_cl", 1000, 0, 1, 48),
    fit4_im("fit4_im", "fit4_im", 2000, 0, 2000, 48),
    fit4_cs("fit4_cs", "fit4_cs", 1000, 0, 1000, 48),
    fit4_cl("fit4_cl", "fit4_cl", 1000, 0, 1, 48),
    fit3Beam_im("fit3Beam_im", "fit3Beam_im", 2000, 0, 2000, 48),
    fit3Beam_cs("fit3Beam_cs", "fit3Beam_cs", 1000, 0, 1000, 48),
    fit3Beam_cl("fit3Beam_cl", "fit3Beam_cl", 1000, 0, 1, 48),
    fit4Beam_im("fit4Beam_im", "fit4Beam_im", 2000, 0, 2000, 48),
    fit4Beam_cs("fit4Beam_cs", "fit4Beam_cs", 1000, 0, 1000, 48),
    fit4Beam_cl("fit4Beam_cl", "fit4Beam_cl", 1000, 0, 1, 48),
    fit4Proton_im("fit4Proton_im", "fit4Proton_im", 2000, 0, 2000, 48),
    fit4Proton_cs("fit4Proton_cs", "fit4Proton_cs", 1000, 0, 1000, 48),
    fit4Proton_cl("fit4Proton_cl", "fit4Proton_cl", 1000, 0, 1, 48),
    fit3BeamProton_im("fit3BeamProton_im", "fit3BeamProton_im", 2000, 0, 2000, 48),
    fit3BeamProton_cs("fit3BeamProton_cs", "fit3BeamProton_cs", 1000, 0, 1000, 48),
    fit3BeamProton_cl("fit3BeamProton_cl", "fit3BeamProton_cl", 1000, 0, 1, 48),
    fit4BeamProton_im("fit4BeamProton_im", "fit4BeamProton_im", 2000, 0, 2000, 48),
    fit4BeamProton_cs("fit4BeamProton_cs", "fit4BeamProton_cs", 1000, 0, 1000, 48),
    fit4BeamProton_cl("fit4BeamProton_cl", "fit4BeamProton_cl", 1000, 0, 1, 48)
{ 
        GHistBGSub::InitCuts(-20, 20, -55, -35);
        GHistBGSub::AddRandCut(35, 55);
}

MyPhysics::~MyPhysics()
{

}

Bool_t	MyPhysics::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();


    TraverseValidEvents();


	return kTRUE;
}

void	MyPhysics::ProcessEvent()
{
    if(etap->GetNParticles()>0)
    {
        Double_t    im  = etap->Particle(0).M();
        Double_t    mm;
        Double_t    sub_im_0    = (etap->SubPhotons(0, 0) + etap->SubPhotons(0, 1)).M();
        Double_t    sub_im_1    = (etap->SubPhotons(0, 2) + etap->SubPhotons(0, 3)).M();
        Double_t    sub_im_2    = (etap->SubPhotons(0, 4) + etap->SubPhotons(0, 5)).M();

        for(int i=0; i<tagger->GetNTagged(); i++)
        {
            mm  = (tagger->GetVectorProtonTarget(i)-etap->Particle(0)).M();
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger->GetTagged_t(i));
        }

        if((sub_im_0>500 && sub_im_0<580) &&
            (sub_im_1>100 && sub_im_1<170) &&
            (sub_im_2>100 && sub_im_2<170))
        {

            fit3.Set(etap->SubPhotons(0, 0), etap->SubPhotons(0, 1), etap->SubPhotons(0, 2), etap->SubPhotons(0, 3), etap->SubPhotons(0, 4), etap->SubPhotons(0, 5));

            if(fit3.Solve()==kTRUE)
            {
                fit3_im.Fill(fit3.GetTotalFitParticle().M());
                fit3_cs.Fill(fit3.GetChi2());
                fit3_cl.Fill(fit3.ConfidenceLevel());
            }

            for(int i=0; i<tagger->GetNTagged(); i++)
            {
                fit4.Set(etap->SubPhotons(0, 0), etap->SubPhotons(0, 1), etap->SubPhotons(0, 2), etap->SubPhotons(0, 3), etap->SubPhotons(0, 4), etap->SubPhotons(0, 5), tagger->GetVectorProtonTarget(i));
                fit3Beam.Set(etap->SubPhotons(0, 0), etap->SubPhotons(0, 1), etap->SubPhotons(0, 2), etap->SubPhotons(0, 3), etap->SubPhotons(0, 4), etap->SubPhotons(0, 5), tagger->GetVectorProtonTarget(i));
                fit4Beam.Set(etap->SubPhotons(0, 0), etap->SubPhotons(0, 1), etap->SubPhotons(0, 2), etap->SubPhotons(0, 3), etap->SubPhotons(0, 4), etap->SubPhotons(0, 5), tagger->GetVectorProtonTarget(i));
                if(protons->GetNParticles()>0)
                {
                    fit4Proton.Set(etap->SubPhotons(0, 0), etap->SubPhotons(0, 1), etap->SubPhotons(0, 2), etap->SubPhotons(0, 3), etap->SubPhotons(0, 4), etap->SubPhotons(0, 5), tagger->GetVectorProtonTarget(i), protons->Particle(0));
                    fit3BeamProton.Set(etap->SubPhotons(0, 0), etap->SubPhotons(0, 1), etap->SubPhotons(0, 2), etap->SubPhotons(0, 3), etap->SubPhotons(0, 4), etap->SubPhotons(0, 5), tagger->GetVectorProtonTarget(i), protons->Particle(0));
                    fit4BeamProton.Set(etap->SubPhotons(0, 0), etap->SubPhotons(0, 1), etap->SubPhotons(0, 2), etap->SubPhotons(0, 3), etap->SubPhotons(0, 4), etap->SubPhotons(0, 5), tagger->GetVectorProtonTarget(i), protons->Particle(0));
                }


                if(fit4.Solve()>0)
                {
                    fit4_im.Fill(fit4.GetTotalFitParticle().M());
                    fit4_cs.Fill(fit4.GetChi2());
                    fit4_cl.Fill(fit4.ConfidenceLevel());
                }
                if(fit3Beam.Solve()>0)
                {
                    fit3Beam_im.Fill(fit3Beam.GetTotalFitParticle().M());
                    fit3Beam_cs.Fill(fit3Beam.GetChi2());
                    fit3Beam_cl.Fill(fit3Beam.ConfidenceLevel());
                }
                if(fit4Beam.Solve()>0)
                {
                    fit4Beam_im.Fill(fit4Beam.GetTotalFitParticle().M());
                    fit4Beam_cs.Fill(fit4Beam.GetChi2());
                    fit4Beam_cl.Fill(fit4Beam.ConfidenceLevel());
                }
                if(protons->GetNParticles()>0)
                {
                    if(fit4Proton.Solve()>0)
                    {
                        fit4Proton_im.Fill(fit4Proton.GetTotalFitParticle().M());
                        fit4Proton_cs.Fill(fit4Proton.GetChi2());
                        fit4Proton_cl.Fill(fit4Proton.ConfidenceLevel());
                    }
                    if(fit3BeamProton.Solve()>0)
                    {
                        fit3BeamProton_im.Fill(fit3BeamProton.GetTotalFitParticle().M());
                        fit3BeamProton_cs.Fill(fit3BeamProton.GetChi2());
                        fit3BeamProton_cl.Fill(fit3BeamProton.ConfidenceLevel());
                    }
                    if(fit4Beam.Solve()>0)
                    {
                        fit4BeamProton_im.Fill(fit4BeamProton.GetTotalFitParticle().M());
                        fit4BeamProton_cs.Fill(fit4BeamProton.GetChi2());
                        fit4BeamProton_cl.Fill(fit4BeamProton.ConfidenceLevel());
                    }
                }

                mm  = (tagger->GetVectorProtonTarget(i)-etap->Particle(0)).M();
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger->GetTagged_t(i));
                if(mm>850 && mm<1025)
                {
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger->GetTagged_t(i));
                }
            }
        }
    }
}

void	MyPhysics::ProcessScalerRead()
{
    //hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    //hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	MyPhysics::Init(const char* configfile)
{
    /*SetConfigFile(configfile);
    Double_t    buf[8];
    std::string config = ReadConfig("Cut-Eta-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Eta-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_eta.SetHistMeson(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for eta physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Eta-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_eta.SetFitMeson(buf[0], buf[1]);
            cout << "Set Cuts for eta fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }


    config = ReadConfig("Cut-Eta-Proton-MinAngle");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf \n", &buf[0]) == 1)
        {
            config = ReadConfig("Cut-Eta-Proton-Coplanarity");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[1], &buf[2]) == 2)
                {
                    hist_eta.SetCheckProton(buf[0], buf[1], buf[2]);
                    cout << "Set Cuts for proton checking in eta data: ";
                    for(int i=0; i<3; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }

    config = ReadConfig("Cut-Eta-Proton-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Eta-Proton-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_eta.SetHistMesonProton(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for eta with proton physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Eta-Proton-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_eta.SetFitMesonProton(buf[0], buf[1]);
            cout << "Set Cuts for eta proton fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }






    config = ReadConfig("Cut-Etap-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Etap-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_etap.SetHistMeson(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for etap physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Etap-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_etap.SetFitMeson(buf[0], buf[1]);
            cout << "Set Cuts for etap fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }

    config = ReadConfig("Cut-Etap-Proton-MinAngle");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf \n", &buf[0]) == 1)
        {
            config = ReadConfig("Cut-Etap-Proton-Coplanarity");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[1], &buf[2]) == 2)
                {
                    hist_etap.SetCheckProton(buf[0], buf[1], buf[2]);
                    cout << "Set Cuts for proton checking in etap data: ";
                    for(int i=0; i<3; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }

    config = ReadConfig("Cut-Etap-Proton-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Etap-Proton-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_etap.SetHistMesonProton(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for etap with proton physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Etap-Proton-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_etap.SetFitMesonProton(buf[0], buf[1]);
            cout << "Set Cuts for etap proton fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }*/

    return kTRUE;
}
