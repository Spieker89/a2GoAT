//////////////////////////////////////////////////////////////////
//GKinFitter
//////////////////////////////////////////////////////////////////

#include "GKinFitter.h"
#include "TDecompLU.h"
#include <iostream>

//ClassImp(GKinFitter);

//-----------------------------------------------------------------------------
GKinFitter::GKinFitter(const Int_t npart, const Int_t ncon){

  fNvar=4;
  fNpart=npart;
  fNpar = fNpart*fNvar;
  fNcon=ncon;
  fNparti=0;
  fNpari=0;
  fNconi=0;

  fmAlpha0.ResizeTo(fNpar,1);
  fmAlpha1.ResizeTo(fNpar,1);
  fmAlpha2.ResizeTo(fNpar,1);
  fmV_Alpha0.ResizeTo(fNpar,fNpar);
  fmV_Alpha0InvPol.ResizeTo(fNpar-fNpart,fNpar-fNpart);
  fmV_Alpha.ResizeTo(fNpar,fNpar);
  fmD.ResizeTo(fNcon,fNpar);
  fmd.ResizeTo(fNcon,1);
  fmlamda.ResizeTo(fNcon,1);
  fmV_D.ResizeTo(fNcon,fNcon);
}

//-----------------------------------------------------------------------------
void GKinFitter::ResetMatrices(){

  fmAlpha0.Zero();
  fmAlpha1.Zero();
  fmAlpha2.Zero();
  fmV_Alpha0.Zero();
  fmV_Alpha0InvPol.Zero();
  fmV_Alpha.Zero();
  fmD.Zero();
  fmd.Zero();
  fmlamda.Zero();
  fmV_D.Zero();
}


Double_t    HHHDet2(const TMatrixD& m)
{
    return m.Determinant();
}
Double_t    HHHDet(const TMatrixD& m)
{
    TMatrixD sub[6];
    Double_t det[6];
    for(int i=0; i<6; i++)
    {
        m.GetSub(4*i, (4*i)+3, 4*i, (4*i)+3, sub[i]);
        sub[i]*= 1e10;
        det[i] = sub[i].Determinant();
        std::cout << det[i] << std::endl;
        //sub[i].Print();
    }

    return 0;
}

//-----------------------------------------------------------------------------
Int_t GKinFitter::Solve(){

  //Solve according to algorithm of Paul Avery:
  //Applied Fitting Theory VI, Formulas for Kinematic Fitting
  //see www.phys.ufl.edu/~avery/fitting.html

  if(fNpart!=fNparti){
    std::cout<<"GKinFitter::Solve() Added wrong number of particles. KinFit not completed"<<std::endl;
    return -1;
  }

  TMatrixD mDT=fmD;
  mDT.T();
  TMatrixD mV_Dinv=fmD*fmV_Alpha0*mDT;
  fmV_D=mV_Dinv;
  TDecompLU lu(fmV_D);
  if(!lu.Decompose()){
    std::cout<<"GKinFitter::Solve() Cannot invert. KinFit not completed"<<std::endl;
    return -1;
  }
  fmV_D.Invert();

  //Double_t det;
  //TMatrixD Einheit = fmV_D*mV_Dinv;
  //if(fNpart == 7) Einheit.Print();

  //Derive langrian multipliers
  fmlamda=fmV_D*fmd;
  //New parameters
  fmAlpha2=fmAlpha1-fmV_Alpha0*mDT*fmlamda;
  //New Covariant matrix
  fmV_Alpha=fmV_Alpha0-fmV_Alpha0*mDT*fmV_D*fmD*fmV_Alpha0;
  //chi2
  TMatrixD mlamdaT(fmlamda);
  mlamdaT.T();
  TMatrixD mchi2(mlamdaT*fmd);
  Cchi2=mchi2[0][0];

  TMatrixD pol0;
  pol0.ResizeTo(fNpar-fNpart,1);
  for(int p=0; p<fNpart; p++)
  {
      if(fmAlpha0[(4*p)+2][0]==0)
         pol0[(3*p)+0][0] = TMath::Pi()/2;
      else
          pol0[(3*p)+0][0] = TMath::ATan(TMath::Sqrt((fmAlpha0[(4*p)][0]*fmAlpha0[(4*p)][0])+(fmAlpha0[(4*p)+1][0]*fmAlpha0[(4*p)+1][0]))/fmAlpha0[(4*p)+2][0]);
      if(fmAlpha0[(4*p)+0][0]==0)
      {
          pol0[(3*p)+1][0] = TMath::Pi()/2;
          if(fmAlpha0[(4*p)+1][0]<0)
            pol0[(3*p)+1][0] += TMath::Pi();
      }
      else
          pol0[(3*p)+1][0] = TMath::ATan(fmAlpha0[(4*p)+1][0]/fmAlpha0[(4*p)][0]);
      pol0[(3*p)+2][0] = fmAlpha0[(4*p)+3][0];
  }

  TMatrixD pol2;
  pol2.ResizeTo(fNpar-fNpart,1);
  for(int p=0; p<fNpart; p++)
  {
      if(fmAlpha2[(4*p)+2][0]==0)
         pol2[(3*p)+0][0] = TMath::Pi()/2;
      else
          pol2[(3*p)+0][0] = TMath::ATan(TMath::Sqrt((fmAlpha2[(4*p)][0]*fmAlpha2[(4*p)][0])+(fmAlpha2[(4*p)+1][0]*fmAlpha2[(4*p)+1][0]))/fmAlpha2[(4*p)+2][0]);
      if(fmAlpha2[(4*p)+0][0]==0)
      {
          pol2[(3*p)+1][0] = TMath::Pi()/2;
          if(fmAlpha2[(4*p)+1][0]<0)
            pol2[(3*p)+1][0] += TMath::Pi();
      }
      else
          pol2[(3*p)+1][0] = TMath::ATan(fmAlpha2[(4*p)+1][0]/fmAlpha2[(4*p)][0]);
      pol2[(3*p)+2][0] = fmAlpha2[(4*p)+3][0];
  }
  TMatrixD diff(pol0);
  diff-=pol2;
  TMatrixD diffT(diff);
  diffT.T();
  //fmV_Alpha0.Print();
  //TMatrixD  fmV_Alpha0Inv(fmV_Alpha0);
  /*std::cout<<fmV_Alpha0Inv.Determinant()<<std::endl;
  TDecompLU lu2(fmV_Alpha0Inv);
  if(!lu2.Decompose()){
    std::cout<<"GKinFitter::Solve() Cannot invert. KinFit not completed"<<std::endl;
    return -1;
  }*/
  //HHHDet(fmV_Alpha0Inv);
  //fmV_Alpha0Inv*=1e10;
  //std::cout<<fmV_Alpha0Inv.Determinant()<<std::endl;
  //fmV_Alpha0Inv.Invert();
  //fmV_Alpha0Inv*=1e-10;
  //fmV_Alpha0InvPol.Print();
  //diffT.Print();
  //diff.Print();
  TMatrixD mvchi2(diffT*fmV_Alpha0InvPol*diff);
  Vchi2=mvchi2[0][0];

  return 1;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddInvMassConstraint(const Double_t Minv)
{
    TLorentzVector ptot(0.0,0.0,0.0,0.0);
    for(Int_t i=0; i<fNpart; i++)
        ptot+=GetInitialParticle(i);


  // d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=ptot.M2()-Minv*Minv;

  // D matrix (derivitives of constraint eqn)
  for(Int_t i=0; i<fNpart; i++){
    //[Cons Number][Var Number]
    fmD[fNconi][0+i*fNvar]=-2*ptot.X();
    fmD[fNconi][1+i*fNvar]=-2*ptot.Y();
    fmD[fNconi][2+i*fNvar]=-2*ptot.Z();
    fmD[fNconi][3+i*fNvar]= 2*ptot.T();
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv){

  //Add invariant mass constraint to subset of particles
  //Np is number of subs particles
  //pid[] contains the particle number i.e the order when they were added

  if(Np>fNpart){
    std::cout<<"GKinFitter::AddSubInvMassConstraint too many particles!"<<std::endl;
    return;
  }

  //Add up the particle 4 vectors
  TLorentzVector ptot(0.0,0.0,0.0,0.0);
  for(Int_t i=0; i<Np; i++){
    ptot+=GetInitialParticle(pid[i]);
  }
  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=ptot.M2()-Minv*Minv;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0;i<Np;i++){
    //[Cons Number][Var Number]
    fmD[fNconi][0+pid[i]*fNvar]=-2*ptot.X();
    fmD[fNconi][1+pid[i]*fNvar]=-2*ptot.Y();
    fmD[fNconi][2+pid[i]*fNvar]=-2*ptot.Z();
    fmD[fNconi][3+pid[i]*fNvar]= 2*ptot.T();
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddTotEnergyConstraint(const Double_t Etot)
{
    TLorentzVector ptot(0.0,0.0,0.0,0.0);
    for(Int_t i=0; i<fNpart; i++)
        ptot+=GetInitialParticle(i);

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=ptot.E()-Etot;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0; i<fNpart; i++){
    //[Cons Number][Var Number]
    fmD[fNconi][0+i*fNvar]=0;
    fmD[fNconi][1+i*fNvar]=0;
    fmD[fNconi][2+i*fNvar]=0;
    fmD[fNconi][3+i*fNvar]=1;
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddTotMomentumConstraint(TVector3 mom)
{
    TLorentzVector ptot(0.0,0.0,0.0,0.0);
    for(Int_t i=0; i<fNpart; i++)
        ptot+=GetInitialParticle(i);

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi+0][0]=ptot.X()-mom.X();
  fmd[fNconi+1][0]=ptot.Y()-mom.Y();
  fmd[fNconi+2][0]=ptot.Z()-mom.Z();

  //D matrix (derivitives of constraint eqn)
  Double_t D[3][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0}};
  for(Int_t i=0; i<fNpart; i++){
    for(Int_t j=0;j<3;j++){
      //[Cons Number][Var Number]
      fmD[fNconi+j][0+i*fNvar]=D[j][0];
      fmD[fNconi+j][1+i*fNvar]=D[j][1];
      fmD[fNconi+j][2+i*fNvar]=D[j][2];
      fmD[fNconi+j][3+i*fNvar]=D[j][3];
    }
  }
  fNconi+=3;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddSubMissMassConstraint(const TLorentzVector Mom, const Int_t Np, const Int_t pid[], const Double_t MissMass){

  // Add missing mass constraint to subset of particles
  // Np is number of subset particles
  // pid[] contains the particle number i.e the order when they were added

  if(Np>fNpart){
    std::cout<<"GKinFitter::AddSubMissMassConstraint too many particles!"<<std::endl;
    return;
  }

  //Add up the particle 4 vectors
  TLorentzVector Ptot(0.0,0.0,0.0,0.0);
  for(Int_t i=0; i<Np; i++){
    Ptot += GetInitialParticle(pid[i]);
  }

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=(Mom-Ptot).M2()-MissMass*MissMass;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0 ;i<Np ;i++)
    {//[Cons Number][Var Number]
    fmD[fNconi][0+pid[i]*fNvar]=-2*(Ptot-Mom).X();
    fmD[fNconi][1+pid[i]*fNvar]=-2*(Ptot-Mom).Y();
    fmD[fNconi][2+pid[i]*fNvar]=-2*(Ptot-Mom).Z();
    fmD[fNconi][3+pid[i]*fNvar]= 2*(Ptot-Mom).T();
  }
    //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddPosKFParticle(GKinFitterParticle kfp)
{
  if(fNparti>fNpart){
    std::cout<<"GKinFitter::AddPosKFParticle already at max particles"<<std::endl;
    return;
  }

  //Add parameters to Alpha0
  fmAlpha0.SetSub(fNpari,0,kfp.GetAlpha());
  fmAlpha1.SetSub(fNpari,0,kfp.GetAlpha());
  //Add error matrix to V_Alpha0
  fmV_Alpha0.SetSub(fNpari,fNpari,kfp.GetVAlpha());
  fmV_Alpha0InvPol.SetSub(3*fNparti,3*fNparti,kfp.GetVAlphaInvPol());

  //increment counters
  fNpari+=fNvar;
  fNparti++;
}


//-----------------------------------------------------------------------------
void GKinFitter::AddNegKFParticle(GKinFitterParticle kfp)
{
  AddPosKFParticle(kfp);
  for(int i=fNpari-fNvar; i<fNpari; i++)
  {
      fmAlpha0[i][0]    = -fmAlpha0[i][0];
      fmAlpha1[i][0]    = -fmAlpha1[i][0];
  }
}

//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetTotalFitParticle()
{
  TMatrixD mtot(fNvar,1);

  //loop over the sub matrices in alpha and add to total
  for(Int_t i=0; i<fNpart; i++){
    mtot+=fmAlpha2.GetSub(i*fNvar,(i+1)*fNvar-1,0,0);
  }
  return TLorentzVector(mtot[0][0],mtot[1][0],mtot[2][0],mtot[3][0]);
}

//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetParticle(Int_t ip){

  //Return the fitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetParticle particle not in fit"<<std::endl;
    return TLorentzVector(0.0,0.0,0.0,0.0);
  }

  TMatrixD mi(fNvar,1);
  mi=fmAlpha2.GetSub(ip*fNvar,(ip+1)*fNvar-1,0,0);
  return TLorentzVector(mi[0][0],mi[1][0],mi[2][0],mi[3][0]);
}

//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetInitialParticle(Int_t ip){

  //Return the unfitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetInitialParticle particle not in fit"<<std::endl;
    return TLorentzVector(0.0,0.0,0.0,0.0);
  }

  TMatrixD mi(fNvar,1);
  mi=fmAlpha1.GetSub(ip*fNvar,(ip+1)*fNvar-1,0,0);
  return TLorentzVector(mi[0][0],mi[1][0],mi[2][0],mi[3][0]);
}

//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetOriginalParticle(Int_t ip){

  //Return the unfitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetInitialParticle particle not in fit"<<std::endl;
    return TLorentzVector(0.0,0.0,0.0,0.0);
  }

  TMatrixD mi(fNvar,1);
  mi=fmAlpha0.GetSub(ip*fNvar,(ip+1)*fNvar-1,0,0);
  return TLorentzVector(mi[0][0],mi[1][0],mi[2][0],mi[3][0]);
}

//-----------------------------------------------------------------------------
void GKinFitter::Debug(){

  std::cout<<"Alpha0 "<<std::endl;
  fmAlpha0.Print();
  std::cout<<"Alpha1 "<<std::endl;
  fmAlpha1.Print();
  std::cout<<"Alpha2 "<<std::endl;
  fmAlpha2.Print();
  std::cout<<"V_Alpha0 "<<std::endl;
  fmV_Alpha0.Print();
  std::cout<<"V_Alpha "<<std::endl;
  fmV_Alpha.Print();
  std::cout<<"d "<<std::endl;
  fmd.Print();
  std::cout<<"D "<<std::endl;
  fmD.Print();
  std::cout<<"V_D "<<std::endl;
  fmV_D.Print();
  std::cout<<"lamda "<<std::endl;
  fmlamda.Print();

  //Check D*deltaAlpha+d=0
  TMatrixD mdelAlpha=fmAlpha2-fmAlpha0;
  TMatrix mCheck1=fmD*mdelAlpha + fmd;
  std::cout<<"delAlpha"<<std::endl;
  mdelAlpha.Print();
  std::cout<<"Check1 "<<std::endl;
  mCheck1.Print();

}










GIterativeKinFitter::GIterativeKinFitter(const Int_t npart, const Int_t ncon) :
    GKinFitter(npart, ncon)
{
}

GIterativeKinFitter::~GIterativeKinFitter()
{

}

void GIterativeKinFitter::AddInvMassConstraint(const Double_t Minv)
{
    conType[fNconi] = ConstraintType_InvMass;
    var[fNconi]    = Minv;

    GKinFitter::AddInvMassConstraint(Minv);
}

void GIterativeKinFitter::AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv)
{
    conType[fNconi] = ConstraintType_SubInvMass;
    nIndices[fNconi]      = Np;
    for(int i=0; i<Np; i++)
        indices[fNconi][i]  = pid[i];
    var[fNconi]    = Minv;

    GKinFitter::AddSubInvMassConstraint(Np, pid, Minv);
}

void GIterativeKinFitter::AddTotEnergyConstraint(const Double_t Etot)
{
    conType[fNconi] = ConstraintType_InvEnergy;
    var[fNconi]  = Etot;

    GKinFitter::AddTotEnergyConstraint(Etot);
}

void GIterativeKinFitter::AddTotMomentumConstraint(const TVector3 mom)
{
    conType[fNconi] = ConstraintType_InvMomentum;
    momentum[fNconi]     = mom;

    GKinFitter::AddTotMomentumConstraint(mom);
}

void GIterativeKinFitter::AddSubMissMassConstraint(const TLorentzVector Mom, const Int_t Np, const Int_t pid[], const Double_t MissMass)
{
    conType[fNconi] = ConstraintType_MisMass;
    beam[fNconi]    = Mom;
    nIndices[fNconi]      = Np;
    for(int i=0; i<Np; i++)
        indices[fNconi][i]  = pid[i];
    var[fNconi]    = MissMass;

    GKinFitter::AddSubMissMassConstraint(Mom, Np, pid, MissMass);
}

Int_t GIterativeKinFitter::Solve()
{
    if(nIter==0)
    {
        if(GKinFitter::Solve()<0)
            return -1;
        nIter++;
        return 1;
    }
    if(nIter==10)
        return -1;

    fmAlpha1    = fmAlpha2;
    TMatrixD    V(fmV_Alpha);
    Double_t    CChiSq  = GetConstraintsChi2();
    Double_t    VChiSq  = GetVariablesChi2();
    oldChiSq = GetChi2();

    ResetConstraints();
    fmD.Zero();
    fmd.Zero();
    fmlamda.Zero();
    fmV_D.Zero();
    for(int i=0; i<fNcon; i++)
    {
        switch(conType[i])
        {
        case ConstraintType_InvMass:
            GKinFitter::AddInvMassConstraint(var[i]);
            break;
        case ConstraintType_SubInvMass:
            GKinFitter::AddSubInvMassConstraint(nIndices[i], indices[i], var[i]);
            break;
        case ConstraintType_InvEnergy:
            GKinFitter::AddTotEnergyConstraint(var[i]);
            break;
        case ConstraintType_InvMomentum:
            GKinFitter::AddTotMomentumConstraint(momentum[i]);
            i+=2;
            break;
        case ConstraintType_MisMass:
            GKinFitter::AddSubMissMassConstraint(beam[i], nIndices[i], indices[i], var[i]);
            break;
        }
    }

    if(GKinFitter::Solve()<0 || oldChiSq<GetChi2())
    {
        fmAlpha2    = fmAlpha1;
        fmV_Alpha   = V;
        Cchi2       = CChiSq;
        Vchi2       = VChiSq;
        return -1;
    }
    else
        nIter++;

    return 1;
}
