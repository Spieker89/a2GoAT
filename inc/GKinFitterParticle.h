#ifndef _GKinFitterParticle_h
#define _GKinFitterParticle_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"

class GKinFitterParticle {

 public:
  GKinFitterParticle();
  GKinFitterParticle(const TLorentzVector lv, const Double_t sig_th, const Double_t sig_ph, const Double_t sig_e);
  GKinFitterParticle(const TLorentzVector *lv);

  // GKinFitterParticle(GKinFitterParticle );
  virtual ~GKinFitterParticle(){}

  void Set4Vector(const TLorentzVector lv)                                                  {flv=lv;SetAlpha();}
  void SetAlpha();//Always call this via Set4Vector
  void SetVAlpha(const TMatrixD V)                                                          {fV_Alpha=V;}
  void SetResolutions(const Double_t sig_th, const Double_t sig_ph, const Double_t sig_E)   {Polar2Cartesian(sig_th,sig_ph,sig_E);}
  TMatrixD GetVAlpha()          const                                                       {return fV_Alpha;}
  TMatrixD GetAlpha()           const                                                       {return fAlpha;}
  TMatrixD GetT();
  TLorentzVector Get4Vector()   const                                                       {return TLorentzVector(fAlpha[0][0], fAlpha[1][0], fAlpha[2][0], fAlpha[3][0]);}
  Int_t GetNVar()               const                                                       {return fNvar;}

  GKinFitterParticle Add(const GKinFitterParticle p1); //returns this+p1
  GKinFitterParticle Subtract(const GKinFitterParticle p1);//returns this-p1

  void Polar2Cartesian(const Double_t sig_th, const Double_t sig_ph, const Double_t sig_e);//Calculate the error matrix in cart. from sig_th etc.

 
 private:
  TMatrixD fAlpha;
  TMatrixD fV_Alpha;
  TMatrixD fT;  //Transformation matrix S. polar ->Cartesian
  TLorentzVector flv;
  Int_t fNvar; //Number of variables =4 for 4 vector
};

#endif
