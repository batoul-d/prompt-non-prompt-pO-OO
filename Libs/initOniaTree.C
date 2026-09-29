#ifndef initOniaTree_C
#define initOniaTree_C

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include "TObjArray.h"

#include <TClonesArray.h>


Int_t           fCurrent; //!current Tree number in a TChain

   Float_t         fMass;
   Float_t         fPt;
   Float_t         fEta;
   Float_t         fPhi;
   Int_t           fSign;

   Float_t         fTauz;
   Float_t         fTauzErr;
   Float_t         fTauxy;
   Float_t         fTauxyErr;

   Float_t         fPosX;
   Float_t         fPosY;
   Float_t         fPosZ;
   Float_t         fMCPosX;
   Float_t         fMCPosY;
   Float_t         fMCPosZ;

   Float_t         fPt1;
   Float_t         fEta1;
   Float_t         fPhi1;
   Int_t           fSign1;
   Float_t         fPt2;
   Float_t         fEta2;
   Float_t         fPhi2;
   Int_t           fSign2;
   UShort_t        fMcMask1;
   UShort_t        fMcMask2;
   Float_t         fChi2MatchMCHMID1;
   Float_t         fChi2MatchMCHMID2;
   Float_t         fChi2MatchMCHMFT1;
   Float_t         fChi2MatchMCHMFT2;
   Float_t         fPtMC1;
   Float_t         fEtaMC1;
   Float_t         fPhiMC1;
   Float_t         fEMC1;
   Float_t         fPtMC2;
   Float_t         fEtaMC2;
   Float_t         fPhiMC2;
   Float_t         fEMC2;
   Float_t         fVx1;
   Float_t         fVy1;
   Float_t         fVz1;
   Float_t         fVt1;
   Float_t         fVx2;
   Float_t         fVy2;
   Float_t         fVz2;
   Float_t         fVt2;

   UInt_t           fMcDecision;
   Int_t 	   fIsAmbig1;
   Int_t 	   fIsAmbig2;

   UShort_t 	   fNumContrib;

   Float_t	   fFwdDcaX1;
   Float_t	   fFwdDcaY1;
   Float_t	   fFwdDcaX2;
   Float_t	   fFwdDcaY2;

   Float_t	   fChi21;
   Float_t	   fChi22;

   Float_t	   fChi2pca;

   Float_t	   fSVertex;

   TBranch        *b_fMass;   //!
   TBranch        *b_fPt;   //!
   TBranch        *b_fEta;   //!
   TBranch        *b_fPhi;   //!
   TBranch        *b_fSign;   //!

   TBranch        *b_fTauz;   //!
   TBranch        *b_fTauzErr;   //!
   TBranch        *b_fTauxy;   //!
   TBranch        *b_fTauxyErr;   //!

   TBranch        *b_fPosX;   //!
   TBranch        *b_fPosY;   //!
   TBranch        *b_fPosZ;   //!
   TBranch        *b_fMCPosX;   //!
   TBranch        *b_fMCPosY;   //!
   TBranch        *b_fMCPosZ;   //!

   TBranch        *b_fPt1;   //!
   TBranch        *b_fEta1;   //!
   TBranch        *b_fPhi1;   //!
   TBranch        *b_fSign1;   //!
   TBranch        *b_fPt2;   //!
   TBranch        *b_fEta2;   //!
   TBranch        *b_fPhi2;   //!
   TBranch        *b_fSign2;   //!
   TBranch        *b_fMcMask1;   //!
   TBranch        *b_fMcMask2;   //!
   TBranch        *b_fChi2MatchMCHMID1;   //!
   TBranch        *b_fChi2MatchMCHMID2;   //!
   TBranch        *b_fChi2MatchMCHMFT1;   //!
   TBranch        *b_fChi2MatchMCHMFT2;   //!
   TBranch        *b_fPtMC1;   //!
   TBranch        *b_fEtaMC1;   //!
   TBranch        *b_fPhiMC1;   //!
   TBranch        *b_fEMC1;   //!
   TBranch        *b_fPtMC2;   //!
   TBranch        *b_fEtaMC2;   //!
   TBranch        *b_fPhiMC2;   //!
   TBranch        *b_fEMC2;   //!
   TBranch        *b_fVx1;   //!
   TBranch        *b_fVy1;   //!
   TBranch        *b_fVz1;   //!
   TBranch        *b_fVt1;   //!
   TBranch        *b_fVx2;   //!
   TBranch        *b_fVy2;   //!
   TBranch        *b_fVz2;   //!
   TBranch        *b_fVt2;   //!


   TBranch        *b_fMcDecision;   //!

   TBranch	  *b_fIsAmbig1;
   TBranch	  *b_fIsAmbig2;

   TBranch	  *b_fNumContrib;

   TBranch	  *b_fFwdDcaX1;
   TBranch	  *b_fFwdDcaX2;
   TBranch	  *b_fFwdDcaY1;
   TBranch	  *b_fFwdDcaY2;

   TBranch	  *b_fChi21;
   TBranch	  *b_fChi22;

   TBranch	  *b_fChi2pca;

   TBranch	  *b_fSVertex;


void initOniaTree(TChain *tree, string TreeName)
{
   std::cout << "[INFO] Initializing TTree " << TreeName.c_str() << std::endl;
   TChain   *fChain;   //!pointer to the analyzed TTree or TChain
   if (!tree) return;
   fChain = tree;

   fCurrent = -1;

   std::cout << "[INFO] Setting Branches " << std::endl;
 if (fChain->GetBranch("fMass"))   fChain->SetBranchAddress("fMass", &fMass, &b_fMass);
 if (fChain->GetBranch("fPt"))   fChain->SetBranchAddress("fPt", &fPt, &b_fPt);
 if (fChain->GetBranch("fEta"))   fChain->SetBranchAddress("fEta", &fEta, &b_fEta);
 if (fChain->GetBranch("fPhi"))   fChain->SetBranchAddress("fPhi", &fPhi, &b_fPhi);
 if (fChain->GetBranch("fSign"))   fChain->SetBranchAddress("fSign", &fSign, &b_fSign);
 if (fChain->GetBranch("fTauz"))  fChain->SetBranchAddress("fTauz", &fTauz, &b_fTauz);
 if (fChain->GetBranch("fTauzErr"))   fChain->SetBranchAddress("fTauzErr", &fTauzErr, &b_fTauzErr);
 if (fChain->GetBranch("fTauxy"))  fChain->SetBranchAddress("fTauxy", &fTauxy, &b_fTauxy);
 if (fChain->GetBranch("fTauxyErr"))  fChain->SetBranchAddress("fTauxyErr", &fTauxyErr, &b_fTauxyErr);

 if (fChain->GetBranch("fPosX"))   fChain->SetBranchAddress("fPosX", &fPosX, &b_fPosX);
 if (fChain->GetBranch("fPosY"))   fChain->SetBranchAddress("fPosY", &fPosY, &b_fPosY);
 if (fChain->GetBranch("fPosZ"))   fChain->SetBranchAddress("fPosZ", &fPosZ, &b_fPosZ);
 if (fChain->GetBranch("fMCPosX"))   fChain->SetBranchAddress("fMCPosX", &fMCPosX, &b_fMCPosX);
 if (fChain->GetBranch("fMCPosY"))   fChain->SetBranchAddress("fMCPosY", &fMCPosY, &b_fMCPosY);
 if (fChain->GetBranch("fMCPosZ"))   fChain->SetBranchAddress("fMCPosZ", &fMCPosZ, &b_fMCPosZ);

 if (fChain->GetBranch("fPt1"))   fChain->SetBranchAddress("fPt1", &fPt1, &b_fPt1);
 if (fChain->GetBranch("fEta1"))   fChain->SetBranchAddress("fEta1", &fEta1, &b_fEta1);
 if (fChain->GetBranch("fPhi1"))   fChain->SetBranchAddress("fPhi1", &fPhi1, &b_fPhi1);
 if (fChain->GetBranch("fSign1"))   fChain->SetBranchAddress("fSign1", &fSign1, &b_fSign1);
 if (fChain->GetBranch("fPt2"))   fChain->SetBranchAddress("fPt2", &fPt2, &b_fPt2);
 if (fChain->GetBranch("fEta2"))   fChain->SetBranchAddress("fEta2", &fEta2, &b_fEta2);
 if (fChain->GetBranch("fPhi2"))   fChain->SetBranchAddress("fPhi2", &fPhi2, &b_fPhi2);
 if (fChain->GetBranch("fSign2"))   fChain->SetBranchAddress("fSign2", &fSign2, &b_fSign2);
 if (fChain->GetBranch("fMcMask1"))   fChain->SetBranchAddress("fMcMask1", &fMcMask1, &b_fMcMask1);
 if (fChain->GetBranch("fMcMask2"))   fChain->SetBranchAddress("fMcMask2", &fMcMask2, &b_fMcMask2);
 if (fChain->GetBranch("fChi2MatchMCHMID1"))   fChain->SetBranchAddress("fChi2MatchMCHMID1", &fChi2MatchMCHMID1, &b_fChi2MatchMCHMID1);
 if (fChain->GetBranch("fChi2MatchMCHMID2"))   fChain->SetBranchAddress("fChi2MatchMCHMID2", &fChi2MatchMCHMID2, &b_fChi2MatchMCHMID2);
 if (fChain->GetBranch("fChi2MatchMCHMFT1"))   fChain->SetBranchAddress("fChi2MatchMCHMFT1", &fChi2MatchMCHMFT1, &b_fChi2MatchMCHMFT1);
 if (fChain->GetBranch("fChi2MatchMCHMFT2"))   fChain->SetBranchAddress("fChi2MatchMCHMFT2", &fChi2MatchMCHMFT2, &b_fChi2MatchMCHMFT2);
 if (fChain->GetBranch("fPtMC1"))   fChain->SetBranchAddress("fPtMC1", &fPtMC1, &b_fPtMC1);
 if (fChain->GetBranch("fEtaMC1"))   fChain->SetBranchAddress("fEtaMC1", &fEtaMC1, &b_fEtaMC1);
 if (fChain->GetBranch("fPhiMC1"))   fChain->SetBranchAddress("fPhiMC1", &fPhiMC1, &b_fPhiMC1);
 if (fChain->GetBranch("fEMC1"))   fChain->SetBranchAddress("fEMC1", &fEMC1, &b_fEMC1);
 if (fChain->GetBranch("fPtMC2"))   fChain->SetBranchAddress("fPtMC2", &fPtMC2, &b_fPtMC2);
 if (fChain->GetBranch("fEtaMC2"))   fChain->SetBranchAddress("fEtaMC2", &fEtaMC2, &b_fEtaMC2);
 if (fChain->GetBranch("fPhiMC2"))   fChain->SetBranchAddress("fPhiMC2", &fPhiMC2, &b_fPhiMC2);
 if (fChain->GetBranch("fEMC2"))   fChain->SetBranchAddress("fEMC2", &fEMC2, &b_fEMC2);
 if (fChain->GetBranch("fVx1"))   fChain->SetBranchAddress("fVx1", &fVx1, &b_fVx1);
 if (fChain->GetBranch("fVy1"))   fChain->SetBranchAddress("fVy1", &fVy1, &b_fVy1);
 if (fChain->GetBranch("fVz1"))   fChain->SetBranchAddress("fVz1", &fVz1, &b_fVz1);
 if (fChain->GetBranch("fVt1"))   fChain->SetBranchAddress("fVt1", &fVt1, &b_fVt1);
 if (fChain->GetBranch("fVx2"))   fChain->SetBranchAddress("fVx2", &fVx2, &b_fVx2);
 if (fChain->GetBranch("fVy2"))   fChain->SetBranchAddress("fVy2", &fVy2, &b_fVy2);
 if (fChain->GetBranch("fVz2"))   fChain->SetBranchAddress("fVz2", &fVz2, &b_fVz2);
 if (fChain->GetBranch("fVt2"))   fChain->SetBranchAddress("fVt2", &fVt2, &b_fVt2);

 if (fChain->GetBranch("fMcDecision"))   fChain->SetBranchAddress("fMcDecision", &fMcDecision, &b_fMcDecision);

 if (fChain->GetBranch("fIsAmbig1"))   fChain->SetBranchAddress("fIsAmbig1", &fIsAmbig1, &b_fIsAmbig1);
 if (fChain->GetBranch("fIsAmbig2"))   fChain->SetBranchAddress("fIsAmbig2", &fIsAmbig2, &b_fIsAmbig2);

 if (fChain->GetBranch("fNumContrib"))   fChain->SetBranchAddress("fNumContrib", &fNumContrib, &b_fNumContrib);

 if (fChain->GetBranch("fFwdDcaX1"))   fChain->SetBranchAddress("fFwdDcaX1", &fFwdDcaX1, &b_fFwdDcaX1);
 if (fChain->GetBranch("fFwdDcaX2"))   fChain->SetBranchAddress("fFwdDcaX2", &fFwdDcaX2, &b_fFwdDcaX2);
 if (fChain->GetBranch("fFwdDcaY1"))   fChain->SetBranchAddress("fFwdDcaY1", &fFwdDcaX1, &b_fFwdDcaY1);
 if (fChain->GetBranch("fFwdDcaY2"))   fChain->SetBranchAddress("fFwdDcaY2", &fFwdDcaX2, &b_fFwdDcaY2);

 if (fChain->GetBranch("fChi21"))   fChain->SetBranchAddress("fChi21", &fChi21, &b_fChi21);
 if (fChain->GetBranch("fChi22"))   fChain->SetBranchAddress("fChi22", &fChi22, &b_fChi22);

 if (fChain->GetBranch("fChi2pca"))   fChain->SetBranchAddress("fChi2pca", &fChi2pca, &b_fChi2pca);
 if (fChain->GetBranch("fSVertex"))   fChain->SetBranchAddress("fSVertex", &fSVertex, &b_fSVertex);

}

void iniBranch(TChain* fChain, string TreeName)
{
  cout << "[INFO] Initializing Branches of " << TreeName.c_str() << endl;
 if (fChain->GetBranch("fMass"))   fChain->GetBranch("fMass")->SetAutoDelete(false);
 if (fChain->GetBranch("fPt"))   fChain->GetBranch("fPt")->SetAutoDelete(false);
 if (fChain->GetBranch("fEta"))   fChain->GetBranch("fEta")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhi"))   fChain->GetBranch("fPhi")->SetAutoDelete(false);
 if (fChain->GetBranch("fSign"))   fChain->GetBranch("fSign")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauz"))  fChain->GetBranch("fTauz")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauzErr"))   fChain->GetBranch("fTauzErr")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauxy"))  fChain->GetBranch("fTauxy")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauxyErr"))  fChain->GetBranch("fTauxyErr")->SetAutoDelete(false);
 if (fChain->GetBranch("fPosX"))   fChain->GetBranch("fPosX")->SetAutoDelete(false);
 if (fChain->GetBranch("fPosY"))   fChain->GetBranch("fPosY")->SetAutoDelete(false);
 if (fChain->GetBranch("fPosZ"))   fChain->GetBranch("fPosZ")->SetAutoDelete(false);
 if (fChain->GetBranch("fMCPosX"))   fChain->GetBranch("fMCPosX")->SetAutoDelete(false);
 if (fChain->GetBranch("fMCPosY"))   fChain->GetBranch("fMCPosY")->SetAutoDelete(false);
 if (fChain->GetBranch("fMCPosZ"))   fChain->GetBranch("fMCPosZ")->SetAutoDelete(false);
 if (fChain->GetBranch("fPt1"))   fChain->GetBranch("fPt1")->SetAutoDelete(false);
 if (fChain->GetBranch("fEta1"))   fChain->GetBranch("fEta1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhi1"))   fChain->GetBranch("fPhi1")->SetAutoDelete(false);
 if (fChain->GetBranch("fSign1"))   fChain->GetBranch("fSign1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPt2"))   fChain->GetBranch("fPt2")->SetAutoDelete(false);
 if (fChain->GetBranch("fEta2"))   fChain->GetBranch("fEta2")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhi2"))   fChain->GetBranch("fPhi2")->SetAutoDelete(false);
 if (fChain->GetBranch("fSign2"))   fChain->GetBranch("fSign2")->SetAutoDelete(false);
 if (fChain->GetBranch("fMcMask1"))   fChain->GetBranch("fMcMask1")->SetAutoDelete(false);
 if (fChain->GetBranch("fMcMask2"))   fChain->GetBranch("fMcMask2")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMID1"))   fChain->GetBranch("fChi2MatchMCHMID1")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMID2"))   fChain->GetBranch("fChi2MatchMCHMID2")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMFT1"))   fChain->GetBranch("fChi2MatchMCHMFT1")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMFT2"))   fChain->GetBranch("fChi2MatchMCHMFT2")->SetAutoDelete(false);
 if (fChain->GetBranch("fPtMC1"))   fChain->GetBranch("fPtMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fEtaMC1"))   fChain->GetBranch("fEtaMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhiMC1"))   fChain->GetBranch("fPhiMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fEMC1"))   fChain->GetBranch("fEMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPtMC2"))   fChain->GetBranch("fPtMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fEtaMC2"))   fChain->GetBranch("fEtaMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhiMC2"))   fChain->GetBranch("fPhiMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fEMC2"))   fChain->GetBranch("fEMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVx1"))   fChain->GetBranch("fVx1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVy1"))   fChain->GetBranch("fVy1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVz1"))   fChain->GetBranch("fVz1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVt1"))   fChain->GetBranch("fVt1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVx2"))   fChain->GetBranch("fVx2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVy2"))   fChain->GetBranch("fVy2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVz2"))   fChain->GetBranch("fVz2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVt2"))   fChain->GetBranch("fVt2")->SetAutoDelete(false);
 if (fChain->GetBranch("fMcDecision"))   fChain->GetBranch("fMcDecision")->SetAutoDelete(false);
 if (fChain->GetBranch("fIsAmbig1"))   fChain->GetBranch("fIsAmbig1")->SetAutoDelete(false);
 if (fChain->GetBranch("fIsAmbig2"))   fChain->GetBranch("fIsAmbig2")->SetAutoDelete(false);
 if (fChain->GetBranch("fNumContrib"))   fChain->GetBranch("fNumContrib")->SetAutoDelete(false);
 if (fChain->GetBranch("fFwdDcaX1"))   fChain->GetBranch("fFwdDcaX1")->SetAutoDelete(false);
 if (fChain->GetBranch("fFwdDcaX2"))   fChain->GetBranch("fFwdDcaX2")->SetAutoDelete(false);
 if (fChain->GetBranch("fFwdDcaY1"))   fChain->GetBranch("fFwdDcaY1")->SetAutoDelete(false);
 if (fChain->GetBranch("fFwdDcaY2"))   fChain->GetBranch("fFwdDcaY2")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi21"))   fChain->GetBranch("fChi21")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi22"))   fChain->GetBranch("fChi22")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2pca"))   fChain->GetBranch("fChi2pca")->SetAutoDelete(false);
 if (fChain->GetBranch("fSVertex"))   fChain->GetBranch("fSVertex")->SetAutoDelete(false);

 if (fChain->GetBranch("fMass"))   fChain->SetBranchStatus("fMass", 1);
 if (fChain->GetBranch("fPt"))   fChain->SetBranchStatus("fPt", 1);
 if (fChain->GetBranch("fEta"))   fChain->SetBranchStatus("fEta", 1);
 if (fChain->GetBranch("fPhi"))   fChain->SetBranchStatus("fPhi", 1);
 if (fChain->GetBranch("fSign"))   fChain->SetBranchStatus("fSign", 1);
 if (fChain->GetBranch("fTauz"))  fChain->SetBranchStatus("fTauz", 1);
 if (fChain->GetBranch("fTauzErr"))   fChain->SetBranchStatus("fTauzErr", 1);
 if (fChain->GetBranch("fTauxy"))  fChain->SetBranchStatus("fTauxy", 1);
 if (fChain->GetBranch("fTauxyErr"))  fChain->SetBranchStatus("fTauxyErr", 1);
 if (fChain->GetBranch("fPosX"))   fChain->SetBranchStatus("fPosX", 1);
 if (fChain->GetBranch("fPosY"))   fChain->SetBranchStatus("fPosY", 1);
 if (fChain->GetBranch("fPosZ"))   fChain->SetBranchStatus("fPosZ", 1);
 if (fChain->GetBranch("fMCPosX"))   fChain->SetBranchStatus("fMCPosX", 1);
 if (fChain->GetBranch("fMCPosY"))   fChain->SetBranchStatus("fMCPosY", 1);
 if (fChain->GetBranch("fMCPosZ"))   fChain->SetBranchStatus("fMCPosZ", 1);
 if (fChain->GetBranch("fPt1"))   fChain->SetBranchStatus("fPt1", 1);
 if (fChain->GetBranch("fEta1"))   fChain->SetBranchStatus("fEta1", 1);
 if (fChain->GetBranch("fPhi1"))   fChain->SetBranchStatus("fPhi1", 1);
 if (fChain->GetBranch("fSign1"))   fChain->SetBranchStatus("fSign1", 1);
 if (fChain->GetBranch("fPt2"))   fChain->SetBranchStatus("fPt2", 1);
 if (fChain->GetBranch("fEta2"))   fChain->SetBranchStatus("fEta2", 1);
 if (fChain->GetBranch("fPhi2"))   fChain->SetBranchStatus("fPhi2", 1);
 if (fChain->GetBranch("fSign2"))   fChain->SetBranchStatus("fSign2", 1);
 if (fChain->GetBranch("fMcMask1"))   fChain->SetBranchStatus("fMcMask1", 1);
 if (fChain->GetBranch("fMcMask2"))   fChain->SetBranchStatus("fMcMask2", 1);
 if (fChain->GetBranch("fChi2MatchMCHMID1"))   fChain->SetBranchStatus("fChi2MatchMCHMID1", 1);
 if (fChain->GetBranch("fChi2MatchMCHMID2"))   fChain->SetBranchStatus("fChi2MatchMCHMID2", 1);
 if (fChain->GetBranch("fChi2MatchMCHMFT1"))   fChain->SetBranchStatus("fChi2MatchMCHMFT1", 1);
 if (fChain->GetBranch("fChi2MatchMCHMFT2"))   fChain->SetBranchStatus("fChi2MatchMCHMFT2", 1);
 if (fChain->GetBranch("fPtMC1"))   fChain->SetBranchStatus("fPtMC1", 1);
 if (fChain->GetBranch("fEtaMC1"))   fChain->SetBranchStatus("fEtaMC1", 1);
 if (fChain->GetBranch("fPhiMC1"))   fChain->SetBranchStatus("fPhiMC1", 1);
 if (fChain->GetBranch("fEMC1"))   fChain->SetBranchStatus("fEMC1", 1);
 if (fChain->GetBranch("fPtMC2"))   fChain->SetBranchStatus("fPtMC2", 1);
 if (fChain->GetBranch("fEtaMC2"))   fChain->SetBranchStatus("fEtaMC2", 1);
 if (fChain->GetBranch("fPhiMC2"))   fChain->SetBranchStatus("fPhiMC2", 1);
 if (fChain->GetBranch("fEMC2"))   fChain->SetBranchStatus("fEMC2", 1);
 if (fChain->GetBranch("fVx1"))   fChain->SetBranchStatus("fVx1", 1);
 if (fChain->GetBranch("fVy1"))   fChain->SetBranchStatus("fVy1", 1);
 if (fChain->GetBranch("fVz1"))   fChain->SetBranchStatus("fVz1", 1);
 if (fChain->GetBranch("fVt1"))   fChain->SetBranchStatus("fVt1", 1);
 if (fChain->GetBranch("fVx2"))   fChain->SetBranchStatus("fVx2", 1);
 if (fChain->GetBranch("fVy2"))   fChain->SetBranchStatus("fVy2", 1);
 if (fChain->GetBranch("fVz2"))   fChain->SetBranchStatus("fVz2", 1);
 if (fChain->GetBranch("fVt2"))   fChain->SetBranchStatus("fVt2", 1);
 if (fChain->GetBranch("fMcDecision"))   fChain->SetBranchStatus("fMcDecision", 1);
 if (fChain->GetBranch("fIsAmbig1"))   fChain->SetBranchStatus("fIsAmbig1", 1);
 if (fChain->GetBranch("fIsAmbig2"))   fChain->SetBranchStatus("fIsAmbig2", 1);
 if (fChain->GetBranch("fNumContrib"))   fChain->SetBranchStatus("fNumContrib", 1);
 if (fChain->GetBranch("fFwdDcaX1"))   fChain->SetBranchStatus("fFwdDcaX1", 1);
 if (fChain->GetBranch("fFwdDcaX2"))   fChain->SetBranchStatus("fFwdDcaX2", 1);
 if (fChain->GetBranch("fFwdDcaY1"))   fChain->SetBranchStatus("fFwdDcaY1", 1);
 if (fChain->GetBranch("fFwdDcaY2"))   fChain->SetBranchStatus("fFwdDcaY2", 1);
 if (fChain->GetBranch("fChi21"))   fChain->SetBranchStatus("fChi21", 1);
 if (fChain->GetBranch("fChi22"))   fChain->SetBranchStatus("fChi22", 1);
 if (fChain->GetBranch("fChi2pca"))   fChain->SetBranchStatus("fChi2pca", 1);
 if (fChain->GetBranch("fSVertex"))   fChain->SetBranchStatus("fSVertex", 1);
}


bool getTChain(TChain *fChain, vector<string> FileNames, string TreeName)
{
  cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;

  for (vector<string>::iterator FileName = FileNames.begin() ; FileName != FileNames.end(); ++FileName){
  TFile *inputFile = TFile::Open(FileName->c_str());
  TIter keyList(inputFile->GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)keyList())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TDirectoryFile")) continue;
    string dir = key->GetName();
    string nb = dir;
    nb.erase(0,9);
    cout << "[INFO] Adding TFile " << FileName->c_str() << dir.c_str() << endl;
    fChain->Add(Form("%s/%s/%s",FileName->c_str(),dir.c_str(),TreeName.c_str()));
  }
  }
  if (!fChain) { cout << "[ERROR] fChain was not created, some input files are missing" << endl; return false; }
  return true;
}



#endif // #ifndef initOniaTree_C
