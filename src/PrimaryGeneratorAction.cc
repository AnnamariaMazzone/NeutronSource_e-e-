//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "HistoManager.hh"
#include "G4RandomDirection.hh"
#include <iostream>
#include <fstream>
G4int dim=1000000;//100000;//5000000;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4double> em(dim), ep(dim), thetam(dim), thetap(dim), phim(dim), phip(dim);
PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction(), fParticleGun(0), fDetector(det)
{
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    
    fParticleGun->SetParticleEnergy(0*eV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    std::ifstream myfile
    ("/Users/annamaria/lavoro/X17/NeutronSource_e+e-/e_0.347_V_17.0_0.00256_5m_2.txt");
    //e_0.347_EM.txt");//e_0.347_V_17.0_0.00256.txt");//_5m.txt
    
    
    //int last=0;
    int ind2=-1;
    int ind=0;
    double em_v,ep_v,thetam_v,thetap_v,phim_v,phip_v;
    if (myfile.is_open())
    {
        while (myfile >> em_v >> ep_v >> thetam_v >> thetap_v >> phim_v >> phip_v)
        {
            ind2++;
           if(ind <1000000){
                em[ind]=em_v;
                ep[ind]=ep_v;
                thetam[ind]=thetam_v;
                phim[ind]=phim_v;
                thetap[ind]=thetap_v;
                phip[ind]=phip_v;
                ind++;
            }
           else continue;
        }
        myfile.close();
    }
}
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    PrimaryGeneratorAction::~PrimaryGeneratorAction()
    {
        delete fParticleGun;
    }
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
    {
        G4int id=anEvent->GetEventID();
        G4ThreeVector Pplus, Pminus;
        G4double Eplus, Eminus, emass, pmass;
        //vertex position uniform within the absorber
        //
        G4double zmax=fDetector->GetContainerLength()+2*fDetector->GetContainerRadius()-2*fDetector->GetContainTappoThickness()-2*fDetector->GetLayerThickness();//dimensione absorber=A
        G4double uz = (zmax)*(G4UniformRand() - 0.5);//=-A/2<x<A/2
        G4double Rmax = fDetector->GetContainerRadius()-fDetector->GetContainThickness()-fDetector->GetLayerThickness();//  3*cm;//
        G4double Rmax2 = Rmax*Rmax;
        G4double R2 = G4UniformRand()*Rmax2;
        G4double R_rnd  = std::sqrt(R2);
        if(uz<=fDetector->GetContainerLength()/2. && uz>=-fDetector->GetContainerLength()/2.)
            //-d/2<x<d/2 d=LUNGHEZZA CILINDRO
        {
            G4double phi = twopi*G4UniformRand();
            G4double ux = R_rnd*std::cos(phi),
            uy = R_rnd*std::sin(phi);
            fParticleGun->SetParticlePosition(G4ThreeVector(ux,uy,uz));
        }
        else //se sono nelle calotte
        {
            G4double cosTheta = G4UniformRand();//cosTheta uniform in [0, pi/2]
            G4double halfLenghtCont=fDetector->GetContainerLength()/2.;
            if(uz<0) {
                cosTheta = -cosTheta;//cosTheta uniform in [pi/2,pi]
                halfLenghtCont= -halfLenghtCont;
            }
            //2*G4UniformRand() - 1;  //cosTheta uniform in [0, pi]
            G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
            G4double phi      = twopi*G4UniformRand();  //phi uniform in [0, 2pi]
            G4ThreeVector ur(R_rnd*sinTheta*std::cos(phi),R_rnd*sinTheta*std::sin(phi),R_rnd*cosTheta+halfLenghtCont);
            
            fParticleGun->SetParticlePosition(ur);
        }
        G4ParticleDefinition* particle
        = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        
        fParticleGun->SetParticleDefinition(particle);
        Eminus=em[id];//20*MeV;//em[id];
        fParticleGun->SetParticleEnergy(Eminus);
        
        G4double cosTheta = cos(thetam[id]*deg);//2*G4UniformRand() - 1.;
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        G4double phi = phim[id]*deg;//twopi*G4UniformRand();
        G4double vx = sinTheta*std::cos(phi),
        vy = sinTheta*std::sin(phi),
        vz = cosTheta;
        Pminus=G4ThreeVector(vx,vy,vz);
//        Pminus=G4RandomDirection();
//        Pminus=G4ThreeVector(1,0,0);
        fParticleGun->SetParticleMomentumDirection(Pminus);
        fParticleGun->GeneratePrimaryVertex(anEvent);
        
        
        G4ParticleDefinition* particle2
        = G4ParticleTable::GetParticleTable()->FindParticle("e+");
        fParticleGun->SetParticleDefinition(particle2);
        
        Eplus=ep[id];//20*MeV;
        fParticleGun->SetParticleEnergy(Eplus);
        
        cosTheta = cos(thetap[id]*deg);//2*G4UniformRand() - 1.;
        sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        phi = phip[id]*deg;//twopi*G4UniformRand();
        vx = sinTheta*std::cos(phi),
        vy = sinTheta*std::sin(phi),
        vz = cosTheta;
        Pplus=G4ThreeVector(vx,vy,vz);
//         Pplus=G4ThreeVector(-1,0,0);
        fParticleGun->SetParticleMomentumDirection(Pplus);
        fParticleGun->GeneratePrimaryVertex(anEvent);
        
              //calcolo massa invariante
               pmass=particle2->GetPDGMass();
              emass=particle->GetPDGMass();
        G4double mom_plus=sqrt(pow(Eplus,2)-pow(pmass,2));
        G4double mom_minus=sqrt(pow(Eminus,2)-pow(emass,2));
              G4double invMass=sqrt(pow(emass,2)+pow(pmass,2)+
              2*(Eplus*Eminus-mom_plus*mom_minus*Pplus.dot(Pminus)));
              G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
              analysisManager->FillH1(15,invMass);
              G4double angleE = Pplus.angle(Pminus);
              analysisManager->FillH1(14,angleE);
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



