//////
////// ********************************************************************
////// * License and Disclaimer                                           *
////// *                                                                  *
////// * The  Geant4 software  is  copyright of the Copyright Holders  of *
////// * the Geant4 Collaboration.  It is provided  under  the terms  and *
////// * conditions of the Geant4 Software License,  included in the file *
////// * LICENSE and available at  http://cern.ch/geant4/license .  These *
////// * include a list of copyright holders.                             *
////// *                                                                  *
////// * Neither the authors of this software system, nor their employing *
////// * institutes,nor the agencies providing financial support for this *
////// * work  make  any representation or  warranty, express or implied, *
////// * regarding  this  software system or assume any liability for its *
////// * use.  Please see the license in the file  LICENSE  and URL above *
////// * for the full disclaimer and the limitation of liability.         *
////// *                                                                  *
////// * This  code  implementation is the result of  the  scientific and *
////// * technical work of the GEANT4 collaboration.                      *
////// * By using,  copying,  modifying or  distributing the software (or *
////// * any work based  on the software)  you  agree  to acknowledge its *
////// * use  in  resulting  scientific  publications,  and indicate your *
////// * acceptance of all terms of the Geant4 Software license.          *
////// ********************************************************************
//////
//////
//////
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
////#include "SteppingAction.hh"
////#include "Run.hh"
////#include "EventAction.hh"
////#include "HistoManager.hh"
////#include "G4ParticleTypes.hh"
////#include "G4RunManager.hh"
////#include "DetectorConstruction.hh"
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////G4double saved_eventID=0;
////G4ThreeVector P_saved;
////G4double  E_saved, mass_saved;
////SteppingAction::SteppingAction(EventAction* event)
////: G4UserSteppingAction(), fEventAction(event),fScoringVolume(0)
////{
////
////}
////
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
////SteppingAction::~SteppingAction()
////{ }
////
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
////void SteppingAction::UserSteppingAction(const G4Step* aStep)
////{
////    // count processes
////    //
////    const G4StepPoint* endPoint = aStep->GetPostStepPoint();
////    const G4VProcess* process   = endPoint->GetProcessDefinedStep();
////    Run* run = static_cast<Run*>(
////                                 G4RunManager::GetRunManager()->GetNonConstCurrentRun());
////    run->CountProcesses(process);
////
////    // energy deposit
////    //
////    if (!fScoringVolume) {
////        const DetectorConstruction* detectorConstruction
////        = static_cast<const DetectorConstruction*>
////        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
////
////        fScoringVolume = detectorConstruction->GetScoringVolume();
////    }
////    G4LogicalVolume* volume
////    = aStep->GetPreStepPoint()->GetTouchableHandle()
////    ->GetVolume()->GetLogicalVolume();
////    const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
////    G4String name   = particle->GetParticleName();
////    // check if we are in scoring volume
////    if (volume != fScoringVolume) return;
////    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
////    if( name== "e-"|| name=="e+" ){
////        if (aStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="conv"){
////            if(saved_eventID==0){
////                saved_eventID= aStep->GetTrack()->GetParentID();
////                mass_saved = particle->GetPDGMass();
////                E_saved=aStep->GetTrack()->GetKineticEnergy();
////                P_saved=aStep->GetTrack()->GetMomentumDirection();
////                //si deve salvare energia e momento della prima
////                return;
////            }
////            if (aStep->GetTrack()->GetParentID()!= saved_eventID){
////                saved_eventID=0;
////                return;
////            }
////            //si prende energia e momento della seconda E CALCOLA IM e rimette a 99999 eventid
////            G4double energy = aStep->GetTrack()->GetKineticEnergy();
////            G4double mass = particle->GetPDGMass();
////            G4ThreeVector P=aStep->GetTrack()->GetMomentumDirection();
////            G4double mom_saved=sqrt(pow(E_saved,2)-pow(mass_saved,2));
////            G4double mom=sqrt(pow(energy,2)-pow(mass,2));
////            G4double invMass=sqrt(pow(mass_saved,2)+pow(mass,2)+
////                                  2*(energy*E_saved-mom_saved*mom*P.dot(P_saved)));
////
////
////            analysisManager->FillH1(17,invMass);
////            G4double angleE = P.angle(P_saved);
////            analysisManager->FillH1(16,angleE);
////            saved_eventID=0;
////        }
////    }
////}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//
//
////
//// ********************************************************************
//// * License and Disclaimer                                           *
//// *                                                                  *
//// * The  Geant4 software  is  copyright of the Copyright Holders  of *
//// * the Geant4 Collaboration.  It is provided  under  the terms  and *
//// * conditions of the Geant4 Software License,  included in the file *
//// * LICENSE and available at  http://cern.ch/geant4/license .  These *
//// * include a list of copyright holders.                             *
//// *                                                                  *
//// * Neither the authors of this software system, nor their employing *
//// * institutes,nor the agencies providing financial support for this *
//// * work  make  any representation or  warranty, express or implied, *
//// * regarding  this  software system or assume any liability for its *
//// * use.  Please see the license in the file  LICENSE  and URL above *
//// * for the full disclaimer and the limitation of liability.         *
//// *                                                                  *
//// * This  code  implementation is the result of  the  scientific and *
//// * technical work of the GEANT4 collaboration.                      *
//// * By using,  copying,  modifying or  distributing the software (or *
//// * any work based  on the software)  you  agree  to acknowledge its *
//// * use  in  resulting  scientific  publications,  and indicate your *
//// * acceptance of all terms of the Geant4 Software license.          *
//// ********************************************************************
////
///// \file SteppingAction.cc
///// \brief Implementation of the SteppingAction class
////
////
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
////#include "SteppingAction.hh"
////#include "Run.hh"
////#include "EventAction.hh"
////#include "HistoManager.hh"
////#include "G4ParticleTypes.hh"
////#include "G4RunManager.hh"
////#include "DetectorConstruction.hh"
////#include "G4VPhysicalVolume.hh"
////#include "G4SteppingManager.hh"
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////G4double saved_eventID=999999999999;
////G4ThreeVector P_saved;
////G4double  E_saved, mass_saved;
////SteppingAction::SteppingAction(EventAction* event)
////: G4UserSteppingAction(), fEventAction(event),fScoringVolume(0)
////{
////
////}
////
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
////SteppingAction::~SteppingAction()
////{ }
////
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
////void SteppingAction::UserSteppingAction(const G4Step* aStep)
////{
////    // count processes
////    //
////    const G4StepPoint* endPoint = aStep->GetPostStepPoint();
////    const G4VProcess* process   = endPoint->GetProcessDefinedStep();
////    Run* run = static_cast<Run*>(
////                                 G4RunManager::GetRunManager()->GetNonConstCurrentRun());
////    run->CountProcesses(process);
////
////    // energy deposit
////    //
////    if (!fScoringVolume) {
////        const DetectorConstruction* detectorConstruction
////        = static_cast<const DetectorConstruction*>
////        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
////
////        fScoringVolume = detectorConstruction->GetScoringVolume();
////    }
////    G4LogicalVolume* volume
////    = aStep->GetPreStepPoint()->GetTouchableHandle()
////    ->GetVolume()->GetLogicalVolume();
////    const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
////    G4String name   = particle->GetParticleName();
////    // check if we are in scoring volume
//////    if (volume != fScoringVolume) return;
//////    // check if it is a primary particle
//////    if(aStep->GetTrack()->GetParentID()!=0)return;
//////    fEventAction->SetIn(true);
//////
//////    G4String volName;
//////    //con post vedo le particelle che escono dal volume con pre quelle che entrano
//////    if(aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
//////        G4double energy = aStep->GetTrack()->GetKineticEnergy();
//////        fEventAction->AddEflow(energy);
//////                    run->ParticleFlux(name,energy);
//////    }
////    if(aStep->GetTrack()->GetParentID()!=0)return;
////    G4String volName=aStep->GetTrack()->GetVolume()->GetName();
////    if (volume == fScoringVolume ){
////        fEventAction->SetIn(true);
////        return;}
////
////    if (aStep->GetTrack()->GetVolume()->GetName()=="World"&&fEventAction->GetIn()){
//////        if( aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
////        G4double energy = aStep->GetTrack()->GetKineticEnergy();
////                fEventAction->AddEflow(energy);
////            run->ParticleFlux(name,energy);
//////        fEventAction->SetIn(false);
////    }
////}
////
//////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
//
////
//// ********************************************************************
//// * License and Disclaimer                                           *
//// *                                                                  *
//// * The  Geant4 software  is  copyright of the Copyright Holders  of *
//// * the Geant4 Collaboration.  It is provided  under  the terms  and *
//// * conditions of the Geant4 Software License,  included in the file *
//// * LICENSE and available at  http://cern.ch/geant4/license .  These *
//// * include a list of copyright holders.                             *
//// *                                                                  *
//// * Neither the authors of this software system, nor their employing *
//// * institutes,nor the agencies providing financial support for this *
//// * work  make  any representation or  warranty, express or implied, *
//// * regarding  this  software system or assume any liability for its *
//// * use.  Please see the license in the file  LICENSE  and URL above *
//// * for the full disclaimer and the limitation of liability.         *
//// *                                                                  *
//// * This  code  implementation is the result of  the  scientific and *
//// * technical work of the GEANT4 collaboration.                      *
//// * By using,  copying,  modifying or  distributing the software (or *
//// * any work based  on the software)  you  agree  to acknowledge its *
//// * use  in  resulting  scientific  publications,  and indicate your *
//// * acceptance of all terms of the Geant4 Software license.          *
//// ********************************************************************
////
///// \file SteppingAction.cc
///// \brief Implementation of the SteppingAction class
////
////
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//#include "SteppingAction.hh"
//#include "Run.hh"
//#include "EventAction.hh"
//#include "HistoManager.hh"
//#include "G4ParticleTypes.hh"
//#include "G4RunManager.hh"
//#include "DetectorConstruction.hh"
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//G4double saved_eventID=999999999999;
//G4ThreeVector P_saved;
//G4double  E_saved, mass_saved;
//SteppingAction::SteppingAction(EventAction* event)
//: G4UserSteppingAction(), fEventAction(event),fScoringVolume(0)
//{
//
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//SteppingAction::~SteppingAction()
//{ }
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void SteppingAction::UserSteppingAction(const G4Step* aStep)
//{
//    // count processes
//    //
//    const G4StepPoint* endPoint = aStep->GetPostStepPoint();
//    const G4VProcess* process   = endPoint->GetProcessDefinedStep();
//    Run* run = static_cast<Run*>(
//                                 G4RunManager::GetRunManager()->GetNonConstCurrentRun());
//    run->CountProcesses(process);
//
//    // energy deposit
//    //
//    if (!fScoringVolume) {
//        const DetectorConstruction* detectorConstruction
//        = static_cast<const DetectorConstruction*>
//        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//
//        fScoringVolume = detectorConstruction->GetScoringVolume();
//    }
//    G4LogicalVolume* volume
//    = aStep->GetPreStepPoint()->GetTouchableHandle()
//    ->GetVolume()->GetLogicalVolume();
//    const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
//    G4String name   = particle->GetParticleName();
//    if (volume != fScoringVolume) return;
//    if(aStep->GetTrack()->GetParentID()!=0)return;//se non sono le coppie generate esci
////    G4String volName=aStep->GetTrack()->GetVolume()->GetName();
//
////    if (aStep->GetTrack()->GetVolume()->GetName()=="World"){
//        if(saved_eventID==999999999999){
//            saved_eventID= G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
//            mass_saved = particle->GetPDGMass();
//            E_saved=aStep->GetTrack()->GetKineticEnergy();
//            P_saved=aStep->GetTrack()->GetMomentumDirection();
//            //si salva energia e momento della prima particella della coppia
//            return;
//        }
//        if (G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() != saved_eventID){
//            saved_eventID=999999999999;
//            //controllo che siano le particelle di una stessa coppia
//            return;
//        }
//        //si prende energia e momento della seconda E CALCOLA MI e rimette a 99999 eventid
//        G4double energy = aStep->GetTrack()->GetKineticEnergy();
//        fEventAction->AddEflow(energy);
//        run->ParticleFlux(name,energy);
//        G4double mass = particle->GetPDGMass();
//        G4ThreeVector P=aStep->GetTrack()->GetMomentumDirection();
//        G4double mom_saved=sqrt(pow(E_saved,2)-pow(mass_saved,2));
//        G4double mom=sqrt(pow(energy,2)-pow(mass,2));
//        G4double invMass=sqrt(pow(mass_saved,2)+pow(mass,2)+
//                              2*(energy*E_saved-mom_saved*mom*P.dot(P_saved)));
//
//        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//        analysisManager->FillH1(17,invMass);
//        G4double angleE = P.angle(P_saved);
//        analysisManager->FillH1(16,angleE);
//        saved_eventID=999999999999;
//
//        G4double edepStep = aStep->GetTotalEnergyDeposit();
//        if (edepStep <= 0.) return;
//        fEventAction->AddEdep(edepStep);
//    }
////}
//
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double saved_eventID=999999999999;
G4ThreeVector P_saved;
G4double  E_saved, mass_saved;
G4double  x_saved,y_saved, z_saved;
bool buona = false;
G4double dE_saved;
SteppingAction::SteppingAction(EventAction* event)
: G4UserSteppingAction(), fEventAction(event),fScoringVolume(0)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    // count processes
    //
    const G4StepPoint* endPoint = aStep->GetPostStepPoint();
    const G4VProcess* process   = endPoint->GetProcessDefinedStep();
    //    if (!process){
    Run* run = static_cast<Run*>(
                                 G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->CountProcesses(process);
    
    // energy deposit
    //
    if (!fScoringVolume) {
        const DetectorConstruction* detectorConstruction
        = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        
        fScoringVolume = detectorConstruction->GetScoringVolume();
    }
    G4LogicalVolume* volume
    = aStep->GetPreStepPoint()->GetTouchableHandle()
    ->GetVolume()->GetLogicalVolume();
    const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
    G4String name   = particle->GetParticleName();
    // check if we are in scoring volume
    if (volume != fScoringVolume) return;
    G4double edepStep = aStep->GetTotalEnergyDeposit();
    if (edepStep <= 0.) return;
    if(particle == G4Positron::Positron())fEventAction->AddEdep(edepStep);//energia depositata nello scoringVolume dai protoni
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

