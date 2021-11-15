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
///// \file electromagnetic/TestEm1/src/PhysicsList.cc
///// \brief Implementation of the PhysicsList class
////
////
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//#include "PhysicsList.hh"
//#include "PhysicsListMessenger.hh"
//
//#include "PhysListEmStandard.hh"
//
//#include "G4EmStandardPhysics.hh"
//#include "G4EmStandardPhysics_option1.hh"
//#include "G4EmStandardPhysics_option2.hh"
//#include "G4EmStandardPhysics_option3.hh"
//#include "G4EmStandardPhysics_option4.hh"
//#include "G4EmStandardPhysicsSS.hh"
//#include "G4EmStandardPhysicsGS.hh"
//#include "G4EmStandardPhysicsWVI.hh"
//#include "G4EmLivermorePhysics.hh"
//#include "G4EmPenelopePhysics.hh"
//#include "G4EmLowEPPhysics.hh"
//
//#include "DetectorConstruction.hh"
//
//#include "G4LossTableManager.hh"
//#include "G4UnitsTable.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4EmParameters.hh"
//
//// particles
//
//#include "G4BosonConstructor.hh"
//#include "G4LeptonConstructor.hh"
//#include "G4MesonConstructor.hh"
//#include "G4BosonConstructor.hh"
//#include "G4BaryonConstructor.hh"
//#include "G4IonConstructor.hh"
//#include "G4ShortLivedConstructor.hh"
//
//#include "G4PhysicsListHelper.hh"
//#include "G4Decay.hh"
//#include "G4RadioactiveDecayBase.hh"
//#include "G4GenericIon.hh"
//#include "G4NuclideTable.hh"
//
//#include "G4ProcessManager.hh"
//#include "StepMax.hh"
//#include "G4Material.hh"
//
//#include "G4Gamma.hh"
//#include "G4Electron.hh"
//#include "G4Proton.hh"
//#include "G4GenericIon.hh"
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//PhysicsList::PhysicsList()
//  : G4VModularPhysicsList(), fEmPhysicsList(nullptr), fEmName(" ")
//{
//  fMessenger = new PhysicsListMessenger(this);
//  SetVerboseLevel(1);
//
//  // EM physics
//  AddPhysicsList("emstandard_opt3");
//
//  // fix lower limit for cut
//  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(10*eV, 1*GeV);
//  SetDefaultCutValue(1*mm);
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//PhysicsList::~PhysicsList()
//{
//  delete fMessenger;
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::ConstructParticle()
//{
//    G4BosonConstructor  pBosonConstructor;
//    pBosonConstructor.ConstructParticle();
//
//    G4LeptonConstructor pLeptonConstructor;
//    pLeptonConstructor.ConstructParticle();
//
//    G4MesonConstructor pMesonConstructor;
//    pMesonConstructor.ConstructParticle();
//
//    G4BaryonConstructor pBaryonConstructor;
//    pBaryonConstructor.ConstructParticle();
//
//    G4IonConstructor pIonConstructor;
//    pIonConstructor.ConstructParticle();
//
//    G4ShortLivedConstructor pShortLivedConstructor;
//    pShortLivedConstructor.ConstructParticle();
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::ConstructProcess()
//{
//  // Transportation
//  //
//  AddTransportation();
//
//  // Electromagnetic physics list
//  //
//  fEmPhysicsList->ConstructProcess();
//
//  // Decay Process
//  //
//  AddDecay();
//
//  // Decay Process
//  //
//  AddRadioactiveDecay();
//
//  // step limitation (as a full process)
//  //
//  AddStepMax();
//
//  // example of Get process
//  auto process = GetProcess("RadioactiveDecay");
//  if (process != nullptr) {
//    G4cout << "\n  GetProcess : " << process->GetProcessName() << G4endl;
//  }
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::AddPhysicsList(const G4String& name)
//{
//  if (verboseLevel>0) {
//    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
//  }
//
//  if (name == fEmName) return;
//
//  if (name == "local") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new PhysListEmStandard(name);
//
//  } else if (name == "emstandard_opt0") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysics();
//
//  } else if (name == "emstandard_opt1") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysics_option1();
//
//  } else if (name == "emstandard_opt2") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysics_option2();
//
//  } else if (name == "emstandard_opt3") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysics_option3();
//
//  } else if (name == "emstandard_opt4") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysics_option4();
//
//  } else if (name == "emstandardSS") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysicsSS();
//
//  } else if (name == "emstandardGS") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysicsGS();
//
//  } else if (name == "emstandardWVI") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmStandardPhysicsWVI();
//
//  } else if (name == "emlivermore") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmLivermorePhysics();
//
//  } else if (name == "empenelope") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmPenelopePhysics();
//
//  } else if (name == "emlowenergy") {
//    fEmName = name;
//    delete fEmPhysicsList;
//    fEmPhysicsList = new G4EmLowEPPhysics();
//
//  } else {
//
//    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
//           << " is not defined"
//           << G4endl;
//  }
//
//  // Em options
//  //
//  G4EmParameters::Instance()->SetBuildCSDARange(true);
//  G4EmParameters::Instance()->SetGeneralProcessActive(false);
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::AddDecay()
//{
//  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
//
//  // Decay Process
//  //
//  G4Decay* fDecayProcess = new G4Decay();
//
//  auto particleIterator=GetParticleIterator();
//  particleIterator->reset();
//  while( (*particleIterator)() ){
//    G4ParticleDefinition* particle = particleIterator->value();
//    if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived())
//      ph->RegisterProcess(fDecayProcess, particle);
//  }
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::AddRadioactiveDecay()
//{
//  G4RadioactiveDecayBase* radioactiveDecay = new G4RadioactiveDecayBase();
//
//  radioactiveDecay->SetARM(true);                //Atomic Rearangement
//
//  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
//  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
//
//  // mandatory for G4NuclideTable
//  //
//  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::AddStepMax()
//{
//  // Step limitation seen as a process
//  StepMax* stepMaxProcess = new StepMax();
//
//  auto particleIterator=GetParticleIterator();
//  particleIterator->reset();
//  while ((*particleIterator)()){
//    G4ParticleDefinition* particle = particleIterator->value();
//    G4ProcessManager* pmanager = particle->GetProcessManager();
//
//    if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
//      pmanager->AddDiscreteProcess(stepMaxProcess);
//  }
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//G4VProcess* PhysicsList::GetProcess(const G4String& processName) const
//{
//  G4ParticleDefinition* particle = G4GenericIon::GenericIon();
//  G4ProcessVector* procList = particle->GetProcessManager()->GetProcessList();
//  G4int nbProc = particle->GetProcessManager()->GetProcessListLength();
//  for (G4int k=0; k<nbProc; k++) {
//    G4VProcess* process = (*procList)[k];
//    if (process->GetProcessName() == processName) return process;
//  }
//  return nullptr;
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
///// \file PhysicsList.cc
///// \brief Implementation of the PhysicsList class
////
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//#include "PhysicsList.hh"
//
//#include "G4SystemOfUnits.hh"
//#include "G4UnitsTable.hh"
//
//#include "HadronElasticPhysicsHP.hh"
//
//#include "G4HadronPhysicsFTFP_BERT_HP.hh"
//#include "G4HadronPhysicsQGSP_BIC.hh"
//#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
//#include "G4HadronInelasticQBBC.hh"
//#include "G4HadronPhysicsINCLXX.hh"
//
//#include "G4IonElasticPhysics.hh"
//#include "G4IonPhysicsXS.hh"
//#include "G4IonPhysicsPHP.hh"
//#include "G4IonINCLXXPhysics.hh"
//
//#include "G4StoppingPhysics.hh"
//#include "GammaNuclearPhysics.hh"
//
//#include "ElectromagneticPhysics.hh"
//#include "G4EmStandardPhysics.hh"
//#include "G4DecayPhysics.hh"
//#include "G4RadioactiveDecayPhysics.hh"
//
//#include "G4Neutron.hh"
//#include "G4ProcessManager.hh"
//#include "G4HadronicInteraction.hh"
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//PhysicsList::PhysicsList()
//:G4VModularPhysicsList(),
// fHadronElastic(nullptr), fHadronInelastic(nullptr),
// fIonElastic(nullptr), fIonInelastic(nullptr),
// fGammaNuclear(nullptr), fElectromagnetic(nullptr),
// fDecay(nullptr), fRadioactiveDecay(nullptr)
//{
//  G4int verb = 0;
//  SetVerboseLevel(verb);
//
//  //add new units
//  //
//  new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);
//  new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
//  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
//
//  // Hadron Elastic scattering
//  fHadronElastic = new HadronElasticPhysicsHP(verb);
//  RegisterPhysics(fHadronElastic);
//
//  // Hadron Inelastic Physics
//  ////fHadronInelastic = new G4HadronPhysicsFTFP_BERT_HP(verb);
//  fHadronInelastic = new G4HadronPhysicsQGSP_BIC(verb);
//  ////fHadronInelastic = new G4HadronPhysicsQGSP_BIC_AllHP(verb);
//  ////fHadronInelastic = new G4HadronInelasticQBBC(verb);
//  ////fHadronInelastic = new G4HadronPhysicsINCLXX(verb);
//  RegisterPhysics(fHadronInelastic);
//
//  // Ion Elastic Physics
//  fIonElastic = new G4IonElasticPhysics(verb);
////  RegisterPhysics(fIonElastic);
//
//  // Ion Inelastic Physics
//  fIonInelastic = new G4IonPhysicsXS(verb);
//  ////fIonInelastic = new G4IonPhysicsPHP(verb)
//  ////fIonInelastic = new G4IonINCLXXPhysics(verb);
////  RegisterPhysics(fIonInelastic);
//
//  // stopping Particles
//  ///RegisterPhysics( new G4StoppingPhysics(verb));
//
//  // Gamma-Nuclear Physics
//  fGammaNuclear = new GammaNuclearPhysics("gamma");
//  RegisterPhysics(fGammaNuclear);
//
//  // EM physics
//  fElectromagnetic = new ElectromagneticPhysics();
//  ////fElectromagnetic = new G4EmStandardPhysics();
//  RegisterPhysics(fElectromagnetic);
//
//  // Decay
//  fDecay = new G4DecayPhysics();
//  RegisterPhysics(fDecay);
//
//  // Radioactive decay
//  fRadioactiveDecay = new G4RadioactiveDecayPhysics();
//  RegisterPhysics(fRadioactiveDecay);
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//PhysicsList::~PhysicsList()
//{ }
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::ConstructProcess()
//{
//  // Transportation first (mandatory)
//  //
//  AddTransportation();
//
//  // Physics constructors
//  //
//  fHadronElastic->ConstructProcess();
//  fHadronInelastic->ConstructProcess();
////  fIonElastic->ConstructProcess();
////  fIonInelastic->ConstructProcess();
//  fGammaNuclear->ConstructProcess();
//  fElectromagnetic->ConstructProcess();
//  fDecay->ConstructProcess();
//  fRadioactiveDecay->ConstructProcess();
//  //  fRadioactiveDecay->SetICM(true);
//
//    // example of GetHadronicModel (due to bug in QGSP_BIC_AllHP)
//  //
//  G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();
//  G4HadronicProcess* process
//       = dynamic_cast<G4HadronicProcess*>(pManager->GetProcess("nCapture"));
//  G4HadronicInteraction* model = process->GetHadronicModel("nRadCapture");
//  if(model) model->SetMinEnergy(19.9*MeV);
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void PhysicsList::SetCuts()
//{
//  SetCutValue(0*mm, "proton");
//  SetCutValue(10*km, "e-");
//  SetCutValue(10*km, "e+");
//  SetCutValue(10*km, "gamma");
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListEmStandard.hh"
#include "PhysListEmLivermore.hh"
#include "PhysListEmPenelope.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
: G4VModularPhysicsList(),
  fEmPhysicsList(0),
  fMessenger(0)
{
  G4LossTableManager::Instance();

  fMessenger = new PhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  fEmName = G4String("penelope");
  fEmPhysicsList = new PhysListEmStandard(fEmName);

  //add new units for cross sections
  //
  new G4UnitDefinition( "mm2/g", "mm2/g","Surface/Mass", mm2/g);
  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PhysicsList::ConstructParticle()
{
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  // Transportation
  //
  AddTransportation();

  // Electromagnetic physics list
  //
  fEmPhysicsList->ConstructProcess();

  // Em options
  //
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetIntegral(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "standard") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmStandard(name);

  } else if (name == "livermore") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmLivermore(name);

  } else if (name == "penelope") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmPenelope(name);

  } else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
