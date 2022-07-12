/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l analyzer.C'("delphes_output.root")'
*/

//------------------------------------------------------------------------------

R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <string.h>
#include <list>
#include <vector>

void analyzer(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  chain.Print();

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton   = treeReader->UseBranch("Photon");
  TClonesArray *branchMET      = treeReader->UseBranch("MissingET");
  TClonesArray *branchTauTag      = treeReader->UseBranch("Jet.TauTag");


  // Book histograms

  //change a
  TH1 *histJetPT = new TH1F("5jet_pt", "jet P_{T}", 100, 0.0, 250.0);
  TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
  TH1 *histDiPhotMass = new TH1F("DiPhotMass", "M_{inv}(phot_{1}, phot_{2})", 100, 40.0, 140.0); //why does this use a particles?
  TH1 *histDiJetMass = new TH1F("DiJetMass", "M_{inv}(jet_{1}, jet_{2})", 100, 40.0, 140.0);
  TH1 *histMET        = new TH1F("MET", "MET", 100, 0.0, 300.0);
  TH1 *histMuPairMass = new TH1F("MuPairmass", "M_{inv}(mu_{1}, mu_{2})", 100, 40.0, 140.0);
  TH1 *histDiJetMassZeDecay = new TH1F("DiJetMass from Z->e e Decay", "M_{inv}(jet_{1}, jet_{2})", 100, 40.0, 140.0);
  TH1 *histDiJetMassZmuDecay = new TH1F("DiJetMass from Z-> mu mu Decay", "M_{inv}(jet_{1}, jet_{2})", 100, 40.0, 140.0);
  TH1 *histDibJetMassZleptonDecay = new TH1F("DibJetMass from Z->lepton Decay", "M_{inv}(bjet_{1}, bjet_{2})", 100, 40.0, 140.0);
  TH1 *histDiTauMassZleptonDecay = new TH1F("DiTauMass from Z->lepton Decay", "M_{inv}(Tau_{1}, Tau_{2})", 100, 40.0, 140.0);

  Double_t weight = 0.0004 * 300 * 1E3 / numberOfEntries;
  Double_t yield = 0;
  Double_t y2 = 0;
  
  std::vector<Jet*> jets;
  std::vector<Jet*> taus;

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    { //why can't I see the end of this {}? 
      //jets.at(0)=branchJet->At(0);

      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      Jet* selectedTau1 = 0;
      Jet* selectedTau2 = 0;
      Jet* selectedJet1 = 0;
      Jet* selectedJet2 = 0;  

      Jet* jet1  = (Jet*) branchJet->At(0);
      Jet* jet2  = (Jet*) branchJet->At(1);

      //Jet *jet1, *jet2, *jet3, *jet4;

      // If event contains at least 1 jet
      // QUESTION: numbers not 2 jets?
      if(branchJet->GetEntries() >= 1)
      	{ 
          // Take first jet
          Jet *jet = (Jet*) branchJet->At(0);

          // Plot jet transverse momentum
          histJetPT->Fill(jet->PT);
      	  // Plot jet transverse momentum

          //histJetPT->Fill(jet1->PT);
          //why is this only supposed to look at top jet?

          //plot mass
          histDiJetMass->Fill(((jet1->P4()) + (jet2->P4())).M());

          // Calculate dijet mass
      	  // Take first two jets
      	  //jet1 *jet = (Jet*) branchJet->At(0);
          //jet2 *jet = (Jet*) branchJet->At(1);    	  

      	  // Print jet transverse momentum
      	  //cout << "Jet pt: "<<jet->PT << endl;
          
      	}

      Electron *elec1, *elec2;

      // If event contains at least 2 electrons
      if(branchElectron->GetEntries() > 1)
	{
	  // Take first two electrons
	  elec1 = (Electron *) branchElectron->At(0);
	  elec2 = (Electron *) branchElectron->At(1);
    double dielectronMass=((elec1->P4()) + (elec2->P4())).M();

	  // Plot their invariant mass
	  histMass->Fill(dielectronMass);

    if (dielectronMass>70 && dielectronMass<110)
    {//plot dijet mass for events with Z boson decay to e+/e-
      histDiJetMassZeDecay->Fill(((jet1->P4()) + (jet2->P4())).M());
      /*if(branchJet->GetEntries()>=4)
      {
        
        for (Int_t j= 0; j <= branchJet->GetEntries(); ++j)
          {
            
            Jet *jet = (Jet*) branchJet->At(j);

            if (jet->NCharged <= 3)
              {
                double deltaR = sqrt(jet->DeltaEta*jet->DeltaEta + jet->DeltaPhi*jet->DeltaPhi);
                if (deltaR < 0.25) 
                {
                  if (selectedTau1 == 0) 
                    {
                      selectedTau1 = jet;
                    }
                  else if (selectedTau2 == 0) 
                    {
                      selectedTau2 = jet;
                    }
                }
       
              }

            if (selectedJet1 == 0 && jet != selectedTau1 && jet != selectedTau2) {
              selectedJet1 = jet;
            }
            else if (selectedJet2 == 0 && jet != selectedTau1 && jet != selectedTau2) {
              selectedJet2 = jet;
            }
                
            if (selectedTau1 != 0 && selectedTau2 != 0 && selectedJet1 != 0 && selectedJet2 != 0) {
              //break; //does break stay in for loop?
              histDiTauMassZleptonDecay->Fill((selectedTau1->P4()+selectedTau2->P4()).M());
            }
          }
      }*/
      /*if (selectedTau1 == 0 || selectedTau2 == 0 || selectedJet1 == 0 || selectedJet2 == 0) 
        {
          // We do not have tau-pair and a jet-pair, so give up
          continue; //jump to beginning of loop
        }*/

        
          // We got tau pair and the jet pair
          // make your pair masses ...
      
    }

      Muon *muon1, *muon2;

      // If event contains at least 2 muons
      if(branchMuon->GetEntries() > 1)
	{
	  // Take first two muons
	  muon1 = (Muon *) branchMuon->At(0);
	  muon2 = (Muon *) branchMuon->At(1);
    double dimuonMass=((muon1->P4()) + (muon2->P4())).M();

	  // Plot their invariant mass
	  histMuPairMass->Fill(dimuonMass);

    //plot dijet mass for events with Z boson decay to e+/e-
    if (dimuonMass>70 && dimuonMass<110)
    {
     histDiJetMassZmuDecay->Fill(((jet1->P4()) + (jet2->P4())).M());
    }
    
	}

      if(branchMET->GetEntries() > 0)
	{
	  MissingET* Met = (MissingET *) branchMET->At(0);
	  histMET->Fill(Met->MET, weight);
	}
      if(branchPhoton->GetEntries() > 1)
	{
	  Photon* phot1 = (Photon *) branchPhoton->At(0);
	  Photon* phot2 = (Photon *) branchPhoton->At(1);
	  if(phot1->PT < 20 | phot2->PT < 20 | abs(phot1->Eta) > 2.5 | abs(phot2->Eta) > 2.5) continue;
	  histDiPhotMass->Fill(((phot1->P4()) + (phot2->P4())).M(), weight);
	}
    
      yield += weight;
      y2 += weight*weight;
//squared?

    }
  
  
  cout << "Event yield:   " << yield << " +/- " << sqrt(y2) << endl;
  cout << "Selection Eff: " << yield / (weight*numberOfEntries) << endl;

  // Show resulting histograms
  //histJetPT->Draw();
  //histMass->Draw();
  //TCanvas * canvas1 = new TCanvas("canvas1");
    TCanvas * canvas2 = new TCanvas("canvas2");
    //canvas1->cd();
    histJetPT->Draw();
    //canvas2->cd();
    //histMass->Draw();
    //canvas2->SaveAs("la.png");

  char outputFile[1024];
  strcpy(outputFile, inputFile);
  *strstr(outputFile, ".root") = '\0';
  strcat(outputFile, "-histos.root");
  TFile *file1 = new TFile(outputFile, "RECREATE");
  histJetPT->Write();
  histMass->Write();
  histMuPairMass->Write();
  histMET->Write();
  histDiPhotMass->Write();
  histDiJetMass->Write();
  histDiJetMassZeDecay->Write();
  histDiJetMassZmuDecay->Write();
  histDibJetMassZleptonDecay->Write();
  histDiTauMassZleptonDecay->Write();
  file1->Close();

  //distinguish b and tau jets
  //what do we call h2/mass?
  //what does source command do?
  //NCharge <= 3
  //check tau tag's N charge
  //(1) n tracks (2)eta  delR=sqrt(del eta ^2 + del phi^2) require <=0.2: plot delR; plot tau pT
  //jet->NCharged <=    dsqrt
  //questions: goals of plots
  /*
  reference jet inside loop
  */
  //
}}

