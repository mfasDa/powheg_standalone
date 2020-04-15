
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
// using namespace Pythia8;

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TClonesArray.h>
#include <TPythia8.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TFile.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh"  // for area support
#endif  // __FJCORE__


//______________________________________________________________________________
void LoadLibs()                                                                 
{                                                                               
  gSystem->Load("libfastjet");                                  
  gSystem->Load("libsiscone");                                  
  gSystem->Load("libsiscone_spherical");                        
  gSystem->Load("libfastjetplugins");                           

  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");
}      

void RunPythia8(Int_t nev  = 5000, Char_t const *foutname = "Pythia8JetSpectra_CT14nlo.root" , Int_t ndeb = 1)
{
  LoadLibs();
  clock_t begin_time = clock();

  TDatime dt;
  static UInt_t sseed = dt.Get();

  if (gSystem->Getenv("CONFIG_SEED")) {
      sseed = atoi(gSystem->Getenv("CONFIG_SEED"));
      cout<<"\nseed for Random number generation is : "<<sseed<<endl;
  }


  const Int_t   chargedOnly =   1;  // charged only or full jets

  const Int_t nR = 2;
  const Float_t Rvals[nR]  ={0.2,0.4}; // Cone radii
  const Float_t Etavals[nR]={1.4,1.0};

  TFile *fout = new TFile(foutname,"RECREATE");
  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TH1F *hFullPtSpecFull[nR]={0};
  TH1F *hFullPtSpecSubPlus[nR]={0};
  TH1F *hFullPtSpecFullEta[nR]={0};
  TH1F *hFullPtSpecSubPlusEta[nR]={0};
  TH1F *hFullPtSpecTrackcut5Gev[nR]={0};
  TH1F *hFullPtSpecTrackcut5GevSub[nR]={0};


  for (Int_t iR = 0; iR < nR; iR++) {
    hFullPtSpecFull[iR] = new TH1F(Form("Full jet cross section_R%.0f with small bin width",10*Rvals[iR]),
					  Form("Full jet cross section R=%.1f with small bin width;P_{T,Ch.jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)",Rvals[iR]),
					  200,0,200);
    hFullPtSpecFull[iR]->Sumw2();

    hFullPtSpecSubPlus[iR] = new TH1F(Form("Subtracted jet cross section_R%.0f with small bin width (plus)",10*Rvals[iR]),
					  Form("Subtracted jet cross section R=%.1f with small bin width (plus);P_{T,Ch.jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)",Rvals[iR]),
					  250,-50,200);
    hFullPtSpecSubPlus[iR]->Sumw2();

    hFullPtSpecFullEta[iR] = new TH1F(Form("Full jet cross section_R%.0f with small bin width Eta05",10*Rvals[iR]),
					  Form("Full jet cross section R=%.1f with small bin width Eta05;P_{T,Ch.jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)",Rvals[iR]),
					  200,0,200);
    hFullPtSpecFullEta[iR]->Sumw2();

    hFullPtSpecSubPlusEta[iR] = new TH1F(Form("Subtracted jet cross section_R%.0f with small bin width (plus) Eta05",10*Rvals[iR]),
					  Form("Subtracted jet cross section R=%.1f with small bin width (plus) Eta05;P_{T,Ch.jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)",Rvals[iR]),
					  250,-50,200);
    hFullPtSpecSubPlusEta[iR]->Sumw2();

    hFullPtSpecTrackcut5Gev[iR] = new TH1F(Form("hFullPtSpecTrackcut5Gev_R%.0f",10*Rvals[iR]),
					  Form("jet cross section R=%.1f with Leading Pt cut at 5Gev;P_{T,Ch.jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)",Rvals[iR]),
					  200,0,200);
    hFullPtSpecTrackcut5Gev[iR]->Sumw2();

    hFullPtSpecTrackcut5GevSub[iR] = new TH1F(Form("hSubPtSpecTrackcut5Gev_R%.0f",10*Rvals[iR]),
					  Form("jet cross section R=%.1f with Leading Pt cut at 5Gev;P_{T,Ch.jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)",Rvals[iR]),
					  200,0,200);
    hFullPtSpecTrackcut5GevSub[iR]->Sumw2();

  }

  // Array of particles
  TClonesArray *particles = new TClonesArray("TParticle", 1000);
  // Create pythia8 object
  TPythia8 *pythia8 = new TPythia8();
 

  vector<fastjet::JetDefinition> jet_def;
  Double_t ghost_maxrap = 6.0;
  fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  for (Int_t iR = 0; iR < nR; iR++) jet_def.push_back(fastjet::JetDefinition(fastjet::antikt_algorithm, Rvals[iR]));

  // Configure
  pythia8->ReadString("Next:numberShowLHA = 1");
  pythia8->ReadString("Next:numberShowInfo = 1");
  pythia8->ReadString("Next:numberShowProcess = 1");
  pythia8->ReadString("Next:numberShowEvent = 1");
  pythia8->ReadString("Main:timesAllowErrors = 10");

  pythia8->ReadString("Init:showChangedSettings = on");
  pythia8->ReadString("Init:showChangedParticleData = off");

  pythia8->ReadString("Beams:frametype = 4");

  pythia8->ReadString("Beams:LHEF = pwgevents.lhe");

  pythia8->ReadString("POWHEG:nFinal = 2");

  pythia8->ReadString("PartonLevel:MPI = on");		//! TEST

  pythia8->ReadString("111:mayDecay  = on");
  pythia8->ReadString("310:mayDecay  = off");
  pythia8->ReadString("3122:mayDecay = off");
  pythia8->ReadString("3112:mayDecay = off");
  pythia8->ReadString("3212:mayDecay = off");
  pythia8->ReadString("3222:mayDecay = off");
  pythia8->ReadString("3312:mayDecay = off");
  pythia8->ReadString("3322:mayDecay = off");
  pythia8->ReadString("3334:mayDecay = off");

  pythia8->ReadString("POWHEG:veto = 1");
  pythia8->ReadString("POWHEG:vetoCount = 1000000");
  pythia8->ReadString("POWHEG:pThard = 2");	//!
  pythia8->ReadString("POWHEG:pTemt = 0");
  pythia8->ReadString("POWHEG:emitted = 3");	//!
  pythia8->ReadString("POWHEG:pTdef = 1");
  pythia8->ReadString("POWHEG:MPIveto = 1");	//!
  pythia8->ReadString("POWHEG:QEDveto = 2");

  pythia8->ReadString("Tune:preferLHAPDF = 2");
  pythia8->ReadString("Tune:pp = 5"); 

  pythia8->ReadString("PDF:pSet = LHAPDF6:CT14nlo");

  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString(Form("Random:seed = %u", sseed%900000000));

  //
  Pythia8::Pythia *pythia = pythia8->Pythia8();

  // Add in user hooks for shower vetoing
  Pythia8::PowhegHooks *powhegHooks = NULL;
  pythia->readString("SpaceShower:pTmaxMatch = 2");
  pythia->readString("TimeShower:pTmaxMatch = 2");
  pythia->readString("MultipartonInteractions:pTmaxMatch = 2");		//!
     
  powhegHooks = new Pythia8::PowhegHooks();
  pythia->setUserHooksPtr((Pythia8::UserHooks*)powhegHooks);

  

  TString PDFused      = pythia->settings.word("PDF:pSet");
  cout<<"\n PDF used is : "<<PDFused<<endl;
  Double_t SumW(0);
 



  // Initialize
  pythia->init();
  //

  // Event loop
  for (Int_t iev = 0; iev < nev; iev++) {

    if (!(iev%1000)) {
      printf(">>>processing ev# %5d / elapsed time: ",iev);
      std::cout << float(clock() - begin_time ) /  CLOCKS_PER_SEC << endl;
      begin_time = clock();
    }

    vector<fastjet::PseudoJet> input_particles;

    pythia8->GenerateEvent();
    // pythia8->EventListing();
    pythia8->ImportParticles(particles,"All");

    Double_t evt_wght =  pythia->info.weight();
    evt_wght *= 1e-9;
    SumW += evt_wght;
    //cout<<" Event weight is : "<<evt_wght*1e9<<endl;
    hNEvent->Fill(0.5,evt_wght);


    Int_t np = particles->GetEntriesFast();

    // Particle loop
    for (Int_t ip = 0; ip < np; ip++) {
      TParticle *part = (TParticle*)particles->At(ip);
      Int_t ist = part->GetStatusCode();

      // Positive codes are final particles.
      if (ist <= 0) continue;

      Int_t pdg = part->GetPdgCode();
      Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();

      if (charge == 0.) continue;

      Float_t eta = part->Eta();
      Float_t pt = part->Pt();

      if (abs(eta) > 0.9) continue;
      if (pt < 0.15) continue;

      input_particles.push_back(fastjet::PseudoJet(part->Px(),part->Py(),part->Pz(),part->Energy())); 

    }

    if (input_particles.size() == 0) {
      //printf("No particle....\n");
      continue;
    }

    for (Int_t iR = 0; iR < nR; iR++) {

      Double_t AreaCut = 0.6 * TMath::Pi() * TMath::Power(Rvals[iR],2);

      fastjet::ClusterSequenceArea clust_seq(input_particles, jet_def[iR], area_def);
      double ptmin = 1.0;
      vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

       if (inclusive_jets.size() == 0) continue;

	Double_t PerpConePhiPlus  = inclusive_jets[0].phi() + TMath::Pi()/2;
        PerpConePhiPlus  = (PerpConePhiPlus>2*TMath::Pi()) ? PerpConePhiPlus - 2*TMath::Pi() : PerpConePhiPlus; // fit to 0 < phi < 2pi



	Double_t PerpConeEta = inclusive_jets[0].eta();
	Double_t PerpConePtPlus(0);

      for (unsigned int j=0; j < input_particles.size(); j++){
	
	Double_t deltaRPlus(0);

	Double_t dPhiPlus = TMath::Abs(input_particles[j].phi() - PerpConePhiPlus);
        dPhiPlus = (dPhiPlus>TMath::Pi()) ? 2*TMath::Pi()-dPhiPlus : dPhiPlus;


	Double_t dEta = TMath::Abs(input_particles[j].eta() - PerpConeEta);

	deltaRPlus  = TMath::Sqrt( TMath::Power( dEta ,2) + TMath::Power( dPhiPlus  ,2) );

	if ( deltaRPlus <= Rvals[iR] )
	   PerpConePtPlus  += input_particles[j].pt();


      }
	
      Double_t PerpConeRhoPlus    = PerpConePtPlus/(TMath::Pi()* TMath::Power(Rvals[iR],2) );

      for (unsigned int i = 0; i < inclusive_jets.size(); i++) {


       if (TMath::Abs(inclusive_jets[i].eta()) > (0.9-Rvals[iR]) ) continue;

	Double_t PtSubPlus    = inclusive_jets[i].perp() - ( PerpConeRhoPlus    * inclusive_jets[i].area()) ;

	hFullPtSpecFull[iR]->Fill(inclusive_jets[i].perp(), evt_wght);
	hFullPtSpecSubPlus[iR]->Fill( PtSubPlus , evt_wght);

        Double_t LeadingPtTrack = 0; 
         
        vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
        for ( unsigned int cns=0; cns < constituents.size(); cns++)
	    if ( constituents[cns].pt() > LeadingPtTrack) LeadingPtTrack = constituents[cns].pt();

        if (LeadingPtTrack > 5){
            hFullPtSpecTrackcut5Gev[iR]->Fill(inclusive_jets[i].perp(), evt_wght);
            hFullPtSpecTrackcut5GevSub[iR]->Fill(PtSubPlus, evt_wght);
	}
        
	if (TMath::Abs(inclusive_jets[i].eta()) < 0.5 ){
	      hFullPtSpecFullEta[iR]->Fill(inclusive_jets[i].perp(), evt_wght);
	      hFullPtSpecSubPlusEta[iR]->Fill( PtSubPlus , evt_wght);
	}

      }
    }
  }

  pythia8->PrintStatistics();

  Double_t sumw = pythia->info.weightSum();
  sumw *= 1e-9;
  Double_t TotalXSec = pythia->info.sigmaGen();

  cout<<"\nTotal Xsec is : "<<TotalXSec<<"  Total Weight is : "<<sumw<<"          "<<SumW<<endl;
  

  for (int iR = 0; iR < nR; iR++){
    
     hFullPtSpecFull[iR]->Scale( (1.0/nev) * (1.0/Etavals[iR]),"width");
     hFullPtSpecSubPlus[iR]->Scale( (1.0/nev) * (1.0/Etavals[iR]),"width");
     hFullPtSpecFullEta[iR]->Scale( (1.0/nev) *  1.0 ,"width");
     hFullPtSpecSubPlusEta[iR]->Scale( (1.0/nev) *  1.0 ,"width");
     hFullPtSpecTrackcut5Gev[iR]->Scale(  (1.0/nev) *  (1.0/Etavals[iR]),"width");
     hFullPtSpecTrackcut5GevSub[iR]->Scale(  (1.0/nev) *  (1.0/Etavals[iR]),"width");

  }

fout->Write();
}


