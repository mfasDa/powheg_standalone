#if defined(__CLING__)
R__ADD_INCLUDE_PATH($PYTHIA_ROOT/include)
R__ADD_INCLUDE_PATH($FASTJET_ROOT/include)
R__LOAD_LIBRARY(libpythia8)
R__LOAD_LIBRARY(libfastjet)
R__LOAD_LIBRARY(libsiscone);
R__LOAD_LIBRARY(libsiscone_spherical);
R__LOAD_LIBRARY(libfastjetplugins);
#else
#include <ROOT/TSeq.hxx>
#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TFile.h>
#endif

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
// using namespace Pythia8;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh" // for area support
#endif                                // __FJCORE__

bool isSelected(const Pythia8::Particle &part, bool selectFinal = true) {
  if(part.id() == 90) return false;
  if(selectFinal && !part.isFinal()) return false;
  return true;
}

 Int_t ImportParticles(Pythia8::Pythia &pythia, TClonesArray &clonesParticles, Option_t *option = "")
 {
    clonesParticles.Clear();
    Int_t nparts=0;
    Int_t i;
    Int_t ioff = 0;
    int numberOfParticles  = pythia.event.size();
    if (pythia.event[0].id() == 90) {
      ioff = -1;
    }
 
    if (!strcmp(option,"") || !strcmp(option,"Final")) {
       for (i = 0; i < numberOfParticles; i++) {
          if (pythia.event[i].id() == 90) continue;
          if (pythia.event[i].isFinal()) {
             new(clonesParticles[nparts]) TParticle(
                                                    pythia.event[i].id(),
                                                    pythia.event[i].isFinal(),
                                                    pythia.event[i].mother1() + ioff,
                                                    pythia.event[i].mother2() + ioff,
                                                    pythia.event[i].daughter1() + ioff,
                                                    pythia.event[i].daughter2() + ioff,
                                                    pythia.event[i].px(),     // [GeV/c]
                                                    pythia.event[i].py(),     // [GeV/c]
                                                    pythia.event[i].pz(),     // [GeV/c]
                                                    pythia.event[i].e(),      // [GeV]
                                                    pythia.event[i].xProd(),  // [mm]
                                                    pythia.event[i].yProd(),  // [mm]
                                                    pythia.event[i].zProd(),  // [mm]
                                                    pythia.event[i].tProd()); // [mm/c]
             nparts++;
          } // final state partice
       } // particle loop
    } else if (!strcmp(option,"All")) {
       for (i = 0; i < numberOfParticles; i++) {
          if (pythia.event[i].id() == 90) continue;
          new(clonesParticles[nparts]) TParticle(
                                                 pythia.event[i].id(),
                                                 pythia.event[i].isFinal(),
                                                 pythia.event[i].mother1() + ioff,
                                                 pythia.event[i].mother2() + ioff,
                                                 pythia.event[i].daughter1() + ioff,
                                                 pythia.event[i].daughter2() + ioff,
                                                 pythia.event[i].px(),       // [GeV/c]
                                                 pythia.event[i].py(),       // [GeV/c]
                                                 pythia.event[i].pz(),       // [GeV/c]
                                                 pythia.event[i].e(),        // [GeV]
                                                 pythia.event[i].xProd(),    // [mm]
                                                 pythia.event[i].yProd(),    // [mm]
                                                 pythia.event[i].zProd(),    // [mm]
                                                 pythia.event[i].tProd());   // [mm/c]
          nparts++;
       } // particle loop
    }
    if(ioff==-1) nparts--;
    return nparts;
 }

void RunPythia8(Char_t const *foutname = "Pythia8JetSpectra_CT14nlo.root", const char *weightname = "main", Int_t nev = -1, Int_t ndeb = 1)
{
  clock_t begin_time = clock();

  TDatime dt;
  static UInt_t sseed = dt.Get();

  if (gSystem->Getenv("CONFIG_SEED"))
  {
    sseed = atoi(gSystem->Getenv("CONFIG_SEED"));
    std::cout << "\nseed for Random number generation is : " << sseed << std::endl;
  }
  
  std::array<Int_t, 5> RVals = {{2, 3, 4, 5, 6}};

  TH1F *hNEvent = new TH1F("hNEvent", "number of events; N", 1, 0., 1.),
       *hSumWeights = new TH1F("hSumWeights", "sum of weights", 1, 0., 1.),
       *hNEventPos = new TH1F("hNEventsPos", "Number of events with positive weight", 1, 0., 1.),
       *hNEventNeg = new TH1F("hNEventsNeg", "Number of events with negative weight", 1, 0., 1.),
       *hSumWeightsPos = new TH1F("hSumWeightsPos", "sum of pos. weights", 1, 0., 1.),
       *hSumWeightsNeg = new TH1F("hSumWeightsNeg", "sum of neg. weights", 1, 0., 1.);
  std::map<int, TH1*> hFullPtSpecFull, hFullPtSpecSubPlus, hFullPtSpecFullEta, hFullPtSpecSubPlusEta, hFullPtSpecTrackcut5Gev, hFullPtSpecTrackcut5GevSub;

  for (auto R : RVals)
  {
    double radius = double(R)/10.;
    hFullPtSpecFull[R] = new TH1F(Form("hFullPtSpecFullR%02d", R),
                                   Form("Full jet cross section R=%.1f with small bin width;P_{T,jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)", radius),
                                   200, 0, 200);
    hFullPtSpecFull[R]->Sumw2();

    hFullPtSpecSubPlus[R] = new TH1F(Form("hSubPtSpecR%02d", R),
                                      Form("Subtracted jet cross section R=%.1f with small bin width (plus);P_{T,jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)", radius),
                                      250, -50, 200);
    hFullPtSpecSubPlus[R]->Sumw2();

    hFullPtSpecFullEta[R] = new TH1F(Form("hFullPtSpecFullR%02dEta05", R),
                                      Form("Full jet cross section R=%.1f with small bin width Eta05;P_{T,jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)", radius),
                                      200, 0, 200);
    hFullPtSpecFullEta[R]->Sumw2();

    hFullPtSpecSubPlusEta[R] = new TH1F(Form("hSubPtSpecR%02dEta05", R),
                                         Form("Subtracted jet cross section R=%.1f with small bin width (plus) Eta05;P_{T,jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)",radius),
                                         250, -50, 200);
    hFullPtSpecSubPlusEta[R]->Sumw2();

    hFullPtSpecTrackcut5Gev[R] = new TH1F(Form("hFullPtSpecTrackcut5Gev_R%02d", R),
                                           Form("jet cross section R=%.1f with Leading Pt cut at 5Gev;P_{T,jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)", radius),
                                           200, 0, 200);
    hFullPtSpecTrackcut5Gev[R]->Sumw2();

    hFullPtSpecTrackcut5GevSub[R] = new TH1F(Form("hSubPtSpecTrackcut5Gev_R%02d", R),
                                              Form("jet cross section R=%.1f with Leading Pt cut at 5Gev;P_{T,jet}(Gev/c);d#sigma/dP_{T}d#eta(mb c/Gev)", radius),
                                              200, 0, 200);
    hFullPtSpecTrackcut5GevSub[R]->Sumw2();
  }

  Double_t ghost_maxrap = 6.0;
  fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  // Configure
  Pythia8::Pythia pythia;
  pythia.readString("Next:numberShowLHA = 1");
  pythia.readString("Next:numberShowInfo = 1");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 1");
  pythia.readString("Main:timesAllowErrors = 10");

  pythia.readString("Init:showChangedSettings = on");
  pythia.readString("Init:showChangedParticleData = off");

  pythia.readString("Beams:frametype = 4");

  pythia.readString("Beams:LHEF = pwgevents.lhe");

  pythia.readString("POWHEG:nFinal = 2");

  pythia.readString("PartonLevel:MPI = on"); //! TEST

  pythia.readString("111:mayDecay  = on");
  pythia.readString("310:mayDecay  = off");
  pythia.readString("3122:mayDecay = off");
  pythia.readString("3112:mayDecay = off");
  pythia.readString("3212:mayDecay = off");
  pythia.readString("3222:mayDecay = off");
  pythia.readString("3312:mayDecay = off");
  pythia.readString("3322:mayDecay = off");
  pythia.readString("3334:mayDecay = off");

  pythia.readString("POWHEG:veto = 1");
  pythia.readString("POWHEG:vetoCount = 1000000");
  pythia.readString("POWHEG:pThard = 2"); //!
  pythia.readString("POWHEG:pTemt = 0");
  pythia.readString("POWHEG:emitted = 3"); //!
  pythia.readString("POWHEG:pTdef = 1");
  pythia.readString("POWHEG:MPIveto = 1"); //!
  pythia.readString("POWHEG:QEDveto = 2");

  pythia.readString("Tune:preferLHAPDF = 2");
  pythia.readString("Tune:pp = 5");

  pythia.readString("PDF:pSet = LHAPDF6:CT14nlo");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %u", sseed % 900000000));

  // Add in user hooks for shower vetoing
  Pythia8::PowhegHooks *powhegHooks = new Pythia8::PowhegHooks();
  pythia.readString("SpaceShower:pTmaxMatch = 2");
  pythia.readString("TimeShower:pTmaxMatch = 2");
  pythia.readString("MultipartonInteractions:pTmaxMatch = 2"); //!
  pythia.setUserHooksPtr((Pythia8::UserHooks *)powhegHooks);

  TString PDFused = pythia.settings.word("PDF:pSet");
  std::cout << "\n PDF used is : " << PDFused << std::endl;
  Double_t SumW(0);

  // Initialize
  pythia.init();
  //

  const double maxeta = 0.7;

  // Event loop
  int iev = 0, acceptedevents = 0;
  while (true)
  {
    if(nev >= 0 && iev >= nev ) {
      std::cout << "Got requested " << iev << " events. Exiting ..." << std::endl;
      break;
    }

    if (!(iev % 1000))
    {
      printf(">>>processing ev# %5d / elapsed time: ", iev);
      std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
      begin_time = clock();
    }

    if(!pythia.next()) {
      std::cout << " No more events, exiting ..." << std::endl;
      break;
    }
    // pythia8->EventListing();

    Double_t evt_wght = pythia.info.getWeightsDetailedValue(weightname);
    if(std::isnan(evt_wght)) {
      std::cerr << "Event with nan-weight found" << std::endl;
      continue;
    }
    evt_wght *= 1e-9;
    SumW += evt_wght;
    //cout<<" Event weight is : "<<evt_wght*1e9<<endl;
    hNEvent->Fill(0.5, 1.);
    hSumWeights->Fill(0.5, evt_wght);
    if(evt_wght>=0){
      hNEventPos->Fill(0.5, 1.);
      hSumWeightsPos->Fill(0.5, evt_wght);
    } else {
      hNEventNeg->Fill(0.5, 1.);
      hSumWeightsNeg->Fill(0.5, evt_wght);
    }
    iev++;
    if(evt_wght < 0) {
      std::cout << "Found event with negative weight - continue" << std::endl;
      continue;    // Skip events with negative weights
    }
    acceptedevents++;

    // Particle loop
    std::vector<fastjet::PseudoJet> input_particles;
    for(auto ipart : ROOT::TSeqI(0, pythia.event.size()))
    {
      auto part = pythia.event[ipart];
      if(!isSelected(part)) continue;

      Float_t charge = part.charge();
      Float_t eta = part.eta();
      Float_t pt = part.pT();

      if (abs(eta) > maxeta)
        continue;
      if (pt < 0.15)
        continue;

      input_particles.push_back(fastjet::PseudoJet(part.px(), part.py(), part.pz(), part.e()));
    }

    if (input_particles.size() == 0)
    {
      //printf("No particle....\n");
      continue;
    }

    for (auto R : RVals)
    {
      double radius = double(R)/10.;

      fastjet::ClusterSequenceArea clust_seq(input_particles, fastjet::JetDefinition(fastjet::antikt_algorithm, radius), area_def);
      double ptmin = 1.0;
      std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

      if (inclusive_jets.size() == 0)
        continue;

      Double_t PerpConePhiPlus = inclusive_jets[0].phi() + TMath::Pi() / 2;
      PerpConePhiPlus = (PerpConePhiPlus > 2 * TMath::Pi()) ? PerpConePhiPlus - 2 * TMath::Pi() : PerpConePhiPlus; // fit to 0 < phi < 2pi

      Double_t PerpConeEta = inclusive_jets[0].eta();
      Double_t PerpConePtPlus(0);

      for (unsigned int j = 0; j < input_particles.size(); j++)
      {

        Double_t deltaRPlus(0);

        Double_t dPhiPlus = TMath::Abs(input_particles[j].phi() - PerpConePhiPlus);
        dPhiPlus = (dPhiPlus > TMath::Pi()) ? 2 * TMath::Pi() - dPhiPlus : dPhiPlus;

        Double_t dEta = TMath::Abs(input_particles[j].eta() - PerpConeEta);

        deltaRPlus = TMath::Sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhiPlus, 2));

        if (deltaRPlus <= radius)
          PerpConePtPlus += input_particles[j].pt();
      }

      Double_t PerpConeRhoPlus = PerpConePtPlus / (TMath::Pi() * TMath::Power(radius, 2));

      for (unsigned int i = 0; i < inclusive_jets.size(); i++)
      {

        if (TMath::Abs(inclusive_jets[i].eta()) > (maxeta - radius))
          continue;

        Double_t PtSubPlus = inclusive_jets[i].perp() - (PerpConeRhoPlus * inclusive_jets[i].area());

        hFullPtSpecFull[R]->Fill(inclusive_jets[i].perp(), evt_wght);
        hFullPtSpecSubPlus[R]->Fill(PtSubPlus, evt_wght);

        Double_t LeadingPtTrack = 0;

        std::vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
        for (unsigned int cns = 0; cns < constituents.size(); cns++)
          if (constituents[cns].pt() > LeadingPtTrack)
            LeadingPtTrack = constituents[cns].pt();

        if (LeadingPtTrack > 5)
        {
          hFullPtSpecTrackcut5Gev[R]->Fill(inclusive_jets[i].perp(), evt_wght);
          hFullPtSpecTrackcut5GevSub[R]->Fill(PtSubPlus, evt_wght);
        }

        if (TMath::Abs(inclusive_jets[i].eta()) < 0.5)
        {
          hFullPtSpecFullEta[R]->Fill(inclusive_jets[i].perp(), evt_wght);
          hFullPtSpecSubPlusEta[R]->Fill(PtSubPlus, evt_wght);
        }
      }
    }
  }

  pythia.stat();

  Double_t sumw = pythia.info.weightSum();
  sumw *= 1e-9;
  Double_t TotalXSec = pythia.info.sigmaGen();

  std::cout << "Generated " << iev << " events" << std::endl;
  std::cout << "\nTotal Xsec is : " << TotalXSec << "  Total Weight is : " << sumw << "          " << SumW << std::endl;

  for (auto R : RVals)
  {
    double radius = double(R)/10.;
    double eta = 2*maxeta - 2*radius;
    hFullPtSpecFull[R]->Scale((1.0 / eta), "width");
    hFullPtSpecSubPlus[R]->Scale((1.0 / eta), "width");
    hFullPtSpecFullEta[R]->Scale(1.0, "width");
    hFullPtSpecSubPlusEta[R]->Scale(1.0, "width");
    hFullPtSpecTrackcut5Gev[R]->Scale((1.0 / eta), "width");
    hFullPtSpecTrackcut5GevSub[R]->Scale((1.0 / eta), "width");
  }

  TFile *fout = new TFile(foutname, "RECREATE");
  hNEvent->Write();
  hNEventPos->Write();
  hNEventNeg->Write();
  hSumWeights->Write();
  hSumWeightsPos->Write();
  hSumWeightsNeg->Write();
  for (auto R : RVals)
  {
    hFullPtSpecFull[R]->Write();
    hFullPtSpecSubPlus[R]->Write();
    hFullPtSpecFullEta[R]->Write();
    hFullPtSpecSubPlusEta[R]->Write();
    hFullPtSpecTrackcut5Gev[R]->Write();
    hFullPtSpecTrackcut5GevSub[R]->Write();
  }
}
