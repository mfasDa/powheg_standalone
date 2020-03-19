#ifndef __CLING__
#include <memory>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TSystem.h>
#endif // __CLING__

std::vector<std::string> getPythiaFiles(const char *pythiafilename) {
    TString dirstring = gSystem->GetFromPipe(Form("find . -name %s", pythiafilename));
    std::unique_ptr<TObjArray> files(dirstring.Tokenize("\n"));
    std::vector<std::string> outputfiles;
    for(auto fname : TRangeDynCast<TObjString>(files.get())) {
        outputfiles.push_back(fname->String().Data());
    }
    return outputfiles;
}

int checkFile(const char *filename){
    std::cout << "Checking file " << filename << std::endl;
    std::unique_ptr<TFile> reader(TFile::Open(filename, "READ"));
    if(!reader.get() || reader->IsZombie()) {
        std::cerr << "File corrupted" << std::endl;
        return 0;
    }
    // get the number of positive weighted events
    auto heventspos = static_cast<TH1 *>(reader->Get("hNEventsPos"))->GetEntries();
    std::cout << "Number of events (pos. weight): " << heventspos << std::endl; 
    // check jet spectrum and for R = 0.2
    auto jetspec = static_cast<TH1 *>(reader->Get("hFullPtSpecFullR02"));
    std::cout << "Entries in jet spectrum for R=0.2: " << jetspec->GetEntries() << ", min " << jetspec->GetMinimum() << ", max " << jetspec->GetMaximum() << std::endl;

    return heventspos;
}

void checkPythiaFiles(const char *pythiafilename){
    std::vector<int> nevents;
    for(auto fl : getPythiaFiles(pythiafilename)) {
        nevents.push_back(checkFile(fl.data()));
    }
    std::cout << "Mean number of events: " << TMath::Mean(nevents.begin(), nevents.end()) << std::endl;
}