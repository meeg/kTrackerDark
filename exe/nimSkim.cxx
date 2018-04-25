#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TH1I.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"

using namespace std;

int main(int argc, char *argv[])
{
    //Initialization of geometry and tracked data
    //GeomSvc* p_geomSvc = GeomSvc::instance();
    //p_geomSvc->init();

    SRawEvent* rawEvent = new SRawEvent();

    TFile* dataFile = new TFile(argv[1], "READ");
    TTree* dataTree = (TTree*)dataFile->Get("save");

    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    dataTree->SetBranchAddress("rawEvent", &rawEvent);
    saveTree->Branch("rawEvent", &rawEvent);

    for(Int_t i = 0; i < dataTree->GetEntries(); ++i)
    {
        dataTree->GetEntry(i);
        if(i % 1000 == 0) cout << i << endl;
        //cout << i << endl;

        if (rawEvent->getTriggerBits()>0 && (rawEvent->getTriggerBits() & (32+64+128)) != 0) {
            saveTree->Fill();
        }
        rawEvent->clear();
        //if (i==100000) break;
    }

    saveTree->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}
