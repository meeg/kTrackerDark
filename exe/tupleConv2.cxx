#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TH1I.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"
#include "JobOptsSvc.h"

using namespace std;

int main(int argc, char *argv[])
{
    JobOptsSvc::instance()->init();

    //Initialization of geometry and tracked data
    GeomSvc* p_geomSvc = GeomSvc::instance();
    p_geomSvc->init();

    SRawEvent* rawEvent = new SRawEvent();
    TClonesArray* tracklets = new TClonesArray("Tracklet");

    TFile* dataFile = new TFile(argv[1], "READ");
    TTree* dataTree = (TTree*)dataFile->Get("save");

    dataTree->SetBranchAddress("rawEvent", &rawEvent);
    dataTree->SetBranchAddress("tracklets", &tracklets);

    Int_t runID, spillID, eventID;
    Int_t triggerBits;
    Int_t nTracks;
    Int_t nHits[1000];
    Double_t chisq[1000];
    Double_t tx[1000];
    Double_t ty[1000];
    Double_t x0[1000];
    Double_t y0[1000];

    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    saveTree->Branch("runID", &runID);
    saveTree->Branch("spillID", &spillID);
    saveTree->Branch("eventID", &eventID);
    saveTree->Branch("triggerBits", &triggerBits);
    saveTree->Branch("nTracks", &nTracks);
    saveTree->Branch("nHits", nHits, "nHits[nTracks]/I");
    saveTree->Branch("chisq", chisq, "chisq[nTracks]/D");
    saveTree->Branch("tx", tx, "tx[nTracks]/D");
    saveTree->Branch("ty", ty, "ty[nTracks]/D");
    saveTree->Branch("x0", x0, "x0[nTracks]/D");
    saveTree->Branch("y0", y0, "y0[nTracks]/D");

    for(Int_t i = 0; i < dataTree->GetEntries(); ++i)
    {
        dataTree->GetEntry(i);
        if(i % 10000 == 0) cout << i << endl;

        //Initial numbers
        runID = rawEvent->getRunID();
        spillID = rawEvent->getSpillID();
        eventID = rawEvent->getEventID();
        triggerBits = rawEvent->getTriggerBits();
        nTracks = 0;

        nTracks = tracklets->GetEntries();
        for(Int_t j = 0 ; j < nTracks; ++j)
        {
            Tracklet* track = (Tracklet*)tracklets->At(j);

            nHits[j] = track->getNHits();
            chisq[j] = track->chisq;
            tx[j] = track->tx;
            ty[j] = track->ty;
            x0[j] = track->x0;
            y0[j] = track->y0;
        }

        if(nTracks > 0) saveTree->Fill();
        tracklets->Clear();
        rawEvent->clear();
    }

    saveFile->cd();
    saveTree->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}