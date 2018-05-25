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
#include <TH2D.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"
#include "JobOptsSvc.h"

using namespace std;

int main(int argc, char *argv[])
{
    //JobOptsSvc::instance()->init();

    //Initialization of geometry and tracked data
    //GeomSvc* p_geomSvc = GeomSvc::instance();
    //p_geomSvc->init();

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
    Int_t nHits_st2[1000];
    Int_t nHits_st3a[1000];
    Int_t nHits_st3b[1000];
    Int_t st2Hits[1000];
    Int_t st3aHits[1000];
    Int_t st3bHits[1000];
    Double_t chisq[1000];
    Double_t tx[1000];
    Double_t ty[1000];
    Double_t x0[1000];
    Double_t y0[1000];

    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    TH2D* rtHist = new TH2D("rt","rt",1000,0,2000,100,-2.0,2.0);

    saveTree->Branch("runID", &runID);
    saveTree->Branch("spillID", &spillID);
    saveTree->Branch("eventID", &eventID);
    saveTree->Branch("triggerBits", &triggerBits);
    saveTree->Branch("nTracks", &nTracks);
    saveTree->Branch("nHits", nHits, "nHits[nTracks]/I");
    saveTree->Branch("nHits_st2", nHits_st2, "nHits_st2[nTracks]/I");
    saveTree->Branch("nHits_st3a", nHits_st3a, "nHits_st3a[nTracks]/I");
    saveTree->Branch("nHits_st3b", nHits_st3b, "nHits_st3b[nTracks]/I");
    saveTree->Branch("st2Hits", st2Hits, "st2Hits[nTracks]/I");
    saveTree->Branch("st3aHits", st3aHits, "st3aHits[nTracks]/I");
    saveTree->Branch("st3bHits", st3bHits, "st3bHits[nTracks]/I");
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
            nHits_st2[j] = 0;
            nHits_st3a[j] = 0;
            nHits_st3b[j] = 0;
            st2Hits[j] = 0;
            st3aHits[j] = 0;
            st3bHits[j] = 0;
            chisq[j] = track->chisq;
            tx[j] = track->tx;
            ty[j] = track->ty;
            x0[j] = track->x0;
            y0[j] = track->y0;
            for (int k=0;k<18;k++) {
                if (track->getSignedHit(k).hit.index != -1) {
                    switch (track->getSignedHit(k).hit.detectorID) {
                        case 1:
                        case 2:
                        case 3:
                        case 4:
                        case 5:
                        case 6:
                            printf("st1 hit!\n");
                            break;
                        case 13:
                        case 14:
                        case 15:
                        case 16:
                        case 17:
                        case 18:
                            //printf("st2 hit, detID %d time %f dist %f\n",track->getSignedHit(k).hit.detectorID,track->getSignedHit(k).hit.tdcTime,track->getSignedHit(k).hit.driftDistance);
                            //if (rawEvent->getTriggerBits()==64)
                                //rtHist->Fill(track->getSignedHit(k).hit.tdcTime,track->getSignedHit(k).hit.driftDistance);
                            //else
                                //rtHist->Fill(track->getSignedHit(k).hit.tdcTime,-track->getSignedHit(k).hit.driftDistance);
                            nHits_st2[j]++;
                            st2Hits[j] += (1 << (track->getSignedHit(k).hit.detectorID-13));
                            break;
                        case 19:
                        case 20:
                        case 21:
                        case 22:
                        case 23:
                        case 24:
                            //printf("st3a hit\n");
                            nHits_st3a[j]++;
                            st3aHits[j] += (1 << (track->getSignedHit(k).hit.detectorID-19));
                            break;
                        case 25:
                        case 26:
                        case 27:
                        case 28:
                        case 29:
                        case 30:
                            //printf("st3b hit\n");
                            nHits_st3b[j]++;
                            st3bHits[j] += (1 << (track->getSignedHit(k).hit.detectorID-25));
                            break;
                    }
                }
            }
        }

        if(nTracks > 0) saveTree->Fill();
        tracklets->Clear();
        rawEvent->clear();
    }

    saveFile->cd();
    rtHist->Write();
    saveTree->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}
