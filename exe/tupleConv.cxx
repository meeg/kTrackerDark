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
#include <TCanvas.h>
#include <TMath.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"
#include "JobOptsSvc.h"

using namespace std;

typedef struct
{
    Int_t nHits;
    Char_t nHits_st1, nHits_st2, nHits_st3a, nHits_st3b;
    Char_t st1Hits;
    Double_t chisq, tx_st1, tx, ty;
    Double_t x_st1, x, y;
    Double_t pz;
} trackdata;

int main(int argc, char *argv[])
{
    JobOptsSvc::instance()->init();

    //Initialization of geometry and tracked data
    GeomSvc* p_geomSvc = GeomSvc::instance();
    p_geomSvc->init();

    SRawEvent* rawEvent = new SRawEvent();
    TClonesArray* tracklets = new TClonesArray("Tracklet");

    vector<trackdata> tracksPos;
    vector<trackdata> tracksNeg;

    TFile* dataFile = new TFile(argv[1], "READ");
    TTree* dataTree = (TTree*)dataFile->Get("save");

    dataTree->SetBranchAddress("rawEvent", &rawEvent);
    dataTree->SetBranchAddress("tracklets", &tracklets);

    Int_t runID, spillID, eventID;
    Int_t triggerBits;
    Int_t nTracks;
    Int_t nHits1, nHits2;
    Char_t nHits1_st1, nHits1_st2, nHits1_st3a, nHits1_st3b;
    Char_t nHits2_st1, nHits2_st2, nHits2_st3a, nHits2_st3b;
    Char_t st1Hits1, st1Hits2;
    Double_t chisq1, chisq2;
    Double_t tx1_st1, tx1, ty1;
    Double_t tx2_st1, tx2, ty2;
    Double_t x1_st1, x1, y1;
    Double_t x2_st1, x2, y2;
    Double_t pz1, pz2;
    Char_t pairNum;

    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    saveTree->Branch("runID", &runID);
    saveTree->Branch("spillID", &spillID);
    saveTree->Branch("eventID", &eventID);
    saveTree->Branch("triggerBits", &triggerBits);
    saveTree->Branch("nTracks", &nTracks);
    saveTree->Branch("nHits1", &nHits1);
    saveTree->Branch("nHits2", &nHits2);
    saveTree->Branch("nHits1_st1", &nHits1_st1);
    saveTree->Branch("nHits1_st2", &nHits1_st2);
    saveTree->Branch("nHits1_st3a", &nHits1_st3a);
    saveTree->Branch("nHits1_st3b", &nHits1_st3b);
    saveTree->Branch("nHits2_st1", &nHits2_st1);
    saveTree->Branch("nHits2_st2", &nHits2_st2);
    saveTree->Branch("nHits2_st3a", &nHits2_st3a);
    saveTree->Branch("nHits2_st3b", &nHits2_st3b);
    saveTree->Branch("st1Hits1", &st1Hits1);
    saveTree->Branch("st1Hits2", &st1Hits2);
    saveTree->Branch("chisq1", &chisq1);
    saveTree->Branch("chisq2", &chisq2);
    saveTree->Branch("tx1_st1", &tx1_st1);
    saveTree->Branch("tx1", &tx1);
    saveTree->Branch("ty1", &ty1);
    saveTree->Branch("x1_st1", &x1_st1);
    saveTree->Branch("x1", &x1);
    saveTree->Branch("y1", &y1);
    saveTree->Branch("tx2_st1", &tx2_st1);
    saveTree->Branch("tx2", &tx2);
    saveTree->Branch("ty2", &ty2);
    saveTree->Branch("x2_st1", &x2_st1);
    saveTree->Branch("x2", &x2);
    saveTree->Branch("y2", &y2);
    saveTree->Branch("pz1", &pz1);
    saveTree->Branch("pz2", &pz2);
    saveTree->Branch("pairNum", &pairNum);

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
        nHits1 = 0;
        nHits2 = 0;
        nHits1_st1 = -1;
        nHits1_st2 = -1;
        nHits1_st3a = -1;
        nHits1_st3b = -1;
        nHits2_st1 = -1;
        nHits2_st2 = -1;
        nHits2_st3a = -1;
        nHits2_st3b = -1;
        st1Hits1 = -1;
        st1Hits2 = -1;
        chisq1 = 999.;
        chisq2 = 999.;
        tx1_st1 = 999.;
        tx1 = 999.;
        ty1 = 999.;
        x1_st1 = 999.;
        x1 = 999.;
        y1 = 999.;
        pz1 = 999.;
        tx2_st1 = 999.;
        tx2 = 999.;
        ty2 = 999.;
        x2_st1 = 999.;
        x2 = 999.;
        y2 = 999.;
        pz2 = 999.;
        pairNum = 0;

        Int_t nTracksTotal = tracklets->GetEntries();

        for(Int_t j = 0 ; j < nTracksTotal; ++j) // count tracks
        {
            Tracklet* track = (Tracklet*)tracklets->At(j);

            if(track->getNHits() < 14) continue;
            if(track->chisq/(track->getNHits() - 5.) > 12.) continue;
            nTracks++;
            trackdata newtrack;
            newtrack.nHits_st1 = 0;
            newtrack.nHits_st2 = 0;
            newtrack.nHits_st3a = 0;
            newtrack.nHits_st3b = 0;
            newtrack.st1Hits = 0;
            for (int k=0;k<18;k++) {
                if (track->getSignedHit(k).hit.index != -1) {
                    switch (track->getSignedHit(k).hit.detectorID) {
                        case 1:
                        case 2:
                        case 3:
                        case 4:
                        case 5:
                        case 6:
                            newtrack.nHits_st1++;
                            newtrack.st1Hits += (1 << (track->getSignedHit(k).hit.detectorID-1));
                            break;
                        case 13:
                        case 14:
                        case 15:
                        case 16:
                        case 17:
                        case 18:
                            newtrack.nHits_st2++;
                            break;
                        case 19:
                        case 20:
                        case 21:
                        case 22:
                        case 23:
                        case 24:
                            newtrack.nHits_st3a++;
                            break;
                        case 25:
                        case 26:
                        case 27:
                        case 28:
                        case 29:
                        case 30:
                            newtrack.nHits_st3b++;
                            break;
                    }
                }
            }
            newtrack.nHits = track->getNHits();
            newtrack.chisq = track->chisq;
            newtrack.tx = track->tx;
            newtrack.ty = track->ty;
            newtrack.x = track->x0;
            newtrack.y = track->y0;
            track->getXZInfoInSt1(newtrack.tx_st1, newtrack.x_st1);
            newtrack.pz = 1./track->invP*sqrt(1. + newtrack.tx_st1*newtrack.tx_st1)/sqrt(1. + newtrack.tx_st1*newtrack.tx_st1 + newtrack.ty*newtrack.ty);

            if(track->getCharge() > 0)
                tracksPos.push_back(newtrack);
            else
                tracksNeg.push_back(newtrack);
        }
        if (tracksPos.empty()) {
            for (Int_t j=tracksNeg.size()-1;j>=0;--j) {
                nHits2 = tracksNeg.at(j).nHits;
                nHits2_st1 = tracksNeg.at(j).nHits_st1;
                nHits2_st2 = tracksNeg.at(j).nHits_st2;
                nHits2_st3a = tracksNeg.at(j).nHits_st3a;
                nHits2_st3b = tracksNeg.at(j).nHits_st3b;
                st1Hits2 = tracksNeg.at(j).st1Hits;
                chisq2 = tracksNeg.at(j).chisq;
                tx2_st1 = tracksNeg.at(j).tx_st1;
                tx2 = tracksNeg.at(j).tx;
                ty2 = tracksNeg.at(j).ty;
                x2_st1 = tracksNeg.at(j).x_st1;
                x2 = tracksNeg.at(j).x;
                y2 = tracksNeg.at(j).y;
                pz2 = tracksNeg.at(j).pz;

                saveTree->Fill();
                pairNum++;
            }
        } else if (tracksNeg.empty()) {
            for (Int_t j=tracksPos.size()-1;j>=0;--j) {
                nHits1 = tracksPos.at(j).nHits;
                nHits1_st1 = tracksPos.at(j).nHits_st1;
                nHits1_st2 = tracksPos.at(j).nHits_st2;
                nHits1_st3a = tracksPos.at(j).nHits_st3a;
                nHits1_st3b = tracksPos.at(j).nHits_st3b;
                st1Hits1 = tracksPos.at(j).st1Hits;
                chisq1 = tracksPos.at(j).chisq;
                tx1_st1 = tracksPos.at(j).tx_st1;
                tx1 = tracksPos.at(j).tx;
                ty1 = tracksPos.at(j).ty;
                x1_st1 = tracksPos.at(j).x_st1;
                x1 = tracksPos.at(j).x;
                y1 = tracksPos.at(j).y;
                pz1 = tracksPos.at(j).pz;

                saveTree->Fill();
                pairNum++;
            }
        } else {
            for (Int_t jPos=tracksPos.size()-1;jPos>=0;--jPos) {
                nHits1 = tracksPos.at(jPos).nHits;
                nHits1_st1 = tracksPos.at(jPos).nHits_st1;
                nHits1_st2 = tracksPos.at(jPos).nHits_st2;
                nHits1_st3a = tracksPos.at(jPos).nHits_st3a;
                nHits1_st3b = tracksPos.at(jPos).nHits_st3b;
                st1Hits1 = tracksPos.at(jPos).st1Hits;
                chisq1 = tracksPos.at(jPos).chisq;
                tx1_st1 = tracksPos.at(jPos).tx_st1;
                tx1 = tracksPos.at(jPos).tx;
                ty1 = tracksPos.at(jPos).ty;
                x1_st1 = tracksPos.at(jPos).x_st1;
                x1 = tracksPos.at(jPos).x;
                y1 = tracksPos.at(jPos).y;
                pz1 = tracksPos.at(jPos).pz;
                for (Int_t jNeg=tracksNeg.size()-1;jNeg>=0;--jNeg) {
                    nHits2 = tracksNeg.at(jNeg).nHits;
                    nHits2_st1 = tracksNeg.at(jNeg).nHits_st1;
                    nHits2_st2 = tracksNeg.at(jNeg).nHits_st2;
                    nHits2_st3a = tracksNeg.at(jNeg).nHits_st3a;
                    nHits2_st3b = tracksNeg.at(jNeg).nHits_st3b;
                    st1Hits2 = tracksNeg.at(jNeg).st1Hits;
                    chisq2 = tracksNeg.at(jNeg).chisq;
                    tx2_st1 = tracksNeg.at(jNeg).tx_st1;
                    tx2 = tracksNeg.at(jNeg).tx;
                    ty2 = tracksNeg.at(jNeg).ty;
                    x2_st1 = tracksNeg.at(jNeg).x_st1;
                    x2 = tracksNeg.at(jNeg).x;
                    y2 = tracksNeg.at(jNeg).y;
                    pz2 = tracksNeg.at(jNeg).pz;

                    saveTree->Fill();
                    pairNum++;
                }
            }
        }

        tracksPos.clear();
        tracksNeg.clear();
        tracklets->Clear();
        rawEvent->clear();
    }

    saveFile->cd();
    saveTree->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}
