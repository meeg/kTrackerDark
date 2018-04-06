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

    dataTree->SetBranchAddress("rawEvent", &rawEvent);
    //dataTree->SetBranchAddress("tracklets", &tracklets);

    vector<int> hitvectors[8];

    TFile* saveFile = new TFile(argv[2], "recreate");
    //TTree* saveTree = new TTree("save", "save");

    saveFile->cd();
    TH2I* fbCorrHists[4];
    char name[500];
    for (int quad=0; quad<4; quad++) {
        sprintf(name, "fbCorrHist_%i", quad);
        fbCorrHists[quad] = new TH2I(name,name,30,0.5,30.5,50,0.5,50.5);
    }
    //TH1I* hitTimeHists[8];
    //for (int quad=0; quad<8; quad++) {
        //sprintf(name, "hitTimeHist_%i", quad);
        //hitTimeHists[quad] = new TH2I(name,name,30,0.5,30.5,50,0.5,50.5);
    //}

    for(Int_t i = 0; i < dataTree->GetEntries(); ++i)
    {
        dataTree->GetEntry(i);
        if(i % 1000 == 0) cout << i << endl;
        //cout << i << endl;

        if (rawEvent->getTriggerBits()>0 && rawEvent->getTriggerBits() & 64 != 0) {
            //if (rawEvent->getTriggerBits()<0 || rawEvent->getTriggerBits() & 64 == 0) continue;

            for(Int_t k = 0; k < rawEvent->getNHitsAll(); ++k)
            {
                Hit h = rawEvent->getHit(k);
                int dpQuadID = h.detectorID-55;
                if(dpQuadID<0 || dpQuadID>7) continue;
                if (h.isInTime())
                    hitvectors[dpQuadID].push_back(h.elementID);

                //if (h.isInTime()) {
                //cout << "in time " << dpQuadID << "  " << h.elementID << endl;
                //} else {
                //cout << "not in time " << dpQuadID << "  " << h.elementID << endl;
                //}
            }
            for (int quad=0;quad<4;quad++) {
                int fID = (quad)%4;
                int bID = quad+4;
                int fHits = hitvectors[fID].size();
                int bHits = hitvectors[bID].size();
                //cout << quad << " " << fHits << " " << bHits << endl;
                for (int fHit=0; fHit<fHits; fHit++) {
                    for (int bHit=0; bHit<bHits; bHit++) {
                        fbCorrHists[quad]->Fill(hitvectors[fID].at(fHit), hitvectors[bID].at(bHit));
                    }
                }
            }

            //saveTree->Fill();
            for (int j=0;j<8;j++) {
                hitvectors[j].clear();
            }
        }
        rawEvent->clear();
        //if (i==100000) break;
    }

    //saveTree->Write();
    saveFile->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}
