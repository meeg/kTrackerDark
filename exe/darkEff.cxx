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

using namespace std;

int main(int argc, char *argv[])
{
    //Initialization of geometry and tracked data
    GeomSvc* p_geomSvc = GeomSvc::instance();
    p_geomSvc->init();

    SRawEvent* rawEvent = new SRawEvent();
    TClonesArray* tracklets = new TClonesArray("Tracklet");

    TFile* dataFile = new TFile(argv[1], "READ");
    TTree* dataTree = (TTree*)dataFile->Get("save");

    dataTree->SetBranchAddress("rawEvent", &rawEvent);
    dataTree->SetBranchAddress("tracklets", &tracklets);

    Plane dpPlanes[8];

    //DPLT1
    dpPlanes[0].detectorID = 55;
    dpPlanes[0].detectorName = "DPLT1";
    dpPlanes[0].planeType = 4;
    dpPlanes[0].nElements = 80;
    dpPlanes[0].spacing = 1.0099 + 0.0181967;
    dpPlanes[0].cellWidth = 1.0099 + 0.0181967;
    dpPlanes[0].xoffset = 0.;
    dpPlanes[0].overlap = 0.;
    dpPlanes[0].angleFromVert = TMath::Pi()/2.;
    dpPlanes[0].x0 = -1.48 + 40.;
    dpPlanes[0].y0 = 0.145 + 3.003*2.54 + dpPlanes[0].nElements*dpPlanes[0].cellWidth/2. + 0.059193;
    dpPlanes[0].z0 = 797.86;
    dpPlanes[0].x1 = dpPlanes[0].x0 - 40.;
    dpPlanes[0].y1 = dpPlanes[0].y0 - dpPlanes[0].nElements*dpPlanes[0].cellWidth/2.;
    dpPlanes[0].x2 = dpPlanes[0].x0 + 40.;
    dpPlanes[0].y2 = dpPlanes[0].y0 + dpPlanes[0].nElements*dpPlanes[0].cellWidth/2.;
    dpPlanes[0].thetaX = -0.005914;
    dpPlanes[0].thetaY = -0.000572;
    dpPlanes[0].thetaZ = 0.005269;
    dpPlanes[0].update();

    //DPLRT1
    dpPlanes[1].detectorID = 56;
    dpPlanes[1].detectorName = "DPRT1";
    dpPlanes[1].planeType = 4;
    dpPlanes[1].nElements = 80;
    dpPlanes[1].spacing = 1.0099 + 0.0181967;
    dpPlanes[1].cellWidth = 1.0099 + 0.0181967;
    dpPlanes[1].xoffset = 0.;
    dpPlanes[1].overlap = 0.;
    dpPlanes[1].angleFromVert = TMath::Pi()/2.;
    dpPlanes[1].x0 = -1.48 - 40.;
    dpPlanes[1].y0 = 0.145 + 3.003*2.54 + dpPlanes[1].nElements*dpPlanes[1].cellWidth/2. - 0.6556508;
    dpPlanes[1].z0 = 797.86;
    dpPlanes[1].x1 = dpPlanes[1].x0 - 40.;
    dpPlanes[1].y1 = dpPlanes[1].y0 - dpPlanes[1].nElements*dpPlanes[1].cellWidth/2.;
    dpPlanes[1].x2 = dpPlanes[1].x0 + 40.;
    dpPlanes[1].y2 = dpPlanes[1].y0 + dpPlanes[1].nElements*dpPlanes[1].cellWidth/2.;
    dpPlanes[1].thetaX = -0.005914;
    dpPlanes[1].thetaY = -0.000572;
    dpPlanes[1].thetaZ = 0.005269;
    dpPlanes[1].update();

    //DPLLB1
    dpPlanes[2].detectorID = 57;
    dpPlanes[2].detectorName = "DPLB1";
    dpPlanes[2].planeType = 4;
    dpPlanes[2].nElements = 80;
    dpPlanes[2].spacing = 1.0099 + 0.0181967;
    dpPlanes[2].cellWidth = 1.0099 + 0.0181967;
    dpPlanes[2].xoffset = 0.;
    dpPlanes[2].overlap = 0.;
    dpPlanes[2].angleFromVert = TMath::Pi()/2.;
    dpPlanes[2].x0 = -1.48 + 40.;
    dpPlanes[2].y0 = 0.145 - 3.003*2.54 - dpPlanes[2].nElements*dpPlanes[2].cellWidth/2. + 1.52231;
    dpPlanes[2].z0 = 797.86;
    dpPlanes[2].x1 = dpPlanes[2].x0 - 40.;
    dpPlanes[2].y1 = dpPlanes[2].y0 - dpPlanes[2].nElements*dpPlanes[2].cellWidth/2.;
    dpPlanes[2].x2 = dpPlanes[2].x0 + 40.;
    dpPlanes[2].y2 = dpPlanes[2].y0 + dpPlanes[2].nElements*dpPlanes[2].cellWidth/2.;
    dpPlanes[2].thetaX = -0.005914;
    dpPlanes[2].thetaY = -0.000572;
    dpPlanes[2].thetaZ = 0.005269;
    dpPlanes[2].update();

    //DPLRB1
    dpPlanes[3].detectorID = 58;
    dpPlanes[3].detectorName = "DPRB1";
    dpPlanes[3].planeType = 4;
    dpPlanes[3].nElements = 80;
    dpPlanes[3].spacing = 1.0099 + 0.0181967;
    dpPlanes[3].cellWidth = 1.0099 + 0.0181967;
    dpPlanes[3].xoffset = 0.;
    dpPlanes[3].overlap = 0.;
    dpPlanes[3].angleFromVert = TMath::Pi()/2.;
    dpPlanes[3].x0 = -1.48 - 40.;
    dpPlanes[3].y0 = 0.145 - 3.003*2.54 - dpPlanes[3].nElements*dpPlanes[3].cellWidth/2. + 0.952712;
    dpPlanes[3].z0 = 797.86;
    dpPlanes[3].x1 = dpPlanes[3].x0 - 40.;
    dpPlanes[3].y1 = dpPlanes[3].y0 - dpPlanes[3].nElements*dpPlanes[3].cellWidth/2.;
    dpPlanes[3].x2 = dpPlanes[3].x0 + 40.;
    dpPlanes[3].y2 = dpPlanes[3].y0 + dpPlanes[3].nElements*dpPlanes[3].cellWidth/2.;
    dpPlanes[3].thetaX = -0.005914;
    dpPlanes[3].thetaY = -0.000572;
    dpPlanes[3].thetaZ = 0.005269;
    dpPlanes[3].update();

    //DPLT2
    dpPlanes[4].detectorID = 59;
    dpPlanes[4].detectorName = "DPLT2";
    dpPlanes[4].planeType = 4;
    dpPlanes[4].nElements = 50;
    dpPlanes[4].spacing = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[4].cellWidth = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[4].xoffset = 0.;
    dpPlanes[4].overlap = 0.;
    dpPlanes[4].angleFromVert = TMath::Pi()/2.;
    dpPlanes[4].x0 = 0.747 + 50.;
    dpPlanes[4].y0 = 0.691 + 3.0083*2.54 + dpPlanes[4].nElements*dpPlanes[4].cellWidth/2. - 2.02706;
    dpPlanes[4].z0 = 1497.88;
    dpPlanes[4].x1 = dpPlanes[4].x0 - 50.;
    dpPlanes[4].y1 = dpPlanes[4].y0 - dpPlanes[4].nElements*dpPlanes[4].cellWidth/2.;
    dpPlanes[4].x2 = dpPlanes[4].x0 + 50.;
    dpPlanes[4].y2 = dpPlanes[4].y0 + dpPlanes[4].nElements*dpPlanes[4].cellWidth/2.;
    dpPlanes[4].thetaX = 0.015964;
    dpPlanes[4].thetaY = -0.000768;
    dpPlanes[4].thetaZ = -0.001297;
    dpPlanes[4].update();

    //DPRT2
    dpPlanes[5].detectorID = 60;
    dpPlanes[5].detectorName = "DPRT2";
    dpPlanes[5].planeType = 4;
    dpPlanes[5].nElements = 50;
    dpPlanes[5].spacing = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[5].cellWidth = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[5].xoffset = 0.;
    dpPlanes[5].overlap = 0.;
    dpPlanes[5].angleFromVert = TMath::Pi()/2.;
    dpPlanes[5].x0 = 0.747 - 50.;
    dpPlanes[5].y0 = 0.691 + 3.0083*2.54 + dpPlanes[5].nElements*dpPlanes[5].cellWidth/2. - 2.04749;
    dpPlanes[5].z0 = 1497.88;
    dpPlanes[5].x1 = dpPlanes[5].x0 - 50.;
    dpPlanes[5].y1 = dpPlanes[5].y0 - dpPlanes[5].nElements*dpPlanes[5].cellWidth/2.;
    dpPlanes[5].x2 = dpPlanes[5].x0 + 50.;
    dpPlanes[5].y2 = dpPlanes[5].y0 + dpPlanes[5].nElements*dpPlanes[5].cellWidth/2.;
    dpPlanes[5].thetaX = 0.015964;
    dpPlanes[5].thetaY = -0.000768;
    dpPlanes[5].thetaZ = -0.001297;
    dpPlanes[5].update();

    //DPLB2
    dpPlanes[6].detectorID = 61;
    dpPlanes[6].detectorName = "DPLB2";
    dpPlanes[6].planeType = 4;
    dpPlanes[6].nElements = 50;
    dpPlanes[6].spacing = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[6].cellWidth = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[6].xoffset = 0.;
    dpPlanes[6].overlap = 0.;
    dpPlanes[6].angleFromVert = TMath::Pi()/2.;
    dpPlanes[6].x0 = 0.747 + 50.;
    dpPlanes[6].y0 = 0.691 - 3.0083*2.54 - dpPlanes[6].nElements*dpPlanes[6].cellWidth/2. + 1.30634;
    dpPlanes[6].z0 = 1497.88;
    dpPlanes[6].x1 = dpPlanes[6].x0 - 50.;
    dpPlanes[6].y1 = dpPlanes[6].y0 - dpPlanes[6].nElements*dpPlanes[6].cellWidth/2.;
    dpPlanes[6].x2 = dpPlanes[6].x0 + 50.;
    dpPlanes[6].y2 = dpPlanes[6].y0 + dpPlanes[6].nElements*dpPlanes[6].cellWidth/2.;
    dpPlanes[6].thetaX = 0.015964;
    dpPlanes[6].thetaY = -0.000768;
    dpPlanes[6].thetaZ = -0.001297;
    dpPlanes[6].update();

    //DPRB2
    dpPlanes[7].detectorID = 62;
    dpPlanes[7].detectorName = "DPRB2";
    dpPlanes[7].planeType = 4;
    dpPlanes[7].nElements = 50;
    dpPlanes[7].spacing = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[7].cellWidth = 1.8999 + 1.32615500000000004e-02;
    dpPlanes[7].xoffset = 0.;
    dpPlanes[7].overlap = 0.;
    dpPlanes[7].angleFromVert = TMath::Pi()/2.;
    dpPlanes[7].x0 = 0.747 - 50.;
    dpPlanes[7].y0 = 0.691 - 3.0083*2.54 - dpPlanes[7].nElements*dpPlanes[7].cellWidth/2. + 1.26992;
    dpPlanes[7].z0 = 1497.88;
    dpPlanes[7].x1 = dpPlanes[7].x0 - 50.;
    dpPlanes[7].y1 = dpPlanes[7].y0 - dpPlanes[7].nElements*dpPlanes[7].cellWidth/2.;
    dpPlanes[7].x2 = dpPlanes[7].x0 + 50.;
    dpPlanes[7].y2 = dpPlanes[7].y0 + dpPlanes[7].nElements*dpPlanes[7].cellWidth/2.;
    dpPlanes[7].thetaX = 0.015964;
    dpPlanes[7].thetaY = -0.000768;
    dpPlanes[7].thetaZ = -0.001297;
    dpPlanes[7].update();

    for(Int_t i = 0; i < 8; ++i)
    {
        cout << dpPlanes[i] << endl;
    }

    cout << dpPlanes[0].intercept(0.0419373, 0.00957994, -20.3776, 9.55376) << "  " << dpPlanes[0].wc << endl;
    cout << dpPlanes[4].intercept(0.0419373, 0.00957994, -20.3776, 9.55376) << "  " << dpPlanes[4].wc << endl;

    Int_t quad;
    Int_t fElementID;
    Int_t bElementID;
    Int_t fElementID_exp;
    Int_t bElementID_exp;
    Double_t fPos, bPos;
    Double_t fPos_exp, bPos_exp;
    Double_t tx, ty, x0, y0;

    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    saveTree->Branch("quad", &quad);
    saveTree->Branch("fElementID", &fElementID);
    saveTree->Branch("bElementID", &bElementID);
    saveTree->Branch("fElementID_exp", &fElementID_exp);
    saveTree->Branch("bElementID_exp", &bElementID_exp);
    saveTree->Branch("fPos", &fPos);
    saveTree->Branch("fPos_exp", &fPos_exp);
    saveTree->Branch("bPos", &bPos);
    saveTree->Branch("bPos_exp", &bPos_exp);
    saveTree->Branch("tx", &tx);
    saveTree->Branch("ty", &ty);
    saveTree->Branch("x0", &x0);
    saveTree->Branch("y0", &y0);

    for(Int_t i = 0; i < dataTree->GetEntries(); ++i)
    {
        dataTree->GetEntry(i);
        if(i % 1000 == 0) cout << i << endl;

        Int_t nTracks = tracklets->GetEntries();
        for(Int_t j = 0 ; j < nTracks; ++j)
        {
            Tracklet* track = (Tracklet*)tracklets->At(j);
            if(track->getNHits() < 16) continue;
            if(track->chisq > 10.) continue;

            tx = track->tx;
            ty = track->ty;
            x0 = track->x0;
            y0 = track->y0;

            //cout << j << "  " << tx << "  " << ty << "  " << x0 << "  " << y0 << endl;

            if(tx > 0.)
            {
                quad = ty > 0 ? 0 : 2;
            }
            else
            {
                quad = ty > 0 ? 1 : 3;
            }
            
            Int_t dID1 = quad;
            Int_t dID2 = quad + 4;

            //cout << "    " << quad << "  " << dID1 << "  " << dID2 << endl;

            fPos_exp = dpPlanes[dID1].intercept(tx, ty, x0, y0) - dpPlanes[dID1].wc;
            bPos_exp = dpPlanes[dID2].intercept(tx, ty, x0, y0) - dpPlanes[dID2].wc;
            if(quad == 0 || quad == 1)
            {
                fElementID_exp = Int_t(dpPlanes[dID1].nElements/2 + fPos_exp/dpPlanes[dID1].cellWidth) + 1;
                bElementID_exp = Int_t(dpPlanes[dID2].nElements/2 + bPos_exp/dpPlanes[dID2].cellWidth) + 1;
            }
            else
            {
                fElementID_exp = Int_t(dpPlanes[dID1].nElements/2 - fPos_exp/dpPlanes[dID1].cellWidth) + 1;
                bElementID_exp = Int_t(dpPlanes[dID2].nElements/2 - bPos_exp/dpPlanes[dID2].cellWidth) + 1;
            }

            if(fElementID_exp > 40) continue;

            //cout << "   " << fPos_exp << "  " << fElementID_exp << "  " << bPos_exp << "  " << bElementID_exp << endl;

            fPos = 9999.;
            bPos = 9999.;
            for(Int_t k = 0; k < rawEvent->getNHitsAll(); ++k)
            {
                Hit h = rawEvent->getHit(k);
                if(h.detectorID-55 != dID1 && h.detectorID-55 != dID2) continue;

                //cout << "    " << k << endl;
                //h.print();

                Double_t pos = (h.elementID - dpPlanes[h.detectorID - 55].nElements/2 + 0.5)*dpPlanes[h.detectorID - 55].cellWidth;
                if(quad == 2 || quad == 3) pos = -pos;

                if(h.detectorID-55 == dID1 && fabs(pos - fPos_exp) < fabs(fPos - fPos_exp))
                {
                    fElementID = h.elementID;
                    fPos = pos;
                }

                if(h.detectorID-55 == dID2 && fabs(pos - bPos_exp) < fabs(bPos - bPos_exp))
                {
                    bElementID = h.elementID;
                    bPos = pos;
                }

                //cout << pos << " " << fPos << "  " << bPos << "  " << fElementID << "  " << bElementID << endl;
            }

            saveTree->Fill();
        }

        rawEvent->clear();
    }

    saveFile->cd();
    saveTree->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}