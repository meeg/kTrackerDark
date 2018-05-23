#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"

using namespace std;

double getDelay(int quadID, int barID) {
    switch(quadID) {
        case 0:
            return 1728.0;
        case 1:
            return 1720.0;
        case 2:
            return 1730.0;
        case 3:
            return (barID>15)?1716.0:1721.0;
        case 4:
            return (barID>32)?1695.0:1703.0;
        case 5:
            return 1697.0;
        case 6:
            return 1704.0;
        case 7:
            return (barID>32)?1695.0:1689.0;
        case 8:
            return 1092;
        case 9:
            return 1091;
        case 10:
            return 1093;
        case 11:
            return 1092;
        default:
            return 1700;
    }
}
int roadHash(int st1, int st2, int h4) {
    return st1*1000 + st2*10 + h4;
}
int main(int argc, char *argv[])
{


    unordered_set<int> roadset;
    ifstream roadsfile("code_400_real_roads.txt");
    if (!roadsfile.is_open()) return 1;
    char line[100];
    while (roadsfile.getline(line,100)) {
        //printf("line=%s\n",line);
        char * tok;
        int st1, st2, h4;
        st1 = atoi(strtok(line," "));
        st2 = atoi(strtok(NULL," "));
        h4 = atoi(strtok(NULL," "));
        int roadhash = roadHash(st1+1, st2+1, h4+1); //VHDL uses zero-indexed vectors, so we need to add 1 to match the elementIDs
        //printf("road=%d\n",roadhash);
        roadset.insert(roadhash);
    }


    //Initialization of geometry and tracked data
    //GeomSvc* p_geomSvc = GeomSvc::instance();
    //p_geomSvc->init();

    SRawEvent* rawEvent = new SRawEvent();

    TFile* dataFile = new TFile(argv[1], "READ");
    TTree* dataTree = (TTree*)dataFile->Get("save");

    dataTree->SetBranchAddress("rawEvent", &rawEvent);
    //dataTree->SetBranchAddress("tracklets", &tracklets);

    set<int> hitsets[12];
    set<int> hitsetsNoTimecut[12];
    set<int> hitsetsOffByOne[12];

    TFile* saveFile = new TFile(argv[2], "recreate");
    //TTree* saveTree = new TTree("save", "save");

    saveFile->cd();
    TH2D* fbCorrHists[4];
    TH2D* fbCorrHistsNIM1[4];
    TH2D* fbCorrHistsNIM3[4];
    TH2D* bhCorrHists[4];
    TH2D* bhCorrHistsNIM1[4];
    TH2D* bhCorrHistsNIM3[4];
    TH2D* fbCorrHistsOffByOne[4];
    TH2D* fbCorrHistsOffByOneNIM1[4];
    TH2D* fbCorrHistsOffByOneNIM3[4];
    TH2D* bhCorrHistsOffByOne[4];
    TH2D* bhCorrHistsOffByOneNIM1[4];
    TH2D* bhCorrHistsOffByOneNIM3[4];
    char name[500];
    for (int quad=0; quad<4; quad++) {
        sprintf(name, "fbCorrHist_%i", quad);
        fbCorrHists[quad] = new TH2D(name,name,30,0.5,30.5,50,0.5,50.5);
        sprintf(name, "fbCorrHistNIM1_%i", quad);
        fbCorrHistsNIM1[quad] = new TH2D(name,name,30,0.5,30.5,50,0.5,50.5);
        sprintf(name, "fbCorrHistNIM3_%i", quad);
        fbCorrHistsNIM3[quad] = new TH2D(name,name,30,0.5,30.5,50,0.5,50.5);
        sprintf(name, "bhCorrHist_%i", quad);
        bhCorrHists[quad] = new TH2D(name,name,8,0.5,8.5,50,0.5,50.5);
        sprintf(name, "bhCorrHistNIM1_%i", quad);
        bhCorrHistsNIM1[quad] = new TH2D(name,name,8,0.5,8.5,50,0.5,50.5);
        sprintf(name, "bhCorrHistNIM3_%i", quad);
        bhCorrHistsNIM3[quad] = new TH2D(name,name,8,0.5,8.5,50,0.5,50.5);

        sprintf(name, "fbCorrHistOffByOne_%i", quad);
        fbCorrHistsOffByOne[quad] = new TH2D(name,name,30,0.5,30.5,50,0.5,50.5);
        sprintf(name, "fbCorrHistOffByOneNIM1_%i", quad);
        fbCorrHistsOffByOneNIM1[quad] = new TH2D(name,name,30,0.5,30.5,50,0.5,50.5);
        sprintf(name, "fbCorrHistOffByOneNIM3_%i", quad);
        fbCorrHistsOffByOneNIM3[quad] = new TH2D(name,name,30,0.5,30.5,50,0.5,50.5);
        sprintf(name, "bhCorrHistOffByOne_%i", quad);
        bhCorrHistsOffByOne[quad] = new TH2D(name,name,8,0.5,8.5,50,0.5,50.5);
        sprintf(name, "bhCorrHistOffByOneNIM1_%i", quad);
        bhCorrHistsOffByOneNIM1[quad] = new TH2D(name,name,8,0.5,8.5,50,0.5,50.5);
        sprintf(name, "bhCorrHistOffByOneNIM3_%i", quad);
        bhCorrHistsOffByOneNIM3[quad] = new TH2D(name,name,8,0.5,8.5,50,0.5,50.5);
    }
    double rfclock = 18.8;
    double tmin=-8.5*rfclock;
    double tmax=7.5*rfclock;
    int nbins=256;
    TH2D* hitTimeHists[12];
    TH2D* hitTimeInTimeHists[12];
    TH2D* hitTimeHistsNIM1[12];
    TH2D* hitTimeHistsNIM3[12];
    TH1D* hitEleHists[12];
    TH1D* hitEleHistsNIM1[12];
    TH1D* hitEleHistsNIM3[12];
    for (int quad=0; quad<12; quad++) {
        sprintf(name, "hitTimeHist_%i", quad);
        if (quad<4) {
            hitTimeHists[quad] = new TH2D(name,name,nbins,tmin,tmax,30,0.5,30.5);
        } else if (quad<8) {
            hitTimeHists[quad] = new TH2D(name,name,nbins,tmin,tmax,50,0.5,50.5);
        } else {
            hitTimeHists[quad] = new TH2D(name,name,nbins,tmin,tmax,8,0.5,8.5);
        }
        sprintf(name, "hitTimeInTimeHist_%i", quad);
        if (quad<4) {
            hitTimeInTimeHists[quad] = new TH2D(name,name,nbins,tmin,tmax,30,0.5,30.5);
        } else if (quad<8) {
            hitTimeInTimeHists[quad] = new TH2D(name,name,nbins,tmin,tmax,50,0.5,50.5);
        } else {
            hitTimeInTimeHists[quad] = new TH2D(name,name,nbins,tmin,tmax,8,0.5,8.5);
        }
        sprintf(name, "hitTimeHistNIM1_%i", quad);
        if (quad<4) {
            hitTimeHistsNIM1[quad] = new TH2D(name,name,nbins,tmin,tmax,30,0.5,30.5);
        } else if (quad<8) {
            hitTimeHistsNIM1[quad] = new TH2D(name,name,nbins,tmin,tmax,50,0.5,50.5);
        } else {
            hitTimeHistsNIM1[quad] = new TH2D(name,name,nbins,tmin,tmax,8,0.5,8.5);
        }
        sprintf(name, "hitTimeHistNIM3_%i", quad);
        if (quad<4) {
            hitTimeHistsNIM3[quad] = new TH2D(name,name,nbins,tmin,tmax,30,0.5,30.5);
        } else if (quad<8) {
            hitTimeHistsNIM3[quad] = new TH2D(name,name,nbins,tmin,tmax,50,0.5,50.5);
        } else {
            hitTimeHistsNIM3[quad] = new TH2D(name,name,nbins,tmin,tmax,8,0.5,8.5);
        }
        sprintf(name, "hitEleHist_%i", quad);
        if (quad<4) {
            hitEleHists[quad] = new TH1D(name,name,30,0.5,30.5);
        } else if (quad<8) {
            hitEleHists[quad] = new TH1D(name,name,50,0.5,50.5);
        } else {
            hitEleHists[quad] = new TH1D(name,name,8,0.5,8.5);
        }
        sprintf(name, "hitEleHistNIM1_%i", quad);
        if (quad<4) {
            hitEleHistsNIM1[quad] = new TH1D(name,name,30,0.5,30.5);
        } else if (quad<8) {
            hitEleHistsNIM1[quad] = new TH1D(name,name,50,0.5,50.5);
        } else {
            hitEleHistsNIM1[quad] = new TH1D(name,name,8,0.5,8.5);
        }
        sprintf(name, "hitEleHistNIM3_%i", quad);
        if (quad<4) {
            hitEleHistsNIM3[quad] = new TH1D(name,name,30,0.5,30.5);
        } else if (quad<8) {
            hitEleHistsNIM3[quad] = new TH1D(name,name,50,0.5,50.5);
        } else {
            hitEleHistsNIM3[quad] = new TH1D(name,name,8,0.5,8.5);
        }
    }
    TH2D* fbQuadHist = new TH2D("fbQuadHist","fbQuadHist",4,-0.5,3.5,4,-0.5,3.5);
    TH2D* fbQuadHistNIM1 = new TH2D("fbQuadHistNIM1","fbQuadHistNIM1",4,-0.5,3.5,4,-0.5,3.5);
    TH2D* fbQuadHistNIM3 = new TH2D("fbQuadHistNIM3","fbQuadHistNIM3",4,-0.5,3.5,4,-0.5,3.5);
    TH2D* bhQuadHist = new TH2D("bhQuadHist","bhQuadHist",4,-0.5,3.5,4,-0.5,3.5);
    TH2D* bhQuadHistNIM1 = new TH2D("bhQuadHistNIM1","bhQuadHistNIM1",4,-0.5,3.5,4,-0.5,3.5);
    TH2D* bhQuadHistNIM3 = new TH2D("bhQuadHistNIM3","bhQuadHistNIM3",4,-0.5,3.5,4,-0.5,3.5);


    TH2D* roadsVsTrig = new TH2D("roadsVsTrig","roadsVsTrig",200,-0.5,199.5,16,-0.5,15.5);
    TH2D* roadsVsTrigNoTimecut = new TH2D("roadsVsTrigNoTimecut","roadsVsTrigNoTimecut",200,-0.5,199.5,16,-0.5,15.5);
    TH2D* nimlikeVsTrig = new TH2D("nimlikeVsTrig","nimlikeVsTrig",200,-0.5,199.5,16,-0.5,15.5);

    double timeCut = rfclock/2;
    int dpTriggerMask = 64;
    int nim1TriggerMask = 32;
    int nim3TriggerMask = 128;
    for(Int_t i = 0; i < dataTree->GetEntries(); ++i)
    {
        dataTree->GetEntry(i);
        if(i % 1000 == 0) cout << i << endl;

        //if (rawEvent->getTriggerBits()>0 && (rawEvent->getTriggerBits() & (dpTriggerMask|nim1TriggerMask|nim3TriggerMask)) != 0) {
        if (rawEvent->getTriggerBits()>0) {
            //if (rawEvent->getTriggerBits()==64) cout << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << rawEvent->getTriggerBits() << endl;
            for(Int_t k = 0; k < rawEvent->getNHitsAll(); ++k) {
                Hit h = rawEvent->getHit(k);
                int dpQuadID = -1;
                int barID = -1;
                //if (h.detectorID >= 41 && h.detectorID < 45) dpQuadID = h.detectorID-33; //H4Y1L = 8, H4Y1R = 9, etc.
                if (h.detectorID >= 55 && h.detectorID <= 62) {
                    dpQuadID = h.detectorID-55;
                    barID = h.elementID;
                }
                else if (h.detectorID == 43) { // H4Y2L
                    //else if (h.detectorID == 41) { // H4Y2L
                    if (h.elementID > 8) {
                        dpQuadID = 8;
                        barID = h.elementID-8;
                    } else {
                        dpQuadID = 10;
                        barID = 9-h.elementID;
                    }
                }
                    else if (h.detectorID == 44) { // H4Y2R
                        //else if (h.detectorID == 42) { // H4Y2R
                        if (h.elementID > 8) {
                            dpQuadID = 9;
                            barID = h.elementID-8;
                        } else {
                            dpQuadID = 11;
                            barID = 9-h.elementID;
                        }
                    }

                    if(dpQuadID<0 || dpQuadID>11) continue;
                    double deltaT = h.tdcTime - getDelay(dpQuadID,barID);

                    //tweak deltaT so the NIM triggers are in time (the delays are tuned for the DP trigger)
                    if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                        deltaT -= 5.0;
                    } else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                        deltaT -= 15.0;
                    }

                    //if (h.isInTime()) {
                    hitsetsNoTimecut[dpQuadID].insert(barID);
                    if (TMath::Abs(deltaT)<timeCut) hitsets[dpQuadID].insert(barID);
                    if (TMath::Abs(deltaT+rfclock)<timeCut) hitsetsOffByOne[dpQuadID].insert(barID);
                    //if (rawEvent->getRunID()==28574 && rawEvent->getEventID() == 10816)
                    //if (rawEvent->getTriggerBits()==64 && TMath::Abs(deltaT)<timeCut) cout << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << rawEvent->getTriggerBits() << " " << dpQuadID << " " << barID << " " << h.tdcTime << " " << deltaT << endl;

                    if (rawEvent->getTriggerBits() == dpTriggerMask) {
                        hitTimeHists[dpQuadID]->Fill(deltaT,barID);
                        //if (isInTime(dpQuadID,barID,h.tdcTime)) {
                        if (h.isInTime()) {
                            hitTimeInTimeHists[dpQuadID]->Fill(deltaT,barID);
                            //hitEleHists[dpQuadID]->Fill(barID);
                        }
                    }
                    if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                        hitTimeHistsNIM1[dpQuadID]->Fill(deltaT,barID);
                        //if (h.isInTime()) hitEleHistsNIM1[dpQuadID]->Fill(barID);
                    }
                    if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                        hitTimeHistsNIM3[dpQuadID]->Fill(deltaT,barID);
                        //if (h.isInTime()) hitEleHistsNIM3[dpQuadID]->Fill(barID);
                    }
                    }

                    int roadBits = 0;
                    int roadBitsNoTimecut = 0;
                    int nimlikeBits = 0;
                    for (int quad=0;quad<4;quad++) {
                        int fID = (quad)%4;
                        int bID = quad+4;
                        int hID = quad+8;

                        for (set<int>::iterator fIt = hitsets[fID].begin(); fIt!=hitsets[fID].end();++fIt) {
                            for (set<int>::iterator bIt = hitsets[bID].begin(); bIt!=hitsets[bID].end();++bIt) {
                                for (set<int>::iterator hIt = hitsets[hID].begin(); hIt!=hitsets[hID].end();++hIt) {
                                    int roadhash = roadHash(*fIt,*bIt,*hIt);
                                    if (roadset.count(roadhash)) {
                                        //cout << "road match: " << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << rawEvent->getTriggerBits() << " " << quad << " " << roadhash << endl;
                                        roadBits |= (1 << quad);
                                    }
                                }
                            }
                        }
                        for (set<int>::iterator fIt = hitsetsNoTimecut[fID].begin(); fIt!=hitsetsNoTimecut[fID].end();++fIt) {
                            for (set<int>::iterator bIt = hitsetsNoTimecut[bID].begin(); bIt!=hitsetsNoTimecut[bID].end();++bIt) {
                                for (set<int>::iterator hIt = hitsetsNoTimecut[hID].begin(); hIt!=hitsetsNoTimecut[hID].end();++hIt) {
                                    int roadhash = roadHash(*fIt,*bIt,*hIt);
                                    if (roadset.count(roadhash)) {
                                        //cout << "road match: " << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << rawEvent->getTriggerBits() << " " << quad << " " << roadhash << endl;
                                        roadBitsNoTimecut |= (1 << quad);
                                    }
                                }
                            }
                        }
                        //if (!hitsets[fID].empty() && !hitsets[bID].empty() && !hitsets[hID].empty()) nimlikeBits |= (1 << quad);
                        //if (!hitsets[fID].empty() && !hitsets[bID].empty()) nimlikeBits |= (1 << quad);
                        //if (!hitsets[fID].empty()) nimlikeBits |= (1 << quad);
                        if (!hitsets[bID].empty()) nimlikeBits |= (1 << quad);
                        //if (!hitsets[hID].empty()) nimlikeBits |= (1 << quad);
                        //if (!hitsetsOffByOne[hID].empty()) nimlikeBits |= (1 << quad);
                        //if ((!hitsets[fID].empty() || !hitsets[bID].empty()) && !hitsets[hID].empty()) nimlikeBits |= (1 << quad);
                    }
                    roadsVsTrig->Fill(rawEvent->getTriggerBits(),roadBits);
                    //if (roadBits) cout << "roadbits: " << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << rawEvent->getTriggerBits() << " " << roadBits << endl;
                    roadsVsTrigNoTimecut->Fill(rawEvent->getTriggerBits(),roadBitsNoTimecut);
                    nimlikeVsTrig->Fill(rawEvent->getTriggerBits(),nimlikeBits);

                    int onlyquadF = -1;
                    int onlyquadB = -1;
                    int onlyquadH = -1;
                    int nhitsF, nhitsB, nhitsH;
                    int nquadsF = 0;
                    int nquadsB = 0;
                    int nquadsH = 0;
                    for (int quad=0;quad<4;quad++) {
                        int fID = (quad)%4;
                        int bID = quad+4;
                        int hID = quad+8;
                        int fHits = hitsets[fID].size();
                        int bHits = hitsets[bID].size();
                        int hHits = hitsets[hID].size();
                        if (fHits>0) {
                            nquadsF++;
                            if (onlyquadF==-1) {
                                onlyquadF = quad;
                                nhitsF = fHits;
                            } else
                                onlyquadF = -2;
                        }
                        if (bHits>0) {
                            nquadsB++;
                            if (onlyquadB==-1) {
                                onlyquadB = quad;
                                nhitsB = bHits;
                            } else
                                onlyquadB = -2;
                        }
                        if (hHits>0) {
                            nquadsH++;
                            if (onlyquadH==-1) {
                                onlyquadH = quad;
                                nhitsH = hHits;
                            } else
                                onlyquadH = -2;
                        }
                    }

                    for (int quad=0;quad<12;quad++) {
                        for (set<int>::iterator it = hitsets[quad].begin(); it!=hitsets[quad].end();++it) {
                            if (rawEvent->getTriggerBits() == dpTriggerMask) {
                                hitEleHists[quad]->Fill(*it);
                            }
                            else if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                                hitEleHistsNIM1[quad]->Fill(*it);
                            }
                            else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                                hitEleHistsNIM3[quad]->Fill(*it);
                            }
                        }
                    }

                    for (int quad=0;quad<4;quad++) {
                        int fID = (quad)%4;
                        int bID = quad+4;
                        int hID = quad+8;
                        //if ((rawEvent->getTriggerBits() & 64) != 0) {
                        //cout << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << quad << " " << fHits << " " << bHits << endl;
                        //}
                        //if (fHits<2 && bHits<2)
                        //if (!hitsets[((quad+1)%4)+8].empty())
                        //if (nquadsB==2 && nquadsH==2)
                        for (set<int>::iterator bIt = hitsets[bID].begin(); bIt!=hitsets[bID].end();++bIt) {
                            for (set<int>::iterator fIt = hitsets[fID].begin(); fIt!=hitsets[fID].end();++fIt) {
                                if (rawEvent->getTriggerBits() == dpTriggerMask) {
                                    fbCorrHists[quad]->Fill(*fIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                                    fbCorrHistsNIM1[quad]->Fill(*fIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                                    fbCorrHistsNIM3[quad]->Fill(*fIt,*bIt);
                                }
                            }
                            for (set<int>::iterator hIt = hitsets[hID].begin(); hIt!=hitsets[hID].end();++hIt) {
                                if (rawEvent->getTriggerBits() == dpTriggerMask) {
                                    bhCorrHists[quad]->Fill(*hIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                                    bhCorrHistsNIM1[quad]->Fill(*hIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                                    bhCorrHistsNIM3[quad]->Fill(*hIt,*bIt);
                                }
                            }
                        }

                        for (set<int>::iterator bIt = hitsetsOffByOne[bID].begin(); bIt!=hitsetsOffByOne[bID].end();++bIt) {
                            for (set<int>::iterator fIt = hitsetsOffByOne[fID].begin(); fIt!=hitsetsOffByOne[fID].end();++fIt) {
                                if (rawEvent->getTriggerBits() == dpTriggerMask) {
                                    fbCorrHistsOffByOne[quad]->Fill(*fIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                                    fbCorrHistsOffByOneNIM1[quad]->Fill(*fIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                                    fbCorrHistsOffByOneNIM3[quad]->Fill(*fIt,*bIt);
                                }
                            }
                            for (set<int>::iterator hIt = hitsetsOffByOne[hID].begin(); hIt!=hitsetsOffByOne[hID].end();++hIt) {
                                if (rawEvent->getTriggerBits() == dpTriggerMask) {
                                    bhCorrHistsOffByOne[quad]->Fill(*hIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                                    bhCorrHistsOffByOneNIM1[quad]->Fill(*hIt,*bIt);
                                }
                                else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                                    bhCorrHistsOffByOneNIM3[quad]->Fill(*hIt,*bIt);
                                }
                            }
                        }

                    }

                    if (onlyquadF>=0 && onlyquadB >= 0) {
                        if (rawEvent->getTriggerBits() == dpTriggerMask) {
                            fbQuadHist->Fill(onlyquadF, onlyquadB);
                        }
                        else if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                            fbQuadHistNIM1->Fill(onlyquadF, onlyquadB);
                        }
                        else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                            fbQuadHistNIM3->Fill(onlyquadF, onlyquadB);
                        }
                        //cout << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << onlyquadF << " " << onlyquadB << " " << nhitsF << " " << nhitsB << endl;
                    }
                    if (onlyquadH>=0 && onlyquadB >= 0) {
                        if (rawEvent->getTriggerBits() == dpTriggerMask) {
                            bhQuadHist->Fill(onlyquadH, onlyquadB);
                        }
                        else if (rawEvent->getTriggerBits() == nim1TriggerMask) {
                            bhQuadHistNIM1->Fill(onlyquadH, onlyquadB);
                        }
                        else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
                            bhQuadHistNIM3->Fill(onlyquadH, onlyquadB);
                        }
                        //cout << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << onlyquadF << " " << onlyquadB << " " << nhitsF << " " << nhitsB << endl;
                    }

                    //saveTree->Fill();
                    /*
                       cout << rawEvent->getRunID() << " " << rawEvent->getEventID() << " " << rawEvent->getTriggerBits();
                       for (int j=0;j<12;j++) {
                       cout << " " << hitsets[j].size();
                       }
                       cout << endl;
                       */

                    for (int j=0;j<12;j++) {
                        hitsets[j].clear();
                        hitsetsNoTimecut[j].clear();
                        hitsetsOffByOne[j].clear();
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
