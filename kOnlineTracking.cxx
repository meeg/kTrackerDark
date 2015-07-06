#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <time.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "ThresholdSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"
#include "KalmanFastTracking.h"
#include "KalmanFitter.h"
#include "VertexFit.h"
#include "MySQLSvc.h"
#include "JobOptsSvc.h"
#include "EventReducer.h"
#include "MODE_SWITCH.h"

using namespace std;
using Threshold::live;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout << "Usage: " << argv[0] << "  <options file>" << endl;
        exit(0);
    }

    //Initialize job options
    JobOptsSvc* jobOptsSvc = JobOptsSvc::instance();
    jobOptsSvc->init(argv[1]);

    //Initialize geometry service with calibration
    GeomSvc::instance()->init();

    //Initialize MySQL service and connect to database, e906-db1 by default
    MySQLSvc* p_mysqlSvc = MySQLSvc::instance();
    p_mysqlSvc->setUserPasswd("production", "qqbar2mu+mu-");
    p_mysqlSvc->connectInput();
    p_mysqlSvc->connectOutput();
    if(!(p_mysqlSvc->initReader() && p_mysqlSvc->initWriter())) exit(EXIT_FAILURE);

    //Data output definition
    int nTracklets;
    SRawEvent* rawEvent = new SRawEvent();
    SRecEvent* recEvent = new SRecEvent();
    TClonesArray* tracklets = new TClonesArray("Tracklet");
    TClonesArray& arr_tracklets = *tracklets;

    //Initialize track finder
#ifdef _ENABLE_KF
    KalmanFilter* filter = new KalmanFilter();
    KalmanFastTracking* fastfinder = new KalmanFastTracking();
#else
    KalmanFastTracking* fastfinder = new KalmanFastTracking(false);
#endif
    VertexFit* vtxfit  = new VertexFit();

    //Initialize event reducer
    TString opt = "aoc";      //turn on after pulse removal, out of time removal, and cluster removal
    if(jobOptsSvc->m_enableTriggerMask) opt = opt + "t";
    if(jobOptsSvc->m_sagittaReducer) opt = opt + "s";
    if(jobOptsSvc->m_updateAlignment) opt = opt + "e";
    EventReducer* eventReducer = new EventReducer(opt);

    //Quality control numbers and plots
    int nEvents_loaded = 0;
    int nEvents_tracked = 0;
    int nEvents_dimuon = 0;
    int nEvents_dimuon_real = 0;

    //the number of events available in DB within requested eventID range (firstEvent through firstEvent+nEvents)
    const int nEvents = p_mysqlSvc->getNEvents();
    cout << Form( "Found %d events in DB in eventID range %d-%d.", nEvents, jobOptsSvc->m_firstEvent, jobOptsSvc->m_firstEvent+jobOptsSvc->m_nEvents );

    //Start tracking
    if(!jobOptsSvc->m_mcMode) p_mysqlSvc->getOutputServer()->Exec(Form("UPDATE summary.production SET kTrackBegin=NOW() WHERE production='%s'", jobOptsSvc->m_outputSchema.c_str()));
    for(int i = 0; i < nEvents; ++i)
    {
        //Read data
        if(!p_mysqlSvc->getNextEvent(rawEvent)) continue;
        ++nEvents_loaded;

        //Do the tracking
        if(live())
        {
            cout << "\r Tracking runID = " << rawEvent->getRunID() << " eventID = " << rawEvent->getEventID() << ", " << (i+1)*100/nEvents << "% finished. ";
            cout << nEvents_tracked*100/nEvents_loaded << "% have at least one track, " << nEvents_dimuon*100/nEvents_loaded << "% have at least one dimuon pair, ";
            cout << nEvents_dimuon_real*100/nEvents_loaded << "% have successful dimuon vertex fit.";
        }

        eventReducer->reduceEvent(rawEvent);
        if(!fastfinder->setRawEvent(rawEvent)) continue;
        ++nEvents_tracked;

        //Output
        arr_tracklets.Clear();
        std::list<Tracklet>& rec_tracklets = fastfinder->getFinalTracklets();
        if(rec_tracklets.empty()) continue;

        recEvent->setRawEvent(rawEvent);
        nTracklets = 0;
        int nPos = 0;
        int nNeg = 0;
        for(std::list<Tracklet>::iterator iter = rec_tracklets.begin(); iter != rec_tracklets.end(); ++iter)
        {
            //iter->print();
            iter->calcChisq();
            new(arr_tracklets[nTracklets++]) Tracklet(*iter);

            SRecTrack recTrack = iter->getSRecTrack();
            recEvent->insertTrack(recTrack);

            iter->getCharge() > 0 ? ++nPos : ++nNeg;
        }
        if(nPos > 0 && nNeg > 0) ++nEvents_dimuon;

        //Perform dimuon vertex fit
        if(vtxfit->setRecEvent(recEvent)) ++nEvents_dimuon_real;

        if(nTracklets > 0)
        {
            p_mysqlSvc->writeTrackingRes(recEvent, tracklets);
        }
        rawEvent->clear();
        recEvent->clear();
    }
    //End tracking

    //push intermediate event subset tables to final tables (if needed)
    bool dropSubsetTables = true; //change to false if there was a processing problem
    p_mysqlSvc->pushToFinalTables( dropSubsetTables );

    //update summary table
    if(!jobOptsSvc->m_mcMode)
    {
        p_mysqlSvc->getOutputServer()->Exec(Form("UPDATE summary.production SET ktracked=1 WHERE production='%s'", jobOptsSvc->m_outputSchema.c_str()));
        p_mysqlSvc->getOutputServer()->Exec(Form("UPDATE summary.production SET kTrackEnd=NOW() WHERE production='%s'", jobOptsSvc->m_outputSchema.c_str()));
    }

    cout << endl;
    cout << "kOnlineTracking ended successfully." << endl;
    cout << "In total " << nEvents_loaded << " events loaded from " << argv[1] << ": " << nEvents_tracked << " events have at least one track, ";
    cout << nEvents_dimuon << " events have at least one dimuon pair, ";
    cout << nEvents_dimuon_real << " events have successful dimuon vertex fit." << endl;

    delete fastfinder;
    delete vtxfit;
    delete eventReducer;
#ifdef _ENABLE_KF
    filter->close();
#endif

    return EXIT_SUCCESS;
}
