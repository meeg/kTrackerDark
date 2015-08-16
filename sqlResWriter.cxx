#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TString.h>

#include "SRecEvent.h"
#include "JobOptsSvc.h"
#include "GeomSvc.h"
#include "MySQLSvc.h"

using namespace std;

int main(int argc, char **argv)
{
    cout << "Uploading file: " << argv[1] << " to sql schema " << argv[2] << endl;

    ///Initialize job option service
    JobOptsSvc* p_jobOptsSvc = JobOptsSvc::instance();

    ///Initialize the geometry service and output file
    GeomSvc* p_geomSvc = GeomSvc::instance();

    ///Intialize the mysql service
    MySQLSvc* p_mysqlSvc = MySQLSvc::instance();
    p_mysqlSvc->setUserPasswd("production", "qqbar2mu+mu-");
    p_mysqlSvc->connectOutput(argv[3], atoi(argv[4]));
    p_mysqlSvc->setOutputSchema(argv[2]);

    ///Retrieve data from file
    TClonesArray* tracklets = new TClonesArray("Tracklet");
    SRecEvent* recEvent = new SRecEvent();

    ///Get trees
    TFile* dataFile = new TFile(argv[1], "READ");
    TTree* configTree = (TTree*)dataFile->Get("config");
    TTree* dataTree = (TTree*)dataFile->Get("save");
    TTree* mixTree = (TTree*)dataFile->Get("save_mix");
    TTree* ppTree = (TTree*)dataFile->Get("save_pp");
    TTree* mmTree = (TTree*)dataFile->Get("save_mm");

    ///Table name and corresponding trees
    TString tableSuffix[4] = {"", "Mix", "PP", "MM"};
    TTree* trees[4] = {dataTree, mixTree, ppTree, mmTree};

    ///Initialize data tables
    if(!p_mysqlSvc->initWriter()) exit(EXIT_FAILURE);
    if(!p_mysqlSvc->initBakWriter()) exit(EXIT_FAILURE);

    ///Write the tracker configuration table
    p_mysqlSvc->writeInfoTable(configTree);

    ///Upload all tables
    for(int i = 0; i < 4; ++i)
    {
        trees[i]->SetBranchAddress("recEvent", &recEvent);

        TClonesArray* trackletArray = NULL;
        if(trees[i]->GetListOfBranches()->Contains("tracklets"))
        {
            trackletArray = tracklets;
            trees[i]->SetBranchAddress("tracklets", &trackletArray);
        }

        p_mysqlSvc->resetWriter();
        int nEvents = trees[i]->GetEntries();
        int printFreq = nEvents/100 > 1 ? nEvents/100 : 1;
        for(int j = 0; j < nEvents; ++j)
        {
            trees[i]->GetEntry(j);
            if(j % printFreq) cout << Form("Uploading %s event, %d percent finished.", tableSuffix[i].Data(), (j+1)*100/nEvents) << endl;

            p_mysqlSvc->writeTrackingRes(tableSuffix[i], recEvent, trackletArray);
        }
        cout << Form("Uploaded %s data events successfully", tableSuffix[i].Data()) << endl << endl;
    }
    cout << "sqlResWriter ends successfully." << endl;

    delete p_mysqlSvc;
    delete p_geomSvc;
    delete p_jobOptsSvc;

    return EXIT_SUCCESS;
}
