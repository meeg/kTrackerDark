kTracker
========

Author: Kun Liu, liuk@fnal.gov

Track reconstruction for E906/SeaQuest experiment

1. Installation

  Before installation, following environmental variables should be set/added:
  * Geant4 (4.9.6p02 or higher): $GEANTSYS
  * ROOT (5.34 or higher): $ROOTSYS
  * BOOST: $BOOST
  * MySQL_dev: $MYSQL_LIB and $MYSQL_INCLUDE
  * PATH: $GEANT4SYS/bin and $ROOTSYS/bin
  * LD_LIBRARY_PATH: $GEANT4SYS/lib and $ROOTSYS/lib
  
  On e906-gat2, there is an example environment setup script at ~liuk/env.sh. 
  
  Also one needs to modify the MODE_SWITCH.h:
  * KTRACKER_ROOT: change to the location of kTracker
  * MYSQL_SERVER: address of mysql server (contains geometry schema and/or data)
  * MYSQL_PORT: port number of mysql server
  * MC_MODE: enable if running MC productions
  * TRIGGER_TRIMING: uncomment to enable using only hodoscope hits that correspond to trigger roads for masking tracks
  * LOAD_ONLINE_ALIGNMENT: enable to use the alignment parameters stored in geometry schema, otherwise ascii files with                             alignment parameters need to be present at KTRACKER_ROOT
  * DIMUON_MODE: change to 0 if one wants to save events without muons pairs
  * KMAG_ON: change to 0 to turn off kMag
  
  Go to TrackExtrapolator and make first, then go back to KTRACKER_ROOT directory and make. Or run './reset.py'
  to do both with one command.
  
2. Executables

  After compilation, following executables should appear:
  * kFastTracking: E866-style track finding with Kalman filter/chi square track fitting
  * kOnlineTracking: track the data in real-time directly from MySQL while it is being decoded, 
                     results are stored in local ROOT file and also pushed back to the database. 
                     One thing to note is online tracking only use single muon vertex finding
  * kVertex: find the single muon/dimuon vertex via vertex fit with Kalman-fitted tracks
  
  If needed, there are several standalone executables that needs be compiled individually. Most likely those
  source files are located at KTRACKER_ROOT/analysis_tools, use './compile analysis_tools/executable_name' to compile:
  * sqlDataReader: reads the data from MySQL and save it in ROOT file
  * sqlMCReader: reads the MC data from MySQL and save it in ROOT file
  * update: update the wire position calucation with new alignment parameters

3. How to use
  
  1. Before running track finding/fitting, data needs to be extracted from MySQL to local ROOT files, 
     and/or apply latest alignment parameters
     * Read data: ./sqlDataReader run_name_in_mysql raw_data
     * Apply alignment parameters: ./update original_root_file updated_root_file (necessary only if there are new alignment parameters available)    
     
  2. With root file containing raw data, one can directly run fast tracking:
     * Fast tracking: ./kFastTracking raw_data raw_data_with_track
  
  3. Alternertively, one can also run online track reconstruction which directly read data from MySQL database
     * Online tracking: ./kOnlineTracking run_name_in_mysql raw_data_with_track
  
  4. After tracks are found, one can run both single muon/dimuon vertex finding to calculate Minv, etc.
     * Vertex finding: ./kVertex raw_data_with_track raw_data_with_vertex
