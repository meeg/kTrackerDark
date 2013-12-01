kTracker
========

Author: Kun Liu, liuk@fnal.gov

Track reconstruction for E906/SeaQuest experiment

1. Installation

  Before installation, following environmental variables should be set/added:
  * Geant4 (4.9.6p02 or higher): $GEANTSYS
  * ROOT (5.34 or higher): $ROOTSYS
  * MySQL_dev: $MYSQL_LIB and $MYSQL_INCLUDE
  * PATH: $GEANT4SYS/bin and $ROOTSYS/bin
  * LD_LIBRARY_PATH: $GEANT4SYS/lib and $ROOTSYS/lib
  
  On e906-gat2, there is an example environment setup script at ~liuk/env.sh.
  
  Also one needs to modify the MODE_SWITCH.h:
  * KTRACKER_ROOT: change to the location of kTracker
  * MYSQL_SERVER: address of mysql server (contains geometry schema and/or data)
  * ALIGNMENT_MODE: enable if running the alignment (it will set both FMag and KMag to 0.)
  * MC_MODE: enable if running MC productions
  * DIMUON_MODE: change to 0 if one wants to save events with less than 2 muons
  * KMAG_ON: change to 0 to turn off kMag
  * ENABLE_KF: comment out if one wants to disable Kalman Filter
  * Other variables
  
  Go to TrackExtrapolator and make first, then go back to KTRACKER_ROOT directory and make. Or run './reset.py'
  to do both with one command.
  
2. Executables

  After compilation, following executables should appear:
  * kSeeder: find track seeds using prop tube and station 3&4 hodo
  * kTracker: find/fit tracks using track seeds produced by kSeeder (slow)
  * kFastTracking: E866-style track finding with Kalman filter/chi square track fitting
  * kOnlineTracking: track the data in real-time directly from MySQL while it is being decoded, 
                     results are stored in local ROOT file and also pushed back to the database. 
                     One thing to note is online tracking only use single muon vertex finding
  * milleAlign: millepede-based detector alignment
  * kVertex: find the single muon/dimuon vertex via vertex fit with Kalman-fitted tracks
  * kVertex_mix: run vertex fit upon randomly combined tracks
  
  If needed, there are several standalone executables that needs be compiled individually. Most likely those
  source files are located at KTRACKER_ROOT/analysis_tools, use './compile analysis_tools/executable_name' to compile:
  * sqlDataReader: reads the data from MySQL and save it in ROOT file
  * sqlMCReader: reads the MC data from MySQL and save it in ROOT file
  * update: update the wire position calucation with new alignment parameters
  * makeRTProfile: produce the R-T profile for drift chambers based on tracking results
  * hodoAlign: find the alignment parameters (shift only) for hodoscopes
  * propAlign: find the alignment parameters (shift only) for prop. tubes

  There are some extra python scripts to run repeated works:
  * runAlignment.py: run the alignment iteratively

3. How to use
  
  1. Before running track finding/fitting, data needs to be extracted from MySQL to local ROOT files, 
     and/or apply latest alignment parameters
     * Read data: ./sqlDataReader run_name_in_mysql raw_data
     * Apply alignment parameters: ./update original_root_file updated_root_file
     
  2. Run complete Kalman filter based track finding and fitting:
     * Find seeds: ./kSeeder raw_data raw_data_with_seed
     * Find tracks: ./kTracker raw_data_with_seed raw_data_with_track
     
  3. Alternatively, one can also directly run fast tracking:
     * Fast tracking: ./kFastTracking raw_data raw_data_with_track
  
  4. One can also run online track reconstruction which directly read data from MySQL database
     * Online tracking: ./kOnlineTracking run_name_in_mysql raw_data_with_track
  
  5. After tracks are found, one can run both single muon/dimuon vertex finding to calculate Minv, etc.
     * With Kalman-fitted tracks: ./kVertex raw_data_with_track raw_data_with_vertex
    
4. Table of naming scheme of hodos and mapping to detectorID
   
   H1B    H1X1  25
   H1T    H1X2  26
   H1L    H1Y1  27
   H1R    H1Y2  28
   H2L    H2Y1  29
   H2R    H2Y2  30
   H2B    H2X1  31
   H2T    H2X2  32
   H3T    H3X1  33
   H3B    H3X2  34
   H4Y1L  H4Y11 35
   H4Y1R  H4Y12 36
   H4Y2L  H4Y21 37
   H4T2R  H4Y22 38
   H4B    H4X1  39
   H4T    H4X2  40
  
   
     
