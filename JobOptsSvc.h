/*! @brief Read a job options file and apply choices at runtime.



  @author Brian G. Tice
  */

#ifndef _JOBOPTSSVC_H
#define _JOBOPTSSVC_H

#include <string>

class JobOptsSvc
{
public:
  ///singlton instance
  static JobOptsSvc* instance();

  ///Initialization with defaults
  bool init( bool forceInit = false );
  ///Initialization using this config file
  bool init(const char* filename);
  ///Close the service and cleanup
  void close();


  ///Return a string with environmental variables expanded
  std::string ExpandEnvironmentals( const std::string& input ) const;

  ///Return the location of PlusTop roads file
  std::string GetRoadsFilePlusTop() const;
  ///Return the location of PlusBottom roads file
  std::string GetRoadsFilePlusBottom() const;
  ///Return the location of MinusTop roads file
  std::string GetRoadsFileMinusTop() const;
  ///Return the location of MinusBottom roads file
  std::string GetRoadsFileMinusBottom() const;

  //@todo should store smart pointers instead of variable length variables

  std::string m_configFile; ///< Name of the config file loaded

  bool m_mcMode;           ///< Running on MC?
  bool m_alignmentMode;    ///< Running in alignment mode?
  bool m_enableTriggerMask;///< Enable hodo masking with trigger road info
  bool m_enableKMag;       ///< Turn kMag on
  bool m_enableOnlineAlignment;   ///< Turn kMag on
  bool m_enableEvaluation; ///< Enable evaluation output
  bool m_thresholdLive;    ///< Live setting for ThresholdSvc

  int m_mySQLPort;    ///< mysql database port
  int m_nEvents;      ///< number of events to process
  int m_firstEvent;   ///< first event to process

  int m_triggerL1;    ///< L1 firmware version (trigger roads)

  int m_thresholdLevel; ///< Level setting for ThresholdSvc

  std::string m_inputFile;  ///< Name of the input file
  std::string m_outputFile; ///< Name of the output file

  std::string m_alignmentFileHodo; ///< Name of hodoscope alignment file
  std::string m_alignmentFileChamber; ///< Name of chamber alignment file
  std::string m_alignmentFileProp; ///< Name of prop tune alignment file
  std::string m_alignmentFileMille; ///< Name of mille alignment file

  std::string m_triggerRepo;  ///< Location of the trigger repository on local disk

  std::string m_calibrationsFile; ///< Name of calibrations file

  std::string m_fMagFile; ///< Name of fMag ascii file
  std::string m_kMagFile; ///< Name of kMag ascii file

  std::string m_geomVersion; ///< Name of geometry version
  std::string m_mySQLServer;  ///< Name of mysql database

  std::string m_mySQLurl;  ///< url of MySQL 

private:
  int m_isInit; ///< Has this service been initialized?
  static JobOptsSvc* p_jobOptsSvc; ///< Singleton pointer
};

#endif
