#include "JobOptsSvc.h"

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <wordexp.h> //to expand environmentals
#include <boost/algorithm/string.hpp> //for strip
#include <boost/lexical_cast.hpp>

using namespace std;

//todo a msg service
//#define DEBUG_JobOptsSvc  //comment out for quiet

bool debug()
{
#ifdef DEBUG_JobOptsSvc
  return true;
#else
  return false;
#endif
}

JobOptsSvc* JobOptsSvc::p_jobOptsSvc = NULL;

JobOptsSvc* JobOptsSvc::instance()
{
  if(debug()) cout << "JobOptsSvc::instance()" << endl;

  if(NULL == p_jobOptsSvc)
    {
      p_jobOptsSvc = new JobOptsSvc;
      
      if(debug()) cout << " creating new JobOptsSvc object" << endl;
      p_jobOptsSvc->init();
    }

  return p_jobOptsSvc;
}

void JobOptsSvc::close()
{
  if(debug()) cout << "   JobOptsSvc::close()" << endl;
  
  if(NULL != p_jobOptsSvc)
    {
      delete p_jobOptsSvc;
    }
  else
    {
      std::cout << "Error: no instance of job options service found! " << std::endl;
    }
}

bool JobOptsSvc::init()
{
  if(debug()) cout << "JobOptsSvc::init()" << endl;

  //hardcoded default standard job options file
  return init("$KTRACKER_ROOT/opts/default.opts");
}

bool JobOptsSvc::init(const char* configfile)
{
  if(debug()) cout << "JobOptsSvc::init( " << configfile << " )" << endl;

  // expand any environment variables in the file name
  m_configFile = ExpandEnvironmentals( configfile );
  if(debug()) cout << "My config file is :" << m_configFile << endl;

  std::ifstream fin( m_configFile.c_str() );
  if(!fin) 
    {
      cout << "JobOptsSvc::init - ERROR - Your config file '" << m_configFile << "' cannot be opened!" << endl;
      throw 1;
    }


  //define the types of options to look for
  map<string,string*> stringOpts;
  stringOpts["InputFile"] = &m_inputFile;
  stringOpts["OutputFile"] = &m_outputFile;

  stringOpts["AlignmentFile_Hodo"] = &m_alignmentFileHodo;
  stringOpts["AlignmentFile_Mille"] = &m_alignmentFileMille;
  stringOpts["AlignmentFile_Prop"] = &m_alignmentFileProp;
  stringOpts["AlignmentFile_Chamber"] = &m_alignmentFileChamber;

  stringOpts["RoadFile_plus_top"] = &m_roadFile_pt;
  stringOpts["RoadFile_plus_bottom"] = &m_roadFile_pb;
  stringOpts["RoadFile_minus_top"] = &m_roadFile_mt;
  stringOpts["RoadFile_minus_bottom"] = &m_roadFile_mb;

  stringOpts["CalibrationsFile"] = &m_calibrationsFile;

  stringOpts["fMagFile"] = &m_fMagFile;
  stringOpts["kMagFile"] = &m_kMagFile;

  stringOpts["MySQL_Server"] = &m_mySQLServer;
  stringOpts["Geometry_Version"] = &m_geomVersion;

  map<string,int*> intOpts;
  intOpts["MySQL_Port"] = &m_mySQLPort;
  intOpts["N_Events"] = &m_nEvents;
  intOpts["FirstEvent"] = &m_firstEvent;

  map<string,bool*> boolOpts;
  boolOpts["MCMode_enable"] = &m_mcMode;
  boolOpts["AlignmentMode_enable"] = &m_alignmentMode;
  boolOpts["TriggerMask_enable"] = &m_enableTriggerMask;
  boolOpts["kMag_enable"] = &m_enableKMag;
  boolOpts["Evaluation_enable"] = &m_enableEvaluation;
  boolOpts["OnlineAlignment_enable"] = &m_enableOnlineAlignment;

  //read the file and find matching options
  string line;
  while(getline(fin, line))
    {
      if(debug()) cout << "line = " << line << endl;
      boost::algorithm::trim(line);

      //skip empty and comment lines
      if(line.empty() || line[0]=='#') continue;

      stringstream ss(line);
      string key, val;
      ss >> key >> val;
      if(val.empty())
	{
	  if(debug()) cout << "JobOptsSvc::init - WARNING - caught std::logic_error on line.  Possible value missing in line: " << line << endl;
	  continue;
	}


      if(debug()) cout << " key = " << key << ", val = " << val << endl;

      //is this a string option?
	{
	  map<string,string*>::iterator it = stringOpts.find(key);
	  if(stringOpts.end() != it)
	    {
	      string expandedVal = ExpandEnvironmentals( val );
	      it->second->assign(expandedVal);
	      if(debug()) cout << " ... assign string key" << endl;
	      continue;
	    }
	}

      //is this an int option?
	{
	  map<string,int*>::iterator it = intOpts.find(key);
	  if(intOpts.end() != it )
	    {
	      *(it->second) = atoi(val.c_str());
	      if(debug()) cout << " ... assign int key" << endl;
	      continue;
	    }
	}

      //is this a bool option?
	{
	  map<string,bool*>::iterator it = boolOpts.find(key);
	  if(boolOpts.end() != it)
	    {
	      *(it->second) = atoi(val.c_str()); //todo:support true/false
	      if(debug()) cout << " ... assign bool key" << endl;
	      continue;
	    }
	}

      if(debug()) cout << " ... key not found.  handle error?" << endl;
    }

  m_mySQLurl = "mysql://" + m_mySQLServer + ":" + boost::lexical_cast<string>(m_mySQLPort);
  return true;
}



std::string JobOptsSvc::ExpandEnvironmentals( const std::string& input ) const
{
  // expand any environment variables in the file name
  wordexp_t exp_result;
  if(wordexp(input.c_str(), &exp_result, 0) != 0)
    {
      //this is a fatal error so throw
      cout << "JobOptsSvc::init - ERROR - Your string '" << input << "' cannot be understood!" << endl;
      throw 1;
    }

  const string output( exp_result.we_wordv[0]);
  if(debug()) cout << "JobOptsSvc::ExpandEnvironmentals - expanded " << input << " -> " << output << endl;

  return output;
}

