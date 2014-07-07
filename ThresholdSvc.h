#ifndef THRESHOLDSVC_H
#define THRESHOLDSVC_H

namespace Threshold
{
  enum Level
  {
    kNone = 0,
    kError,
    kWarning,
    kInfo,
    kDebug,
    kVerbose
  };

  class ThresholdSvc
  {
    public:
      ///Singleton gettor
      static ThresholdSvc& Get();

      ///Set the level
      void SetLevel( Level level );
      ///Set live mode
      void SetLive( bool live );

      ///Is threshold at or above error
      bool DoError() const;
      ///Is threshold at or above warning
      bool DoWarning() const;
      ///Is threshold at or above info
      bool DoInfo() const;
      ///Is threshold at or above debug
      bool DoDebug() const;
      ///Is threshold at or above verbose
      bool DoVerbose() const;


      ///Do we want to do stuff for live mode?
      bool DoLive() const;

    private:
      ///Default constructor (private)
      ThresholdSvc();

      /// No copy
      ThresholdSvc(ThresholdSvc const&);              // Don't Implement
      /// No assignment
      void operator=(ThresholdSvc const&); // Don't implement

      Level m_level;
      bool m_live;
  };

  //short cuts
  bool live();
  bool error();
  bool warning();
  bool info();
  bool debug();
  bool verbose();
}

#endif
