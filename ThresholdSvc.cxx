#include "ThresholdSvc.h"

using namespace Threshold;

//Singleton gettor
ThresholdSvc& ThresholdSvc::Get()
{
    static ThresholdSvc instance;
    return instance;
}

//Constructor
ThresholdSvc::ThresholdSvc() :
    m_level( kInfo ),
    m_live(false)
{ }


//settors and gettors
void ThresholdSvc::SetLevel( Level level )
{
    m_level = level;
}
void ThresholdSvc::SetLive( bool live )
{
    m_live = live;
}


//check if thresholds are passed
bool ThresholdSvc::DoError()   const
{
    return kError <= m_level;
}
bool ThresholdSvc::DoWarning() const
{
    return kWarning <= m_level;
}
bool ThresholdSvc::DoInfo()    const
{
    return kInfo <= m_level;
}
bool ThresholdSvc::DoDebug()   const
{
    return kDebug <= m_level;
}
bool ThresholdSvc::DoVerbose() const
{
    return kVerbose <= m_level;
}

bool ThresholdSvc::DoLive() const
{
    return m_live;
}


//shortcuts for ostream
namespace Threshold
{
bool live()
{
    return ThresholdSvc::Get().DoLive();
}
bool error()
{
    return ThresholdSvc::Get().DoError();
}
bool warning()
{
    return ThresholdSvc::Get().DoWarning();
}
bool info()
{
    return ThresholdSvc::Get().DoInfo();
}
bool debug()
{
    return ThresholdSvc::Get().DoDebug();
}
bool verbose()
{
    return ThresholdSvc::Get().DoVerbose();
}
}
