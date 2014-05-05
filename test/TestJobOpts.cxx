#include "../kTrackerServices/JobOptsSvc.h"

int main()
{
  JobOptsSvc* jobOptsSvc = JobOptsSvc::instance();
  jobOptsSvc->init();

  jobOptsSvc->close();

  return 0;
}
