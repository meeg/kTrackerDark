#include "../JobOptsSvc.h"

#include <iostream>

using namespace std;

int main()
{
  JobOptsSvc* jobOptsSvc = JobOptsSvc::instance();

  //try functions
  cout << jobOptsSvc->GetRoadsFilePlusTop() << endl;
  cout << jobOptsSvc->GetRoadsFilePlusBottom() << endl;
  cout << jobOptsSvc->GetRoadsFileMinusTop() << endl;
  cout << jobOptsSvc->GetRoadsFileMinusBottom() << endl;

  jobOptsSvc->close();

  return 0;
}
