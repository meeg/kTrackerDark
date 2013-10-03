#ifndef TubeField_H
#define TubeField_H 1

#include "globals.hh"
#include "G4ElectricField.hh"
#include <cmath>
#include <iostream>
#include <cstdlib>

class TubeField : public G4ElectricField
{
  public:
    TubeField();
    ~TubeField();

    void GetFieldValue(const double Point[3], double *Efield) const;
    void SetOuterRadius(double);
    void SetTubeRadius(double);
    void SetZPlanes(double, double, double, double);

    // OuterRadius is outer radius of the gas.  tubeRadius is outer radius of the metal tube.
    double voltage, outerRadius, innerRadius, tubeRadius;
    double derivedConstant;
    double z1v, z2v, z1h, z2h;
};

#endif
