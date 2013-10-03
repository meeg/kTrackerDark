#include "TubeField.hh"

using namespace std;

TubeField::TubeField()
{
  voltage = 1850*volt;
  innerRadius = 20*micrometer;
}

TubeField::~TubeField()
{
}

void TubeField::SetOuterRadius(double outRad)
{
  outerRadius = outRad;
  derivedConstant = 2*pi*epsilon0*voltage/log(outerRadius/innerRadius);
}

void TubeField::SetTubeRadius(double tubRad)
{
  tubeRadius = tubRad;
}

void TubeField::SetZPlanes(double zz1v, double zz2v, double zz1h, double zz2h)
{
  z1v = zz1v;
  z2v = zz2v;
  z1h = zz1h;
  z2h = zz2h;
}

void TubeField::GetFieldValue(const double Point[3], double *Bfield) const
{
  Bfield[0] = 0.;
  Bfield[1] = 0.;
  Bfield[2] = 0.;
  Bfield[3] = 0.;
  Bfield[4] = 0.;
  Bfield[5] = 0.;

  double px = Point[0];
  double py = Point[1];
  double pz = Point[2];
  bool v = false;

  if (abs(z1v-pz) < tubeRadius*2.0)
  {
    pz = pz - z1v;
    v = true;
  }
  else if (abs(z2v-pz) < tubeRadius*2.0)
  {
    pz = pz - z2v;
    v = true;
  }
  else if (abs(z1h-pz) < tubeRadius*2.0)
  {
    pz = pz - z1h;
    v = false;
  }
  else if (abs(z2h-pz) < tubeRadius*2.0)
  {
    pz = pz - z2h;
    v = false;
  }
  else
    return;

  if (v == true)
  {
    double tempx1 = fmod(abs(px - tubeRadius/2.0), tubeRadius*2.0);
    if (tempx1 > tubeRadius)
      tempx1 = fmod(abs(px - tubeRadius/2.0), tubeRadius);
    else
      tempx1 = tubeRadius - fmod(abs(px - tubeRadius/2.0), tubeRadius);
    double tempx2 = fmod(abs(px + tubeRadius/2.0), tubeRadius*2.0);
    if (tempx2 > tubeRadius)
      tempx2 = fmod(abs(px + tubeRadius/2.0), tubeRadius);
    else
      tempx2 = tubeRadius - fmod(abs(px + tubeRadius/2.0), tubeRadius);
    double tempz1 = pz + tubeRadius*sqrt(3)/2;
    double tempz2 = pz - tubeRadius*sqrt(3)/2;

    double tempr1 = sqrt(tempx1*tempx1 + tempz1*tempz1);
    double tempr2 = sqrt(tempx2*tempx2 + tempz2*tempz2);

    if (tempr1 > outerRadius && tempr2 > outerRadius)
      return;
    if (tempr1 < outerRadius && tempr2 < outerRadius)
      return;

    if (tempr1 < outerRadius)
    {
      px = tempx1;
      pz = tempz1;
      if (px > 36.5*tubeRadius || px < -35.5*tubeRadius)
        return;
    }
    if (tempr2 < outerRadius)
    {
      px = tempx2;
      pz = tempz2;
      if (px > 35.5*tubeRadius || px < -36.5*tubeRadius)
        return;
    }
    if (py > 36.5*tubeRadius || py < -36.5*tubeRadius)
      return;

    double r = sqrt(px*px+pz*pz);
    if (r < innerRadius)
      r = innerRadius;

    double strength = derivedConstant/(2*pi*epsilon0*r);

    Bfield[3] = strength*px/r;
    Bfield[4] = 0.;
    Bfield[5] = strength*pz/r;
  }
  else
  {
    double tempy1 = fmod(abs(py - tubeRadius/2.0), tubeRadius*2.0);
    if (tempy1 > tubeRadius)
      tempy1 = fmod(abs(py - tubeRadius/2.0), tubeRadius);
    else
      tempy1 = tubeRadius - fmod(abs(py - tubeRadius/2.0), tubeRadius);
    double tempy2 = fmod(abs(py + tubeRadius/2.0), tubeRadius*2.0);
    if (tempy2 > tubeRadius)
      tempy2 = fmod(abs(py + tubeRadius/2.0), tubeRadius);
    else
      tempy2 = tubeRadius - fmod(abs(py + tubeRadius/2.0), tubeRadius);

    double tempz1 = pz + tubeRadius*sqrt(3)/2;
    double tempz2 = pz - tubeRadius*sqrt(3)/2;

    double tempr1 = sqrt(tempy1*tempy1 + tempz1*tempz1);
    double tempr2 = sqrt(tempy2*tempy2 + tempz2*tempz2);

    if (tempr1 > outerRadius && tempr2 > outerRadius)
      return;
    if (tempr1 < outerRadius && tempr2 < outerRadius)
      return;

    if (tempr1 < outerRadius)
    {
      py = tempy1;
      pz = tempz1;
      if (py > 36.5*tubeRadius || py < -35.5*tubeRadius)
        return;
    }
    if (tempr2 < outerRadius)
    {
      py = tempy2;
      pz = tempz2;
      if (py > 35.5*tubeRadius || py < -36.5*tubeRadius)
        return;
    }
    if (px > 36.5*tubeRadius || px < -36.5*tubeRadius)
      return;

    double r = sqrt(py*py+pz*pz);
    if (r < innerRadius)
      r = innerRadius;

    double strength = derivedConstant/(2*pi*epsilon0*r);

    Bfield[3] = 0.;
    Bfield[4] = strength*py/r;
    Bfield[5] = strength*pz/r;
  }
}
