#ifndef SIMSTAT_H
#define SIMSTAT_H

#include <vector>
#include <cmath>
#include <ctime>

#include "pfm.h"
#include "countpars.h"
#include "mytimer.h"

class CSimStat
{
 public:
  CSimStat(CPfm &opfmA,CPfm &opfmB);
  double trun; //in s
  double smax;
  double ssum;
  int imax;
  bool bimaxp; //is on palindrome
 private:
  CTimer otimer;
  CCountPars opars;
  long double precision;
};

#endif
