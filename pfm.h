#ifndef PFM_H
#define PFM_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>

#include "pfm_helper.h"
#include "convolution.h"
#include "stringutils.h"

//Exceptions
class EThresholdNotSet{};

class CPfm {
 public:
  double nseqs;
  double ic_over;
  int n_over;
  std::string sid;
  int nlen;
  std::vector<int> mpssm;
  std::vector<double> vbg;
  bool bregularize;
  CPfm(std::string sfn,double gc,bool bregularize);
  CPfm(std::ifstream &f,double gc,bool bregularize);
  CPfm(CPfm& opfmA,CPfm& opfmB,int ishift,bool bcompl);
  int get_t() throw (EThresholdNotSet);
  double get_alpha(int m=-1) throw (EThresholdNotSet);
  double get_beta() throw (EThresholdNotSet);
  void adjust_t(std::string smethod,double tp=-23880) throw(EThresholdNotSet);
  //returns information content
  double get_ic(int istart=0,int nr=-1);
  //returns gc content
  double get_gc();
  void reinit(double igc, bool bregularize);
 private:
  double ic;
  double alpha;
  double beta;
  bool bt; //threshold set or not?
  int t;
  double gc;
  std::vector<double> mpcm;
  std::vector<long double> mpwm;
  void read_pcm(std::ifstream &f);
  void pcm2pssm(bool bregularize);
  friend std::ostream& operator<< (std::ostream& strm, const CPfm& obj);
};


#endif

