//Class returns gamma's for the statistics

#ifndef COUNTPARS_H
#define COUNTPARS_H

#include <vector>
#include <iostream>

#include "pfm.h"
#include "convolution.h"

class CCountPars
{
 public:
  CCountPars(CPfm &opfm,double seqgc=-1);
  CCountPars(CPfm &opfmA,CPfm &opfmB,double seqgc=-1);
  void get_pars();
  double alpha;
  std::vector<double> gamma;
  std::vector<double> gammap;
  void do_enlarge_left();
 private:
  int tA;
  int tB;
  std::vector<double> vbg;
  std::vector<int> mA;
  std::vector<int> mB;
  void print_Qsub(std::vector<double> &Qsub,int nrangeA,int nrangeB);
  void get_gamma(std::vector<int> &mA,int tA,std::vector<int> &mB,int tB,std::vector<double> &vbg,std::vector<double> &mygamma);
  void enlarge_right(std::vector<int> &m,int n);
  void enlarge_left(std::vector<int> &m,int n);
};

#endif
