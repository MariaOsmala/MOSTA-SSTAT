#include "simstat.h"

using namespace std;

//imax is always B shifted against A, and bimaxp is true
//if B has to be the reverse complementary pfm.
CSimStat::CSimStat(CPfm &opfmA,CPfm &opfmB) : otimer(),opars(opfmA,opfmB),precision(.000000001)
{
  //get timer
  otimer.start();
  
  //enlarge matrix to the left
  opars.do_enlarge_left();
  //get parameter
  opars.get_pars();

  //get length of matrices
  int nA=opfmA.nlen;
  int nB=opfmB.nlen;

  ssum=0;
  smax=0;

  //get alpha
  double alpha;
  double beta;
  if (nA>=nB)
    {
      alpha = opfmB.get_alpha();
      beta = opfmA.get_alpha();
    } else 
    {
      alpha = opfmA.get_alpha();
      beta = opfmB.get_alpha();
    }

  //constant definition and variable declarations
  long double epsilon = .0000000001;
  //read parameters from CCountPars and get similarity measures
  for (int i=0;i<nA+nB-1;i++)
    {
      //      cerr << "i=" << i << "\t" << opars.gamma[i]*opars.alpha << "\t" << opars.gammap[i]*opars.alpha << endl;
      ssum += opars.gamma[i];
      if (opars.gamma[i]>smax)
	{
	  smax = opars.gamma[i];
	  imax = i;
	  //no palindrome
	  bimaxp = 0;
	}
      //for the complementary strand
      ssum += opars.gammap[i];
      if (opars.gammap[i]>smax)
	{
	  smax = opars.gammap[i];
	  imax = i;
	  //palindromic hit
	  bimaxp=1;
	}
    }
  //is palindrome?
  //get correct imax position
  int inull =0;
  int nv=0;
  int nu=0;
  if (nA>=nB)
    {
      imax = imax-nB+1;
      nv = nB;
      nu = nA;
      inull = nB-1;
    } else 
    {
      imax = -(imax-nA+1);
      nv = nA;
      nu = nB;
      inull = nA-1;
      if (bimaxp)
	imax = -nB+nA-imax;
    }

  //get covariance, v is smaller (corresp. to alpha!)
  long double Guv=0;
  long double Gvu=0;
  long double cov=0;
  for (int i=0;i<nA+nB-1;i++)
    {
      //      cerr << i << "\t" << opars.gamma[i] << "\t" << opars.gammap[i] << endl;
      if (i<=inull)
	{
	  Gvu += opars.gamma[i];
	  Gvu += opars.gammap[i];
	}
      if (i>=inull)
	{
	  Guv += opars.gamma[i];
	  Guv += opars.gammap[i];
	}
    }

  Guv *= opars.alpha;
  Gvu *= opars.alpha;
  
  cov += Guv + Gvu;
  cov -= 2 * (nu+nv) * alpha * beta;
  cov -= alpha * ((opars.gamma[inull]+opars.gammap[inull])-2*beta);
  cov *= 2;
  //  cerr << "Guv=\t" << Guv << endl;
  //cerr << "Gvu=\t" << Gvu << endl;
  //cerr << "pi_u=\t" << beta << endl;
  //cerr << "pi_v=\t" << alpha << endl;
  //cerr << "inull=\t" << inull << endl;

  //adjust similarity with alpha
  int np = (nA+nB-1)*2;
  ssum = log(ssum/np/alpha);
  //  ssum = log(ssum/alpha);
  smax = log(smax/alpha);

  ssum = cov;

  //stop time measuring
  trun = otimer.get_time();
}
