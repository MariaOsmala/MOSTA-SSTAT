#include <cstdlib>
#include "costat.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall costat <gc> <transfac-file> <threshold-method> <threshold-parameter> [<window-size>] [bregularize]\n");
  //check command line
  if (argc<4)
    {
      cout << cswrongargs;
      return(1);
    }

  bool bregularize=1;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> bregularize;
    }

  //window size
  long nseq = -1;
  if (argc>5)
    {
      istringstream istr5(argv[5]);
      istr5 >> nseq;
    }

  //get gc content
  double gc;
  istringstream istrgc(argv[1]);
  istrgc >> gc;
  if ((gc<=0) || (gc>=1))
    {
      cerr << "Error in input: given gc content " << gc << " is not between 0 and 1.\n";
      exit(1);
    }

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //get threshold parameter
  double tp;
  istringstream istr4(argv[4]);
  istr4 >> tp;

  //matrix
  CPfmLoader vopfm(string(argv[2]),gc,bregularize);
  
  //adjust matrices
  for (int i=0;i<vopfm.size();i++)
    {
      //create PFM object
      vopfm[i].adjust_t(stmethod,tp);
    }

  //pair-wise similarity probabilities
  if (nseq>0)
    cout << "matrixA\tmatrixB\tp\n";
  else
    cout << "matrixA\tmatrixB\trA\trB\trAB\talphaA\talphaB\n";
  for (int i=0;i<vopfm.size();i++)
    {
      for (int j=i+1;j<vopfm.size();j++)
	{
	  //create object for the cooccurences
	  CCoocStat ocstat(vopfm[i],vopfm[j],gc);
	  if (nseq>0)
	    {
	      long double p = ocstat.calc_winprob(nseq);
	      //output
	      cout << vopfm[i].sid << "\t" << vopfm[j].sid << "\t" << p << endl;
	    } 
	  else
	    cout << vopfm[i].sid << "\t" << vopfm[j].sid << "\t" << ocstat.rA << "\t" << ocstat.rB << "\t" << ocstat.rAB << "\t" << ocstat.alphaA << "\t" << ocstat.alphaB << endl;
	}
    }

  return(0);
}
