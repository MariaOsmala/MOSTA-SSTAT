#include <cstdlib>
#include "bsanno.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall bsanno <sequence file> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<bregularize>] [<gc-content for global option>]\n");
  //check command line
  if (argc<5)
    {
      cout << cswrongargs;
      return(1);
    }

  //read parameters

  bool bregularize=1;
  if (argc>5)
    {
      istringstream istr(argv[5]);
      istr >> bregularize;
    }

  double gc=-1;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> gc;

      if ((gc<=0) || (gc>=1))
	{
	  cerr << "Error in input: given gc content " << gc << " is not between 0 and 1.\n";
	  exit(1);
	}

    }

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //get threshold parameter
  double tp;
  istringstream istr4(argv[4]);
  istr4 >> tp;

  //PFM, set gc content arbitrarily since it'll be reinit!
  CPfmLoader vopfm(string(argv[2]),.5,bregularize);

  cout << CBSAnnotator::get_header();

  //iterate over sequences
  CSequences oseqs(argv[1]);
  while (!oseqs.eof())
    {
      CSequence oseq = oseqs.next();

      //iterate over matrices
      for (int i=0;i<vopfm.size();i++)
	{
	  double igc=gc;
	  if (igc<0)
	    igc = oseq.get_gc();

	  //create PFM object
	  vopfm[i].reinit(igc,bregularize);
	  vopfm[i].adjust_t(stmethod,tp);

	  //create annotator
	  CBSAnnotator obs(vopfm[i]);
	  //annotate
	  obs.annotate(oseq);
	  
	  cout << obs;
	}
    }
  return(0);
}
