
#include "scluster.h"
#include "pfm_helper.h"
#include "clustermatrix.h"

using namespace std;

int main(int argc, char* argv[])
{
  //check command line
  if (argc<5)
    {
      cout << "Insufficient parameters\nCall scluster <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<use-sge>] [<p=.95>]\n";
      return(1);
    }

  //read parameters
  //get gc content
  double gc;
  istringstream istr1(argv[1]);
  istr1 >> gc;
  if ((gc<=0) || (gc>=1))
    {
      cerr << "Error in input: given gc content " << gc << " is not between 0 and 1.\n";
      exit(1);
    }

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //and corresponding parameter
  double t;
  istringstream istr4(argv[4]);
  istr4 >> t;

  //use SGE engine?
  bool bsge = false;
  if (argc>5)
    {
      istringstream istr5(argv[5]);
      istr5 >> bsge;
    }

  //define p?
  double p = .95;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> p;
    }

   //get matrices
  string sflist(argv[2]);

  //prepare matrix to hold similarities
  CClusterMatrix ocm(sflist,gc,stmethod,t,bsge);
  ocm.do_cluster(p);

  //create file for matrix output
  stringstream sfmatout;
  sfmatout << sflist << ".cluster";
  ofstream fmatout(sfmatout.str().c_str());  
  ocm.print_pfms(fmatout);
  fmatout.close();
  
  return(0);
}
