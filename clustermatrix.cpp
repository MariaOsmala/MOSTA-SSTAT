#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "clustermatrix.h"

using namespace std;

CClusterMatrix::CClusterMatrix(string &asflist,double &gc,string &astmethod,double &at, bool absge) : CSimilarityMatrix(asflist,gc,astmethod,at),bsge(absge)
{
}

//executes clustering: get similarity matrix and merge successivley.
// parameter:
//   p: percent quantile for smax_limit
void CClusterMatrix::do_cluster(double p)
{
  //call inherited method to fill matrix
  fill_matrix((int)bsge-2,true);

  //get highest similarity and quantiles
  double smax =0;
  int imax = -1;
  int jmax = -1;
  vector<double> vdistr;
  vdistr.reserve(n*(n-1)/2);
  for (int i=0;i<n;i++)
    {
      for (int j=i+1;j<n;j++)
	{
	  vdistr.push_back(m[i*nm+j]);
	  //save if maximum
	  if (m.at(i*nm+j)>smax)
	    {
	      smax = m[i*nm+j];
	      imax = i;
	      jmax = j;
	    }
	}
    }
  //sort distribution
  sort(vdistr.begin(),vdistr.end());
  //get p% quantile
  double smax_limit = max(vdistr[(int)(p*vdistr.size())]-.00000001,0.0);
  cerr << "smax_limit=" << smax_limit << endl;

  //create file for matrix output
  stringstream sfmatrices;
  sfmatrices << sflist << ".matrices";
  ofstream fmatrices(sfmatrices.str().c_str());  

  //prepare member vector
  vector< vector<int> > vmembers;
  vmembers.reserve(n);
  int koffset = n;

  //print header
  cout << "matrixA\tmatrixB\tQA\tQB\ticA\ticB\tSmax\timax\tbimaxp\tQ\tic\n";
  //now merge maximum, get highest similarity and iterate k-1 times
  while (smax>smax_limit)
    {
      //merge maximum pair
      CPfm opfm_new(vopfm[imax],vopfm[jmax],mpos[imax*nm+jmax],(bool)mpos[jmax*nm+imax]);

      cerr << "trying: " << vopfm[imax].sid << " with " << vopfm[jmax].sid << endl;

      //get members of new cluster
      vector<int> vm;
      if (imax<koffset)
	vm.push_back(imax);
      else
	{
	  for (int km=0;km<vmembers[imax-koffset].size();km++)
	    vm.push_back(vmembers[imax-koffset].at(km));
	}
      if (jmax<koffset)
	vm.push_back(jmax);
      else
	{
	  for (int km=0;km<vmembers[jmax-koffset].size();km++)
	    vm.push_back(vmembers[jmax-koffset].at(km));
	}
      
      //adjust matrix
      opfm_new.adjust_t(stmethod,t);
      //check similarity of members to new representative
      cerr << "\tmember\tSmax\n";
      int nmember = vm.size();
      vector<bool> vbs(nmember,false);
      #pragma omp parallel default(shared)
      {
        #pragma omp for
	for (int km=0;km<nmember;km++)
	  {
	    int i = vm[km];
	    CSimStat osstat(vopfm[i],opfm_new);
	    cerr << "\t" << vopfm[i].sid << "\t" << osstat.smax << endl;
	    vbs[km] = osstat.smax>smax_limit;
	  }
      }
      //check whether criterium is fulfilled
      bool bnolimit=true;
      for (int km=0;km<nmember;km++)
	bnolimit = bnolimit && vbs[km];

      //check increase of information content
      if (bnolimit)
	{
	  cerr << "\t--> MERGED!\n";
	  //output
	  cout << vopfm[imax].sid << "\t" << vopfm[jmax].sid << "\t";
	  cout << vopfm[imax].get_beta() << "\t" << vopfm[jmax].get_beta() << "\t";
	  cout << vopfm[imax].get_ic() << "\t" << vopfm[jmax].get_ic() << "\t";
	  cout << smax << "\t" << mpos[imax*nm+jmax] << "\t" << mpos[jmax*nm+imax] << "\t";
	  
	  //add new matrix
	  int knew=n;
	  n++;
	  vopfm.push_back(opfm_new);
	  bused[knew] = true;

	  //add members
	  vmembers.push_back(vm);


	  //remove old matrices
	  bused[imax] = false;
	  bused[jmax] = false;

	  //output matrix into .matrices file
	  fmatrices << vopfm[knew];

	  //output beta, ic, ...
	  cout << vopfm[knew].get_beta() << "\t";
	  cout << vopfm[knew].get_ic() << endl;

	  //compute similarity to new pfm
          #pragma omp parallel default(shared)
	  {
            #pragma omp for
	    for (int i=0;i<knew;i++)
	      {
		if (bused[i])
		  {
		    CSimStat osstat(vopfm[i],vopfm[knew]);
		    m[i*nm+knew] = osstat.smax;
		    mpos[i*nm+knew] = osstat.imax;
		    mpos[knew*nm+i] = (int)osstat.bimaxp;
		  }
	      }
	  }
	} else
	{
	  //set smax to negative value since ic decreases
	  m[imax*nm+jmax] = smax_limit-1;
	}

      //NOW, find maximum
      smax = smax_limit;
      for (int i=0;i<n;i++)
	{
	  if (bused[i])
	    {
	      for (int j=i+1;j<n;j++)
		{
		  if (bused[j])
		    {
		      //save if maximum
		      if (m[i*nm+j]>smax)
			{
			  smax = m[i*nm+j];
			  imax = i;
			  jmax = j;
			}
		    }
		}
	    }
	}
    }
  fmatrices.close();
}

//prints transfac files for non-merged matrices
// parameter:
//   osout: output stream
void CClusterMatrix::print_pfms(ostream& osout)
{
  //iterate over matrices
  for (int i=0;i<bused.size();i++)
    {
      //print them if not merged
      if (bused[i])
	osout << vopfm[i];
    }
}
