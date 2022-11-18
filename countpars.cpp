#include "countpars.h"
#include "pfm_helper.h"

using namespace std;

CCountPars::CCountPars(CPfm &opfm,double seqgc)
{
  //initialize variables
  tA = opfm.get_t();
  tB = opfm.get_t();
  if (seqgc>=0)
    {
      vbg.reserve(4);
      vbg.push_back((1-seqgc)/2);
      vbg.push_back(seqgc/2);
      vbg.push_back(seqgc/2);
      vbg.push_back((1-seqgc)/2);
    } else 
    {
      vbg = opfm.vbg;
    }
  mA = opfm.mpssm;
  //for extension to two matrices
  mB = opfm.mpssm;
  alpha = opfm.get_alpha();
}

//ATTENTION: we use bg distribution from A!
CCountPars::CCountPars(CPfm &opfmA,CPfm &opfmB,double seqgc)
{
  if (seqgc>=0)
    {
      vbg.reserve(4);
      vbg.push_back((1-seqgc)/2);
      vbg.push_back(seqgc/2);
      vbg.push_back(seqgc/2);
      vbg.push_back((1-seqgc)/2);
    }
  //initialize variables
  if (opfmA.nlen>=opfmB.nlen)
    { 
      tA = opfmA.get_t();
      tB = opfmB.get_t();
      if (seqgc<0)
	vbg = opfmA.vbg;
      mA = opfmA.mpssm;
      mB = opfmB.mpssm;
      alpha = opfmA.get_alpha();
    } else // swap matrices such that A is larger than B
    {
      tA = opfmB.get_t();
      tB = opfmA.get_t();
      if (seqgc<0)
	vbg = opfmA.vbg;
      mA = opfmB.mpssm;
      mB = opfmA.mpssm;
      alpha = opfmB.get_alpha();
    }
}

// enlarges bigger matrix to the left
// we need that for getting negative shift for similarity measure
void CCountPars::do_enlarge_left()
{
  enlarge_left(mA,mA.size()+mB.size()-4);
}

// enlarges given matrix to the left
// parameter:
//   n: final size!
void CCountPars::enlarge_left(vector<int> &m,int n)
{
  if (m.size()<n)
    {
      int nloc = m.size();
      vector<int> mt=m;
      m.clear();
      m.reserve(n);
      for (int i=0;i<(n-nloc);i++)
	m.push_back(0);
      for (int i=0;i<nloc;i++)
	m.push_back(mt[i]);
    }
}
// enlarges given matrix to the right
// parameter:
//   n: final size!
void CCountPars::enlarge_right(vector<int> &m,int n)
{
  if (m.size()<n)
    {
      int nloc = m.size();
      m.reserve(n);
      for (int i=0;i<(n-nloc);i++)
	{
	  m.push_back(0);
	}
    }
}

void CCountPars::get_pars()
{
  //get reverse complementary matrix
  vector<int> mBc;
  pfm_reversecompl(mB,mBc);

  //make matrices equally large
  int nlen = mA.size();
  enlarge_right(mB,nlen);
  enlarge_right(mBc,nlen);

  get_gamma(mA,tA,mB,tB,vbg,gamma);
  get_gamma(mA,tA,mBc,tB,vbg,gammap);
}

void CCountPars::get_gamma(vector<int> &mA,int tA,vector<int> &mB,int tB,vector<double> &vbg,vector<double> &mygamma)
{
  int nlen = mA.size()/4;

  //get min and max for each row of matrix
  vector<int> vminA;
  vector<int> vmaxA;
  pfm_minmax(mA,vminA,vmaxA);
  vector<int> vminB;
  vector<int> vmaxB;
  pfm_minmax(mB,vminB,vmaxB);
  
  //get differences between min and max
  vector<int> vdiff;
  vdiff.reserve(nlen);
  for (int i=0;i<nlen;i++)
    vdiff.push_back(vminA[i]-vmaxA[i]); //returns negative differences for sort!
  //get new order of traversal
  vector<int> vi; //contains value-sorted indices
  pfm_sortindices(vdiff,vi);

  //get vmin and vmax for new coordinates
  vector<int> xvminA;
  vector<int> xvmaxA;
  xvminA.reserve(nlen);
  xvmaxA.reserve(nlen);
  vector<int> xvminB;
  vector<int> xvmaxB;
  xvminB.reserve(nlen);
  xvmaxB.reserve(nlen);
  for (int i=0;i<nlen;i++)
    {
      xvminA.push_back(vminA[vi[i]]);
      xvmaxA.push_back(vmaxA[vi[i]]);
      xvminB.push_back(vminB[vi[i]]);
      xvmaxB.push_back(vmaxB[vi[i]]);
    }

  //initialize (before first step)
  //get most positive and most negative score
  int xvminiA = 0;
  int xvmaxiA = 0;
  int tvminiB = 0;
  int tvmaxiB = 0;
  for (int i=0;i<nlen;i++)
    {
      xvminiA += xvminA[i];
      xvmaxiA += xvmaxA[i];
      tvminiB += xvminB[i];
      tvmaxiB += xvmaxB[i];
    }  
  vector<int> xvminiB(nlen,tvminiB);
  vector<int> xvmaxiB(nlen,tvmaxiB);
  
  //minimal/maximal score from last step
  int noldloA = 0;
  int noldupA = 0;
  //initialize vectors holding min/max for B
  vector<int> voldloB(nlen,0);
  vector<int> voldupB(nlen,0);

  //declare and initialize Q
  vector<vector<vector<double> > > Q(nlen,vector<vector<double> >(2,vector<double>(4,0)));
  for (int ishift=0;ishift<nlen;ishift++)
    {
      Q[ishift][0][0] = 1;
    }
  

  //remember old/new index per shift for invariant traversals
  vector<int> vold(nlen,0);
  vector<int> vnew(nlen,1);
  //start iteration
  for (int i=0;i<nlen;i++)
    {
      int noldrangeA = noldupA-noldloA+2;
      //get lower/upper bound for scores
      xvmaxiA -= xvmaxA[i];
      xvminiA -= xvminA[i];
      int nupA = min(noldupA+xvmaxA[i],tA-xvminiA);
      int nloA = max(noldloA+xvminA[i],tA-xvmaxiA);
      int nrangeA = nupA-nloA+2; //plus 1 for certain s>t

      //iterate over all shifts
      for (int ishift=0;ishift<nlen;ishift++)
	{
	  //position on mB
	  int ipos = vi[i]-ishift;
	  int iold = vold[ishift];
	  int inew = vnew[ishift];

	  //get lower/upper bound for scores
	  int noldrangeB = voldupB[ishift] - voldloB[ishift] + 2;
	  int nupB;
	  int nloB;
	  int nrangeB;
	  if (ipos>=0)
	    {
	      xvmaxiB[ishift] -= vmaxB[ipos];
	      xvminiB[ishift] -= vminB[ipos];
	      nupB = min(voldupB[ishift]+vmaxB[ipos],tB-xvminiB[ishift]);
	      nloB = max(voldloB[ishift]+vminB[ipos],tB-xvmaxiB[ishift]);
	      nrangeB = nupB-nloB+2; //plus 1 for surely s>t
	    } else 
	    {
	      nupB = voldupB[ishift];
	      nloB = voldloB[ishift];
	      nrangeB = noldrangeB;
	    }

	  //initialize Q
	  Q[ishift][inew].clear();
	  Q[ishift][inew].resize(nrangeA*nrangeB,0);

	  //iterate over alphabet
	  for (int j=0;j<4;j++)
	    {
	      //new score for letter j in A/B
	      int snewA = mA[vi[i]*4+j];
	      int snewB = (ipos >=0) ? mB[ipos*4+j] : 0;
	      //probability of the letter
	      double p = vbg[j];
	      //iterate over old indices A
	      for (int ioldA=0;ioldA<noldrangeA;ioldA++)
		{
		  //score for A
		  int sA = ioldA+noldloA+snewA;
		  //is it a certain hit
		  if (ioldA==noldrangeA-1)
		    sA = nupA+1;
		  //check whether this score has to be considered
		  if (sA>=nloA)
		    {
		      //if it is a certain hit set to max value
		      if (sA>nupA)
			sA = nupA+1;
		      //iterate over old indices B
		      for (int ioldB=0;ioldB<noldrangeB;ioldB++)
			{
			  double pold = Q[ishift][iold][ioldA*noldrangeB+ioldB];
			  if (pold>0)
			    {
			      //score for B
			      int sB = ioldB+voldloB[ishift]+snewB;
			      //was it a certain B hit?
			      if (ioldB==noldrangeB-1)
				sB = nupB+1;
			      //do we have to consider this hit?
			      if (sB>=nloB)
				{
				  //if it is a certain hit set to max value
				  if (sB>nupB)
				    sB = nupB+1;
				  //now we add the prob for this score to Q
				  Q[ishift][inew][(sA-nloA)*nrangeB+(sB-nloB)] += pold*p;
				}
			    }
			}
		    }
		}
	    }
	  //save bounds for next position
	  voldupB[ishift] = nupB;
	  voldloB[ishift] = nloB;
	  //swap indices
	  vold[ishift] = 1-iold;
	  vnew[ishift] = 1-inew;
	}
      //save bounds for next iteration
      noldupA = nupA;
      noldloA = nloA;
    }
  //merge scores=t and scores>t
  CConvolution oconv;
  mygamma = vector<double>(nlen,0);
  for (int ishift=0;ishift<nlen;ishift++)
    {
      int iact = vold[ishift];
      int nrangeB=voldupB[ishift]-voldloB[ishift]+2;
      //compute convolution
      if (ishift>0)
	oconv.convolute(mB,vbg,nlen-ishift);
      double nsum=0;
      for (int iB=0;iB<nrangeB;iB++)
	{
	  //Q is only a 2 x nrangeB matrix
	  double pold=Q[ishift][iact][iB] + Q[ishift][iact][nrangeB+iB];
	  double p = (iB<nrangeB-1) ? oconv.pvalue(tB-(iB+voldloB[ishift])) : 1;
	  mygamma[ishift] += pold*p;
	}
    }
  //adjust gamma
  for (int ishift=0;ishift<nlen;ishift++)
    mygamma[ishift] /= alpha;
}

//prints Q matrix
// parameter:
//   Qsub: Q submatrix
//   nrangeA: range of A scores
//   nrangeB: the same for B
void CCountPars::print_Qsub(vector<double> &Qsub,int nrangeA,int nrangeB)
{
  cout << "*** Qsub ***" << endl;
  for (int jA=0;jA<nrangeA;jA++)
    {
      for (int jB=0;jB<nrangeB;jB++)
	{
	  int i = jA*nrangeB+jB;
	  double p = Qsub[i];
	  cout << "\t(" << jA;
	  if (jA==nrangeA-1)
	    cout << "*";
	  cout << "," << jB;
	  if (jB==nrangeB-1)
	    cout << "*";
	  cout << "):" << p;
	}
      cout << endl;
    }
  cout << "### Qsub ###\n";
} 
