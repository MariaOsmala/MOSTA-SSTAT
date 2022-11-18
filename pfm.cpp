#include <cstdlib>
#include <cstring>
#include "pfm.h"

using namespace std;

//internal function declarations:
class ELineNotFound{};

string mygetline(ifstream &f,char* sstart) throw (ELineNotFound);

//change gc content of PFM
// parameter:
//   igc: new gc content
void CPfm::reinit(double igc, bool abregularize)
{
  //initialize variables
  bregularize = abregularize;
  gc = igc;
  alpha = -1;
  beta = -1;
  bt = 0;
  ic_over = -1;
  n_over=-1;

  //calculate background distribution
  vbg.reserve(4);
  vbg.push_back((1-gc)/2);
  vbg.push_back(gc/2);
  vbg.push_back(gc/2);
  vbg.push_back((1-gc)/2);

  //create PSSM
  pcm2pssm(bregularize);
}

//my fancy constructor
// parameter:
//   asfn: filename
//   agc: gc content
//   at: threshold
CPfm::CPfm(string asfn,double agc,bool abregularize) : bregularize(abregularize),alpha(-1),beta(-1),bt(0),ic(-1),ic_over(-1),n_over(-1),nseqs(0)
{
  //read file
  ifstream f;
  f.open(asfn.c_str());
  if (!f)
    {
      cerr << "Error in pfm.cpp::CPfm: Cannot open file " 
	   << asfn 
	   << " for input.\n";
      exit(1);
    }
  read_pcm(f);
  f.close();
  nlen = mpcm.size()/4;
  //get number of sequences
  for (int i=0;i<nlen;i++)
    for (int j=0;j<4;j++)
      nseqs += mpcm[i*4+j];
  nseqs /= nlen;
  //initialize
  reinit(agc,bregularize);
}

//my fancy constructor based on a filestream
// parameter:
//   f: input filestream
//   agc: gc content
//   at: threshold
CPfm::CPfm(ifstream &f,double agc,bool abregularize) : bregularize(abregularize),alpha(-1),beta(-1),bt(0),ic(-1),ic_over(-1),n_over(-1),nseqs(0)
{
  gc = agc;
  //calculate background distribution
  vbg.reserve(4);
  vbg.push_back((1-gc)/2);
  vbg.push_back(gc/2);
  vbg.push_back(gc/2);
  vbg.push_back((1-gc)/2);
  //read file
  read_pcm(f);
  nlen = mpcm.size()/4;
  //get number of sequences
  for (int i=0;i<nlen;i++)
    for (int j=0;j<4;j++)
      nseqs += mpcm[i*4+j];
  nseqs /= nlen;
  //create PSSM
  pcm2pssm(bregularize);
}

//another constructor to merge to PFMs
// parameter:
//   opfmA: first Pfm
//   opfmB: second Pfm
//   ishift: shift
//   bcompl: take complement?
CPfm::CPfm(CPfm& opfmA,CPfm& opfmB,int ishift,bool bcompl)
{
  int nA = opfmA.nlen;
  int nB = opfmB.nlen;
  int nnew = (ishift<0) ? nA-ishift : nB+ishift;
  //compute background vectors
  vector<int> vbgA(4,0);
  vector<int> vbgB(4,0);
  for (int j=0;j<4;j++)
    {
      vbgA[j] = (int)(opfmA.vbg[j]*opfmA.nseqs);
      vbgB[j] = (int)(opfmB.vbg[j]*opfmB.nseqs);
    }

  //compute new mpcm
  mpcm.reserve(nnew*4);
  //consider reverse complement
  vector<double> mpcm_add;
  if (bcompl)
    {
      pfm_reversecompl(opfmB.mpcm,mpcm_add);
    }
  else 
    {
      mpcm_add = opfmB.mpcm;
    }
  
  //remember overlapping part
  vector<bool> boverlap;
  boverlap.reserve(nnew);

  //left flank
  for (int i=ishift;i<0;i++)
    {
      for (int j=0;j<4;j++)
	mpcm.push_back(mpcm_add[(i-ishift)*4+j]+vbgA[j]);
      boverlap.push_back(false);
    }
  //remember overlapping part
  int istartover = -1;
  n_over=0;
  //overlapping part
  for (int i=0;i<nA;i++)
    {
      int iposB = i-ishift;
      if ((iposB>=0) && (iposB<nB))
	{
	  //remember start position of overlap
	  if (istartover<0)
	    istartover = boverlap.size();
	  //fill PCM
	  for(int j=0;j<4;j++)
	    {
	      mpcm.push_back(opfmA.mpcm[i*4+j]+mpcm_add[iposB*4+j]);
	    }
	  n_over++;
	  boverlap.push_back(true);
	}
      else
	{
	  for(int j=0;j<4;j++)
	    mpcm.push_back(opfmA.mpcm[i*4+j]+vbgB[j]);
	  boverlap.push_back(false);
	}
      
    }
  //right flank
  for (int i=nA-ishift;i<nB;i++)
    {
      for (int j=0;j<4;j++)
	mpcm.push_back(mpcm_add[i*4+j]+vbgA[j]);
      boverlap.push_back(false);
    }
  
  //initialize variables
  sid = opfmA.sid+string(":")+opfmB.sid;
  nlen = mpcm.size()/4;
  gc = opfmA.get_gc();
  //get number of sequences
  nseqs = 0;
  for (int i=0;i<nlen;i++)
    for (int j=0;j<4;j++)
      nseqs += mpcm[i*4+j];
  nseqs /= nlen;
  //calculate background distribution
  vbg.reserve(4);
  vbg.push_back((1-gc)/2);
  vbg.push_back(gc/2);
  vbg.push_back(gc/2);
  vbg.push_back((1-gc)/2);
  //create PSSM
  pcm2pssm(opfmA.bregularize);

  //get overlapping ic
  ic_over = get_ic(istartover,n_over);

  //set alpha to -1 because threshold has to be adjusted first
  alpha = -1;
  beta = -1;
  bt = 0;
  ic = -1;
}

//returns gc content
double CPfm::get_gc()
{
  return(gc);
}

//adjusts threshold
// parameter:
//   smethod: 'balanced', 'typeI', 'typeII', 'typeIext'
//   tp: parameter for the method
void CPfm::adjust_t(string smethod,double tp) throw(EThresholdNotSet)
{
  alpha = -1;
  beta = -1;
  if ((smethod!="balanced") && (tp==-23880))
    throw EThresholdNotSet();
  CConvolution oconv;
  if (smethod=="threshold")
    {
      t = (int)tp;
    } else if (smethod == "typeI")
    {      
      for (int i=0;i<nlen;i++)
	oconv.convolute(mpssm,vbg,i);
      t = oconv.threshold(-log(1-tp)/500);
      alpha = oconv.pvalue(t);
    } else if (smethod == "typeII")
    {
      for (int i=0;i<nlen;i++)
	{
	  vector<double> v(4,0);
	  for (int j=0;j<4;j++)
	    v[j] = mpwm[i*4+j];
	  oconv.convolute(mpssm,v,i);
	}
      t = oconv.threshold(1-tp);
      beta = 1-oconv.pvalue(t);
    } else if (smethod == "balanced")
    {
      for (int i=0;i<nlen;i++)
	oconv.convolute(mpssm,vbg,i);
      CConvolution oconvB;
      for (int i=0;i<nlen;i++)
	{
	  vector<double> v(4,0);
	  for (int j=0;j<4;j++)
	    v[j] = mpwm[i*4+j];
	  oconvB.convolute(mpssm,v,i);
	}
      t = oconv.balanced(oconvB,500);
      alpha = oconv.pvalue(t);
      beta = 1-oconvB.pvalue(t);
    } else if (smethod == "typeIext")
    {
      for (int i=0;i<nlen;i++)
	oconv.convolute(mpssm,vbg,i);
      CConvolution oconvB;
      for (int i=0;i<nlen;i++)
	{
	  vector<double> v(4,0);
	  for (int j=0;j<4;j++)
	    v[j] = mpwm[i*4+j];
	  oconvB.convolute(mpssm,v,i);
	}
      t = oconv.balanced(oconvB,500);
      if (oconv.pvalue_region(t,500)>tp)
	t=oconv.threshold(-log(1-tp)/500);
      alpha = oconv.pvalue(t);
      beta = 1-oconvB.pvalue(t);
    } else throw EThresholdNotSet();
  bt=1;
}

double CPfm::get_alpha(int m) throw (EThresholdNotSet)
{
  if (!bt)
    throw EThresholdNotSet();
  double alphalocal;
  if (alpha<0)
    {
      CConvolution oconv;
      for (int i=0;i<nlen;i++)
	oconv.convolute(mpssm,vbg,i);
      alphalocal = oconv.pvalue(t);
      alpha = alphalocal;
    } else 
    {
      alphalocal = alpha;
    }
  if (m>0)
    {
      alphalocal = 1-exp(-m * alpha);
    }
  return(alphalocal);
}

double CPfm::get_beta() throw (EThresholdNotSet)
{
  if (!bt)
    throw EThresholdNotSet();
  if (beta<0)
    {
      CConvolution oconvB;
      for (int i=0;i<nlen;i++)
	{
	  vector<double> v(4,0);
	  for (int j=0;j<4;j++)
	    v[j] = mpwm[i*4+j];
	  oconvB.convolute(mpssm,v,i);
	}
      beta = oconvB.pvalue(t);
    }
  return(beta);
}

int CPfm::get_t() throw (EThresholdNotSet)
{
  if (!bt)
    throw EThresholdNotSet();
  return(t);
}

// reads transfac file
// parameter:
//   sfn: filename
void CPfm::read_pcm(ifstream &f)
{
  //Search for the ID tag
  try
    {
      sid = mygetline(f,"ID");
    } catch (ELineNotFound) 
    {
      cerr << "Error in read_pcm: Cannot find line starting with ID.\n";
    }
  //Search for the P0 Tag
  string sp0;
  try
    {
      sp0 = mygetline(f,"P0");
    } catch (ELineNotFound) 
    {
      cerr << "Error in read_pcm: Cannot find line starting with P0.\n";
      exit(1);
    }
  //read content of matrix
  string s;
  int i=1;
  while (i>0)
    {
      getline(f,s);
      if (s.substr(0,2)=="XX")
	{
	  if (i==0)
	    {
	      cerr << "Error in read_pcm: Cannot find first line 01 of transfac file.\n";
	      exit(1);
	    } else i = -2;
	}
      if (i>0)
	{
	  //split into vector and assign to matrix
	  istringstream istr(s);
	  for (int j=0;j<5;j++)
	    {
	      double k;
	      istr >> k;
	      //skip line number
	      if (j>0)
		mpcm.push_back(k);
	    }
	}
      i++;
    }
  //read until next '\\' or eof
  while (!(f.eof() || (s.substr(0,2)=="//")))
    getline(f,s);
}

void CPfm::pcm2pssm(bool bregularize)
{
  //some initialization
  double epsilon = .05;
  //not very nice but otherwise I get a compiler error...
  vector<long double> mpwma(mpcm.size(),0);
  mpwm = mpwma;
  //use Sven's method or standard?
  if (bregularize)
    {
      //prepare the score distribution (summarized over all positions)
      vector<double> vc(4,0);
      vector<long double> vp(4,0);
      double ncount=0;
      double add=0;
      for (int i = 0; i < nlen; i++) 
	{
	  for (int j = 0; j < 4; j++) 
	    {
	      vc[j] += mpcm[i*4+j];
	    }
	}
      ncount = vc[0] + vc[1] + vc[2] + vc[3];
      if (vc[0] == 0 || vc[1] == 0 || vc[2] == 0 || vc[3] == 0) 
	{
	  add = 1;
	}
      
      for (int j = 0; j < 4; j++) 
	{
	  vp[j] = ( vc[j] + add/4.0 ) / ( ncount + add );
	}
      
      //regularize the matrix
      
      double prec = 0.0001;
      double cutoff = 1.5;
      
      // calculate raw_p_mat
      vector<long double> raw_p_mat(mpcm.size(),0);
      
      for (int i = 0; i < nlen; i++) 
	{
	  double nc = mpcm[i*4] + mpcm[i*4 + 1] + mpcm[i*4 + 2] + mpcm[i*4 + 3];
	  if (nc == 0) 
	    {
	      for (int j = 0; j < 4; j++) 
		{
		  raw_p_mat[i*4+j] = vp[j];
		}
	    } else 
	    {
	      for (int j = 0; j < 4; j++) 
		{
		  raw_p_mat[i*4+j] = (long double)mpcm[i*4 + j]/nc;
		}
	    }
	}

      //regularize
      for (int i = 0; i < nlen; i++) 
	{
	  long double delta = 0;
	  double nc = 0;
	  for (int j = 0; j < 4; j++) 
	    {
	      if (mpcm[i*4 + j] > 0) 
		{
		  delta += raw_p_mat[i*4+j] * log(raw_p_mat[i*4+j]/vp[j]);
		}
	      nc += mpcm[i*4+j];
	    }
	  delta *= 2 * nc;
	  
	  if (delta <= cutoff) 
	    {
	      for (int j=0;j<4;j++)
		mpwm[i*4+j] = vp[j];
	    } else 
	    {
	      long double deltaw = delta;
	      long double w = 0.5;
	      long double v = 0.25;
	      while (abs(delta - cutoff - deltaw) > prec) 
		{
		  deltaw = 0;
		  for (int j = 0; j < 4; j++) 
		    {
		      deltaw += ((1 - w)*raw_p_mat[i*4 + j] + w*vp[j]) * log(((1 - w)*raw_p_mat[i*4+j] + w*vp[j])/vp[j]);
		    }
		  deltaw *= 2 * nc;
		  if (deltaw >= delta - 1.5) 
		    {
		      w += v;
		    } else 
		    {
		      w -= v;
		    }
		  v *= 0.5;
		}
	      for (int j = 0; j < 4; j++) 
		{
		  mpwm[i*4+j] = (1 - w)*raw_p_mat[i*4+j] + w*vp[j];
		}
	    }
	}
    } else 
    {
      double xreg = 0.01;
      for (int i = 0; i < nlen; i++) 
	{
	  double nc =0;
	  //get sum per position
	  for (int j = 0; j < 4; j++) 
	    {
	      nc += mpcm[i*4+j];
	    }
	  //get weights
	  for (int j = 0; j < 4; j++) 
	    {
	      mpwm[i*4+j] = ((double)mpcm[i*4+j]/(double)nc+xreg)/(1+4*xreg);
	    }
	}
    }
  
  //calculate scoring matrix
  mpssm.clear();
  mpssm.reserve(mpcm.size());
  for (int i = 0; i < nlen; i++) 
    {
      for (int j = 0; j < 4; j++) 
	{
	  long double x= log(mpwm[i*4+j]/vbg[j])/epsilon;
	  //correct rounding
	  long double xadd=.5;
	  if (x<0)
	    xadd *= -1;
	  mpssm.push_back((int)(x+xadd));
	}
    }
}

//returns information content average over columns
double CPfm::get_ic(int istart,int nr)
{
  double iclocal=0;
  if ((istart==0) && (nr<0) && (ic>=0))
    {
      //return ic from alst computation
      iclocal=ic;
    } else
    {
      //init region
      if (nr<0)
	nr = nlen-istart;

      iclocal = 0;
      //iterate over positions
      for (int i=istart;i<istart+nr;i++)
	{
	  //iterate over bases
	  for (int j=0;j<4;j++)
	    {
	      if (mpwm[i*4+j]>0)
		iclocal += mpwm[i*4+j] * log(mpwm[i*4+j]/vbg[j]);
	    }
	}
      //save if ic from whole PFM
      if ((istart==0) && (nr<0))
	ic = iclocal;
    }
  return(iclocal);
}

//funtion to printout PFM in transfac format
ostream& operator<<(ostream& strm,const CPfm& obj)
{
  strm << "ID\t" << obj.sid << endl;
  strm << "XX\nP0\tA\tC\tG\tT\n";
  for (int i=0;i<(obj.mpcm.size()/4);i++)
    {
      strm << i;
      for (int j=0;j<4;j++)
	strm << "\t" << obj.mpcm[4*i+j];
      strm << endl;
    }
  strm << "XX\n//\n";
  return strm;
}

string mygetline(ifstream &f,char* sstart) throw (ELineNotFound)
{
  string s;
  while (getline(f,s) && (s.substr(0,2)!=sstart)) {}
  if (s.substr(0,2)!=sstart)
    throw ELineNotFound();
  string sout(s.substr(strlen(sstart)+1));
  chomp(sout);
  return(sout);
}
