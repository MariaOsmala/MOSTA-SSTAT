
#include "sequences.h"

using namespace std;

//initialize static variables

//create map for accepted letters
pair<char, char> a[] =
{
pair<char, char>('A', 'A'),
pair<char, char>('C', 'C'),
pair<char, char>('G', 'G'),
pair<char, char>('T', 'T'),
pair<char, char>('N', 'N'),
pair<char, char>('X', 'N'),
};
map<char, char> CSequence::dtrans(a, a + sizeof(a) / sizeof(a[0]));

//create map for complement
pair<char, char> b[] =
{
pair<char, char>('A', 'T'),
pair<char, char>('C', 'G'),
pair<char, char>('G', 'C'),
pair<char, char>('T', 'A'),
pair<char, char>('N', 'N'),
};
map<char, char> CSequence::dcompl(b, b + sizeof(b) / sizeof(b[0]));

//definition of base transformation
// parameter:
//   cc: character to transform
//   brev: use the complement?
char CSequence::basetrans(char cc,bool brev) throw(EWrongFormat)
{
  char c = toupper(cc);
  //make map to get complement letters
 
  char cr='Z';

  if ((c!=' ') && (c!='-') && (c!='\\'))
    {
      map<char, char>::iterator cur  = dtrans.find(c);

      //is letter allowed?
      if (cur == dtrans.end())
	{
	  stringstream serr("");
	  serr << "Sequence contains a prohibited letter: ";
	  serr << "'" << c << "'.\n";
	  throw(EWrongFormat(serr.str()));
	}
      cr = (*cur).first;
      
      if (brev)
	cr = dcompl[cr];
    }
  return(cr);
}

CSequence::CSequence(string& asid,const string& asseq) throw(EWrongFormat) : sseq(""),srseq(""),sid(asid),gc(0)
{
  
  //get gc content
  long int ngc=0;
  long int nat=0;

  //transform strings
  for (int i=0;i<asseq.length();i++)
    {
      
      char c = basetrans(asseq[i],false);
      if (c!='Z')
	sseq += c;

      char cr = basetrans(asseq[asseq.length()-i-1],true);
      if (c!='Z')
	srseq += cr;
      
      if ((c=='G') || (c=='C'))
	ngc++;
      if ((c=='A') || (c=='T'))
	nat++;
    }

  gc = (double)ngc/(ngc+nat);
}

string CSequence::get_sequence(bool breverse) const
{
  if (breverse)
    {
      return(srseq);
    } else 
    {
      return(sseq);
    }
}

double CSequence::get_gc()
{
  return(gc);
}

string CSequence::get_sid() const
{
  return(sid);
}

//class CSequences

CSequences::CSequences(char* sf) throw(EFileNotFound,EWrongFormat) : fin(sf),beof(true)
{
  //could we open the file?
  if (!fin)
    throw EFileNotFound(string("CSequences::CSequences: File "+string(sf)+" could not be opened.\n"));
  
  //read first id
  string s;
  if (getline(fin,s))
    {
      if (s.substr(0,1)!=">")
	throw EWrongFormat(string("CSequences::CSsequences: File "+string(sf)+" is not in correct FASTA format - not starting with '>'.\n"));
      //name of next gene
      sid = s.substr(1);
      chomp(sid);
      beof=false;
    } else throw EWrongFormat(string("CSequences::CSequences: File "+string(sf)+" seems to be empty.\n"));  
}

CSequences::~CSequences()
{
  fin.close();
}

bool CSequences::eof()
{
  return(beof);
}

//make CSequence object containing next sequence
CSequence CSequences::next()
{
  //read next sequence from file
  string s;
  stringstream sseq (stringstream::in | stringstream::out);
  while (getline(fin,s) && (s.substr(0,1)!=">"))
    {
      chomp(s);
      sseq << s;
    }

  //create sequence object
  CSequence o(sid,sseq.str());

  //check whether we are at the end of the file
  if (s.substr(0,1)==">")
    {
      sid = s.substr(1);
      chomp(sid);
    } else 
    {
      beof = true;
    }
  return(o);
}

//funtion to printout sequence in fasta format
ostream& operator<<(ostream& strm,const CSequence& obj)
{
  strm << ">" << obj.sid << endl;
  strm << obj.get_sequence() << endl;
  return strm;
}
