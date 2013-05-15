#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include "capmon/feature.h"
#include "capmon/gff_validation.h"
#include <fstream>
#include "boost/regex.hpp"
#include "boost/regex.hpp"
#include <stdexcept>
#include <sys/stat.h> 
#include<time.h>   

/*
* g++ -std=c++0x capmon/gff_validation.cpp -O2 -c -o  gff_validation_mini.o
* g++ -Wall -Wextra -ansi -pedantic -Werror -g -static-libstdc++ -O -std=c++0x converter.cpp gff_validation_mini.o -lboost_regex -o coord
*/

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::runtime_error;

#define EXEC "coordthingy"
#define VERSION "0.0.0"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

const static char * CORRECT_USAGE = 
"Exec: " EXEC
"\nVersion: " VERSION " "
"\nBuild: " __DATE__ " @ " __TIME__
"\nWith gcc: " __VERSION__
"\nWith glibc: "  TOSTRING(__GLIBC__) "."  TOSTRING(__GLIBC_MINOR__)
"\n"
"\nDaniel S. T. Hughes. dsth@ebi.ac.uk, dsth@cantab.net\n"
"\nUsage: ./" EXEC " [COMMAND]\n"
"\n\tgenome2trans\t\tconvert..."
"\n\tgenome2cds\t\tconvert..."
"\n\tgenome2pep\t\tconvert..."
"\n\ttrans2genome\t\tconvert..."
"\n\ttrans2trans\t\tconvert...";

typedef unsigned char uchar;
typedef unsigned int uint;

typedef boost::regex regex;
typedef boost::smatch smatch;

typedef struct _VCF { // no need?!?
    string chr;
    uint pos;
    // string id;
    // char ref;
    string ref; // Each base must be one of A,C,G,T,N. Bases should be in uppercase. Multiple bases are permitted. The value in the POS field refers to the position of the first base in the String.
    string alt; // ALT comma separated list of alternate non-reference alleles called on at least one of the samples. Options are base Strings made up of the bases A,C,G,T,N, or an angle-bracketed ID String (”<ID>”).
    // qual, filter, info
} VCF;

typedef struct {
    string trans;
    uint start; // 0-base-start
    uint end;
} BED;

struct feature_named_min : public feature_min {
    // explicit feature_ext(const feature_min& fm) : feature_min(fm) {}
    feature_named_min() = delete;
    feature_named_min(const feature_min& fm) = delete;
    // ~feature_ext () {}
    feature_named_min(std::string n, uint a, uint b ) : feature_min(a,b), _name(n) {}
    std::string id() const { return _name; } 
  private :
    std::string _name;
};

typedef std::pair<bool,std::vector<feature_named_min>> feat_vec;

/// genomic coords always absolute/transcript strand dependent?!?

//    std::ifstream in(filename);
//        if(in == 0) throw runtime_error("problem opening gff file : " + std::string(filename));
//
static regex bed_regex("([^#\\t]+)\\t(\\d+)\\t(\\d+).*");
static regex vcf_regex("([^#\\t]+)\\t(\\d+)\\t[^\\t]+\\t([^\\t]+)\\t([^\\t]+)\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+.*");

std::ifstream& getvcf(std::ifstream& infs, VCF* vcf) { // prolly should be externally visible?!?

    smatch match_obj; // should prolly allocate just the once?!?
    string str;
    // delim
    // The delimiting character. The operation of extracting successive characters is stopped when this character is read.
    // char widen ( char c ) const; Widen character // Converts c to its locale equivalent.
    
    //// duh - need to itirate through lines till match!?!
    
    while(getline(infs, str, '\n')) {

        //cout << "the line"<<str<<"\n";
        if(regex_match(str,match_obj,vcf_regex)) {
        // cout << "we have a match"<<str<<"\n";
            vcf->chr=match_obj[1];
            vcf->pos=std::stoi(std::string(match_obj[2]));
            vcf->ref=match_obj[3];
            vcf->alt=match_obj[4];
            return infs;
        }
    } 
    
    return infs;

}

std::ifstream& getbed(std::ifstream& infs, BED* bed) { // prolly should be externally visible?!?
    smatch match_obj; // should prolly allocate just the once?!?
    string str;
    while(getline(infs, str, '\n')) {
        if(regex_match(str,match_obj,bed_regex)) {
            bed->trans=match_obj[1];
            bed->start=std::stoi(std::string(match_obj[2]));
            bed->end=std::stoi(std::string(match_obj[3]));
            return infs;
        }
    } 
    return infs;
}

// template<typename _CharT, typename _Traits, typename _Alloc>
//    inline basic_istream<_CharT, _Traits>&
// getline(basic_istream<_CharT, _Traits>& __is,
// basic_string<_CharT, _Traits, _Alloc>& __str)
// { return getline(__is, __str, __is.widen('\n')); }
// 
//  template<>
//      basic_istream<char>&
//          getline(basic_istream<char>& __in, basic_string<char>& __str,
    //                      char __delim);
//


class Object {
    struct ObjectConcept {  
        virtual ~ObjectConcept() {} 
        virtual bool has_concept_of_some_property() const = 0;
        virtual std::string name() const = 0;
    };

    template<typename T> struct ObjectModel : ObjectConcept {
        ObjectModel( const T& t ) : object( t ) {}
        virtual ~ObjectModel() {} // not needed it's already virtual in base?!?
        virtual bool has_concept_of_some_property() const { return object.has_some_property_blah(); } 
        virtual std::string name() const { return typeid(object).name(); }
      private:
        T object; // the actual stored type
    };
    std::shared_ptr<ObjectConcept> object; 
  public:
    template<typename T> Object(const T& obj) : object( new ObjectModel<T>( obj ) ) {}
    std::string name() const { return object->name(); }
    bool has_concept_of_some_property() const { return object->has_concept_of_some_property(); }
};

struct feature_converter : public feature_min {

//// have a conversion ctor from feature_min!?!?!

    explicit feature_converter(const feature_min& fm) : feature_min(fm) {}

    feature_converter() = delete;
    ~feature_converter () {}
    feature_converter(uint a, uint b) : feature_min(a,b), tstart(0), tend(0) {} // just forward args and default initialise?!?
    // feature_min(uint a, uint b) : _gstart(a), _gend(b), tstart(0), tend(0) {} 
    uint tstart;
    uint tend;
};

// should prolly keep transcript start/stop - i.e. to make sure the transcript has exons spanning full extension?!?
class transcript {

    transcript()=delete;
  private:

    uchar _strand;
    uint _gpstart;
    uint _gpend;
    uint _peptide_transcript_offset;
int _peptide_transcript_length;
    std::vector<feature_converter> exons;
    /// seems horrible to store cds too but alternative it to iterate through all of them for all transcritps for no reason?!?
    std::vector<feature_converter> cds;

  public:

    transcript(uchar s) : _strand(s), _gpstart(0), _gpend(0), _peptide_transcript_offset(0), _peptide_transcript_length(0) { if (_strand!=1&&_strand!=0) 
      throw runtime_error("must be 1 or 0"); };

    void addexon(const feature_converter& f) { exons.push_back(std::move(f)); }
    void addcds(const feature_converter& f) { cds.push_back(std::move(f)); }

    void printexons() {
        for_each(exons.begin(),exons.end(), [](feature_converter& f) { cout << "transcript: "<< f.gstart() <<" "<<f.gend() << " " << f.tstart << " " << f.tend<< "\n"; });
    }

    int convertt2g(uint i) { // not using negative values as message so uint fine - but all a bit risky

        if (!exons.at(0).tstart) { // initialised to 0?
            _sortexons();
            _fill_in_coords();

            if(cds.size()>0) _fill_in_peptide();
        } 

        // cout <<"will convert " << i<<"\n";
        // int tmp=1;
        for(auto it=exons.begin(); it!=exons.end(); it++/*,++tmp*/) {

        // cout << "checking i="<<i<<" against "<<it->tstart<<" "<<it->tend<<"\n";
            if(it->tstart<=i&&it->tend>=i) {
                // cout << "we have our exon " << it->tstart << "-"<<it->tend<<"\n";
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
//     cout << " (exon:"<<tmp<<") ";
                if (_strand) return it->gstart()+i-it->tstart;
                else return it->gend()-i+it->tstart;
            }


        }
        
        throw runtime_error("argument must fall within range of transcript");
        
    }

    // silly uint convertg2t(uint i) {
    int convertg2t(uint i) {

        /// should be factored - clrearly only invoked for first mapping call to any transcript - peptide stage not necessarily essential?!?
        if (!exons.at(0).tstart) {
            _sortexons();
            _fill_in_coords();

            if(cds.size()>0) _fill_in_peptide();
            // do we fill in cds details here?!? - treat them in the same way or just give exons extra parameter?!?
            // and/or phase to allow them to deal with peptide location?!?
            //
            // just find peptide absolute genomic positions?!?
            // do we do any atual qc on the models?!?
        } 


        // cout <<"will convert " << i<<"\n";
        /// how to bother with introns?!? - could be horribly lazy and use two loops?!?
        for(auto it=exons.begin(); it!=exons.end(); it++) {

       //   cout << "checking i="<<i<<" against "<<it->gstart()<<" "<<it->gend()<<"\n";
            if(it->gstart()<=i&&it->gend()>=i) {
                // cout << "we have our exon " << it->tstart << "-"<<it->tend<<"\n";
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
                if (_strand) return i-it->gstart()+it->tstart;
                else return it->gend()-i+it->tstart;
            }


        }

        auto it2=exons.begin();
        std::vector<feature_converter>::iterator it1;
        int intron = 0;
        while (1) {

            // auto it1=it2++;
            it1=it2++;
            intron++;
            if (it2==exons.end()) break;
// cout << "using "<<it2->gend()+1 <<" - "<<it1->gstart()-1 << "\n";
        
            if((_strand&&it1->gend()+1<=i&&it2->gstart()-1>=i)||(!_strand&&it2->gend()-1<=i&&it1->gstart()+1>=i)) {
                cout << "# g2t : [it's in intron " << intron << "] : " << i;// : "; //  << 
                // it1->gend()+1 << "-" << it2->gstart()-1<<"\n";
                // _strand?(it1->gend()+1):(it2->gend()-1) << "-" << _strand?(it2->gstart()-1):(it1->gstart()+1) << "\n";
                return 0;
                // return 0;
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
            }
        }
        
        throw runtime_error("argument must fall within range of transcript : " + std::to_string(i));

    }

    //r is there any point in peptide->genome?!?
    //r these prolly ought to be handled with overloading or a switch parameter and not separate function calls?!?

    int convertg2cds(uint i) {
    // duh: uint convertg2cds(uint i) {

        //r can either do same mapping pairs but for cds or can convert to transcript coords and adjust with offset/length?!?

        /// should be factored - clrearly only invoked for first mapping call to any transcript - peptide stage not necessarily essential?!?
        if (!exons.at(0).tstart) {
            _sortexons();
            _fill_in_coords();
            if(cds.size()>0) _fill_in_peptide();
            // just find peptide absolute genomic positions?!?
        } 

        // cout <<"will convert " << i<<"\n";
        /// how to bother with introns?!? - could be horribly lazy and use two loops?!?
        for(auto it=exons.begin(); it!=exons.end(); it++) {

       //   cout << "checking i="<<i<<" against "<<it->gstart()<<" "<<it->gend()<<"\n";
            if(it->gstart()<=i&&it->gend()>=i) {
                // cout << "we have our exon " << it->tstart << "-"<<it->tend<<"\n";
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
                if (_strand) {
                    //int ti = i-it->gstart()+it->tstart;
                    // clearly can do in terms of transcript or pepptide pos?!?
                    int tp = i-it->gstart()+it->tstart-_peptide_transcript_offset;
                    //cout << "\ntrans="<<ti<<"\npep="<<tp<<"\noffset="<<_peptide_transcript_offset<<"\nlength="<<_peptide_transcript_length<<"\n";
                    if(tp<0) {
                        cout << "# g2c : [it's in 5' utr] ";
                        return -5;
                    } else if (tp>_peptide_transcript_length) {
                        cout << "# g2c : [it's in 3' utr] ";
                        return -3;
                    } else {
                        return tp;
                    }
                } else { 
                    // int ti = it->gend()-i+it->tstart;
                    int tp = it->gend()-i+it->tstart-_peptide_transcript_offset;
                    // cout << "\ntrans="<<ti<<"\npep="<<tp<<"\noffset="<<_peptide_transcript_offset<<"\nlength="<<_peptide_transcript_length<<"\n";
                    if(tp<0) {
                        cout << "# g2c : [it's in 5' utr] ";
                        return -5;
                    } else if (tp>_peptide_transcript_length) {
                        cout << "# g2c : [it's in 3' utr] ";
                        return -3;
                    } else {
                        return tp;
                    }
                }
            }


        }

        /// do we care about this?!?
        auto it2=exons.begin();
        std::vector<feature_converter>::iterator it1;
        int intron = 0;
        while (1) {

            // auto it1=it2++;
            it1=it2++;
            intron++;
            if (it2==exons.end()) break;
// cout << "using "<<it2->gend()+1 <<" - "<<it1->gstart()-1 << "\n";
        
            if((_strand&&it1->gend()+1<=i&&it2->gstart()-1>=i)||(!_strand&&it2->gend()-1<=i&&it1->gstart()+1>=i)) {
                cout << "# g2c : [it's in intron " << intron << "] : " << i;// : "; //  << 
                // it1->gend()+1 << "-" << it2->gstart()-1<<"\n";
                // _strand?(it1->gend()+1):(it2->gend()-1) << "-" << _strand?(it2->gstart()-1):(it1->gstart()+1) << "\n";
                return 0;
                // return 0;
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
            }
        }
        
        throw runtime_error("argument must fall within range of transcript : " + std::to_string(i));

    }

    int convertg2p(uint i) {
    // duh: uint convertg2cds(uint i) {

        //r can either do same mapping pairs but for cds or can convert to transcript coords and adjust with offset/length?!?

        /// should be factored - clrearly only invoked for first mapping call to any transcript - peptide stage not necessarily essential?!?
        if (!exons.at(0).tstart) {
            _sortexons();
            _fill_in_coords();
            if(cds.size()>0) _fill_in_peptide();
            // just find peptide absolute genomic positions?!?
        } 

        // cout <<"will convert " << i<<"\n";
        /// how to bother with introns?!? - could be horribly lazy and use two loops?!?
        for(auto it=exons.begin(); it!=exons.end(); it++) {

       //   cout << "checking i="<<i<<" against "<<it->gstart()<<" "<<it->gend()<<"\n";
            if(it->gstart()<=i&&it->gend()>=i) {
                // cout << "we have our exon " << it->tstart << "-"<<it->tend<<"\n";
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
                if (_strand) {
                    //int ti = i-it->gstart()+it->tstart;
                    // clearly can do in terms of transcript or pepptide pos?!?
                    int tp = i-it->gstart()+it->tstart-_peptide_transcript_offset;
                    //cout << "\ntrans="<<ti<<"\npep="<<tp<<"\noffset="<<_peptide_transcript_offset<<"\nlength="<<_peptide_transcript_length<<"\n";
                    if(tp<0) {
                        cout << "# g2p : [it's in 5' utr] ";
                        return -5;
                    } else if (tp>_peptide_transcript_length) {
                        cout << "# g2p : [it's in 3' utr] ";
                        return -3;
                    } else {
                    //r do we want remainder - i.e. either as base offset within codon OR by using a float?!?
                        // return tp%3;
                        return (tp/3);
                    }
                } else { 
                    // int ti = it->gend()-i+it->tstart;
                    int tp = it->gend()-i+it->tstart-_peptide_transcript_offset;
                    // cout << "\ntrans="<<ti<<"\npep="<<tp<<"\noffset="<<_peptide_transcript_offset<<"\nlength="<<_peptide_transcript_length<<"\n";
                    if(tp<0) {
                        cout << "# g2p : [it's in 5' utr] ";
                        return -5;
                    } else if (tp>_peptide_transcript_length) {
                        cout << "# g2p : [it's in 3' utr] ";
                        return -3;
                    } else {
                        // return tp%3;
                        return (tp/3);
                    }
                }
            }


        }

        /// do we care about this?!?
        auto it2=exons.begin();
        std::vector<feature_converter>::iterator it1;
        int intron = 0;
        while (1) {

            // auto it1=it2++;
            it1=it2++;
            intron++;
            if (it2==exons.end()) break;
// cout << "using "<<it2->gend()+1 <<" - "<<it1->gstart()-1 << "\n";
        
            if((_strand&&it1->gend()+1<=i&&it2->gstart()-1>=i)||(!_strand&&it2->gend()-1<=i&&it1->gstart()+1>=i)) {
                cout << "# g2p : [it's in intron " << intron << "] : " << i;// : "; //  << 
                // it1->gend()+1 << "-" << it2->gstart()-1<<"\n";
                // _strand?(it1->gend()+1):(it2->gend()-1) << "-" << _strand?(it2->gstart()-1):(it1->gstart()+1) << "\n";
                return 0;
                // return 0;
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
            }
        }
        
        throw runtime_error("argument must fall within range of transcript : " + std::to_string(i));

    }
  private:
    // template<typename T> inline void doit(T& it, T& eit) {

    void _sortexons() {
        // could use std::function?!?
        // cout << "we have strand="<<_strand<<"\n";
        if (_strand) sort(exons.begin(),exons.end(),[](feature_converter a, feature_converter b) { return (a.gstart()<b.gstart()); });
        else sort(exons.begin(),exons.end(),[](feature_converter a, feature_converter b) { return (a.gstart()>b.gstart()); });
    }

    void _sortcds() {
        // if (_strand) 
            sort(cds.begin(),cds.end(),[](feature_converter a, feature_converter b) { return (a.gstart()<b.gstart()); });
        // else sort(cds.begin(),cds.end(),[](feature_converter a, feature_converter b) { return (a.gstart()>b.gstart()); });
    }

    void _fill_in_coords() {
        int i=0,last_te=0;
        for (auto it=exons.begin() ; it!=exons.end(); it++,i++) { // for ( ; it!=eit; it++,i++) {
            // cout <<"i="<<i<< " " << it->gstart() << "\n";
            if(i==0) {
                it->tstart=1; // already 0 initialised tstart?!?!
                // use ternaries e.g. it->tend = strand ? it->gend-it->gstart+1 : ...?!?
                it->tend=it->gend()-it->gstart()+1;
            } else {
                it->tstart=last_te+1;
                it->tend=it->gend()-it->gstart()+it->tstart; 

            }
            last_te=it->tend;
        }
    }

    void _fill_in_peptide() {

        // if we're not gonna do model qc no actual reason to order the cds?!?
        _sortcds();
/// this checking should be an option?!?
        /* do we really care about this?!? - just use first and last?!?
        for (auto cit=cds.begin() ; cit!=cds.end(); cit++) { 
            
            cout << "check cds " << cit->gstart() << "-" << cit->gend() << "\n";
            bool found=false;
        
            for (auto eit=exons.begin() ; eit!=exons.end(); eit++) { 
                if (cit->gstart()>=eit->gstart()&&cit->gend()<=eit->gend()) {
                    cout << " heya!!!\n";
                    found=true;
                }
            }
            if (!found) throw std::string("unable to find corresponding exon - this gff model doesn't make sense!?!");
        }
        */

        feature_converter ff = cds.front(), fb = cds.back();

        ////// all we really want is peptide transcript offset and peptide length!?!?
        ///// all else follows from that?!?
        //// offset is trivial - length is then just full transcript length minus other offset!?!?

        // _pstart = cds.front().gstart();
        // _pend = cds.back().gend();
        // sorted 5'->3' 
        _gpstart = ff.gstart();
        _gpend = fb.gend();

        /// can clearly use convertg2t - but it there really any point - i.e. if we have corresponding first/last exons?!?
//cout << "peptide 5' start is "<< _gpstart << "\n";
//cout << "peptide 3' end is "<< _gpend << "\n";

//        printexons();

        for (auto eit=exons.begin() ; eit!=exons.end(); eit++) { 
            ///// do 5'
            if (ff.gstart()>=eit->gstart()&&ff.gend()<=eit->gend()) {
//                cout << "checkin "<<eit->gstart()<<"-"<<eit->gend()<<"\n";
              // always the 5' start - will be transcript start if strand is positive?!?
              if (_strand==1) _peptide_transcript_offset = eit->gend()-_gpstart+1;
              // if (_strand) _peptide_transcript_length=_gpstart-eit->gstart()+1;
              if (_strand==0) _peptide_transcript_length= eit->tend+_gpstart-eit->gend()-1;
                
            }
            ///// do 3'
            if (fb.gstart()>=eit->gstart()&&fb.gend()<=eit->gend()) {
              // if (_strand) _peptide_transcript_offset = _gpstart-eit->gstart()+1;
              if (_strand) _peptide_transcript_length= eit->tend+_gpend-eit->gend()-1;
              if (!_strand) _peptide_transcript_offset = eit->gend()-_gpend+1;
              // if (_strand) _peptide_transcript_length=eit->tend-(eit->gend()-_gpend+1);
            }
              //_peptide_transcript_length = _strand==1? eit->gend()-_gpend+1:_gpstart-eit->gstart()+1;
        }
//cout << "peptide transcript offset is "<< _peptide_transcript_offset << "\n";
//cout << "peptide transcript length is "<< _peptide_transcript_length << "\n";
        if (_peptide_transcript_offset==0) throw runtime_error("unable to find corresponding exon for most 5' cds - this gff model doesn't make sense!?!");
        if (_peptide_transcript_length==0) throw runtime_error("unable to find corresponding exon for most 3' cds - this gff model doesn't make sense!?!");

    }

};

struct transcript_holder {
  public:
    transcript* tptr;
    // transcript* operator*() { return tptr; }
    transcript_holder(uchar a) {
        tptr = new transcript(a);
    }
    ~transcript_holder() {
        delete tptr;
    }
    /*
    void addexon(const feature_converter& f) { exons.push_back(std::move(f)); }
    void addcds(const feature_converter& f) { cds.push_back(std::move(f)); }

    void printexons() {
        for_each(exons.begin(),exons.end(), [](feature_converter& f) { cout << "transcript: "<< f.gstart() <<" "<<f.gend() << " " << f.tstart << " " << f.tend<< "\n"; });
    }

    uint convertt2g(uint i) {
    */
};

class iter {
     // this is not the way to do it - i.e. wanna i.e. too many dereferences?!?
// have it store forward or reverse iterator?!?
};

// inline void doit(feature_converter& f,int &i, int& last_te) {
// }

template<class Compare> 
  std::vector<feature_named_min>::iterator bvfm_search(
    std::vector<feature_named_min>::iterator first, 
    std::vector<feature_named_min>::iterator last, 
    // int pos, // want signed here for additional messages?!?
    Compare comp) {
// std::vector<feature_named_min>::iterator bvfm_search(std::vector<feature_named_min>::iterator first, std::vector<feature_named_min>::iterator last, int pos, bool (*fn)(int a, int b)) {
    //r really not much point in using template machinery for everything when not actually using template machinery?!?
    // std::vector<feature_named_min>::iterator it, first = pvfm->begin();
    std::vector<feature_named_min>::iterator it;
    // std::vector<feature_named_min>::iterator it, tmp = first;
    std::iterator_traits<std::vector<feature_named_min>::iterator>::difference_type count, step;
    count = std::distance(first,last);
    // count = std::distance(first,pvfm->end());
    while (count > 0) {
        it = first;
//cout << "at : " << (first-tmp) << " checking " << it->gstart() << "\n";
        step = count / 2;
        std::advance(it, step);
        if (comp(*it)) {
        // if (it->gstart()<pos) {
            first = ++it;
            count -= step + 1;
        } else count = step;
    }
    return first;
}

class COORD_CACHE {

  public:

  std::map<const string,feat_vec> ordered_cache;
  // std::map<const string,std::pair<bool,std::vector<feature_min>>> ordered_cache;
    // static std::map<string,std::vector<feature_ext>> ordered_cache;
 //    gff_holder* gff;

    // void grab_transcripts_by_location(gff_holder& blah, string chr, uint pos) {

    COORD_CACHE() = delete;
    COORD_CACHE(gff_holder* blah) /* : gff(blah) {} ;*/ {

    /*
        cout << "here we need to make transcripts accessible by chr - put shared pointer into gff_holder?!?\n"
          " sort them all or sort if first time seeing it - i.e. just copy all into vector value in map with chr key\n"
          " then sort as required - i.e. have a pair to check if sorted?! - this is strand independent so use feature min and copy?\n";
        */

//cout << "we have " << blah->mrna_coords.size() << " transcripts to organise\n";

        int i=1;
        for (auto it = blah->mrna_coords.begin() ; it != blah->mrna_coords.end() ; it++,++i) {

            // cout << it->second.scfname() << "\n";
            // cout << it->first << "\n";
            // ordered_cache

//r use feat_named_min or pair<feat_min,name>

            //A1 check if key exists
            //A2 if not - construct vector with single value of new feature_min
            //A3 bung in to std::pair<bool,std::vector<feature_min>> bvfm with false and rvalue ref for vector - i.e. static_cast<std::vector<feature_min>>(feature_min...
            //A4 then insert temporary directly into ordered_cache.insert(std::p.... (it->second.scfname(),bvfm);
            //B else if it exists access vector and do push_back!?!
            if(ordered_cache.count(it->second.scfname())==0) {
//                cout << "["<<i<<"]"<< "first time we've seen " << it->second.scfname() << " called " << it->first << " from "<< it->second.parent() << "\n";
                // should prolly make a single statement?!?
                std::vector<feature_named_min> vtmp(1,feature_named_min(it->first,it->second.gstart(),it->second.gend()));
                // std::pair<bool,std::vector<feature_named_min>>(false,static_cast<std::vector<feature_named_min>&&>(vtmp));
                ordered_cache.insert(std::pair<std::string,feat_vec>(it->second.scfname(),feat_vec(false,static_cast<std::vector<feature_named_min>&&>(vtmp))));
            } else {
//                cout << "["<<i<<"]"<<"already seen " << it->second.scfname() << " called " << it->first << " from "<< it->second.parent() << "\n";
                feat_vec* pfv = &ordered_cache.find(it->second.scfname())->second;
                pfv->first = false;
                pfv->second.push_back(feature_named_min(it->first,it->second.gstart(),it->second.gend()));
            }
            /*
            ordered_cache.insert(
              std::pair<std::string, std::pair<bool,feature_min> > (
              it->second.scfname(),
                std::pair<bool,feature_min> (
                  false,feature_min(it->second.gstart(),it->second.gend())
                )
              )
            );
            */
        }
    }

    std::vector<string> grab_transcripts_by_location(string chr, uint pos);
    // void grab_transcripts_by_location(string chr, uint pos, std::vector<string>&);



};

/// prolly ought to use T&& return more than accepting T& just for clarity?!?
std::vector<string> COORD_CACHE::grab_transcripts_by_location(string chr, uint pos) {
// void COORD_CACHE::grab_transcripts_by_location(string chr, uint pos, std::vector<string>& list) {

// cout << "this is where we grab the coords\n";
//    auto it=gff->mrna_coords.find(chr);
    /// should use count?!?


        if(ordered_cache.count(chr)==0) throw std::runtime_error("bugger");
        // if(ordered_cache.find(chr)==ordered_cache.end()) throw std::runtime_error("bugger");

            auto it = ordered_cache.find(chr);
            if(it->second.first) {
//cout <<"already ordered\n";
            } else {
//cout << "this is an unsorted vector\n";
                std::vector<feature_named_min>* pvfm = &it->second.second;
//cout << "we have a list of " << pvfm->size() << " transcripts to sort\n";

                //// the precise sort functor is partly determined by the exact search?!?
                //cout << "here = "<<pvfm->at(0).gstart() << "\n";
                //cout << "here = "<<pvfm->front().gstart() << "\n";
                //cout << "here = "<<pvfm->back().gstart() << "\n";
                std::sort(pvfm->begin(),pvfm->end(),[](feature_named_min a,feature_named_min b){ return a.gstart()<b.gstart(); }); // forget the return without ->bool and you will have implicit cast of void->bool - i.e. downcast is ilegal - thankfully?!?
// cout << "first transcript starts at "<<pvfm->front().gstart() << "\n";
// cout << "last transcripts starts at "<<pvfm->back().gstart() << "\n";
                it->second.first=true;
            }
                // Unlike upper_bound, this function returns an iterator to the element also if it compares equivalent to value and not only if it compares greater.
                // use a dummy struct?!? auto low=std::lower_bound (pvfm->begin(), pvfm->end(), [](int a.gstart(), int p){ return  a<p });
                // auto low=std::lower_bound (pvfm->begin(), pvfm->end(),   a<p });

            std::vector<feature_named_min>* pvfm = &it->second.second; // really no need for this as sorting shouldn't invalidate iterators

///r put in lower_bound for start and upper_bound for end thus getting full range of transcripts to return with O(log(N))?!?
///r factor out the below?!? put in multiple transcripts over region - before and after - to test?!?


//// clearly want a relatively bizare condition?!?

// int i = 0;
// for_each(pvfm->begin(),pvfm->end(),[&i](feature_named_min& poo) { if(i<20) cout << "[" << i++ << "] " << poo.gstart() << "-" << poo.gend() << "\n"; });

            auto first = bvfm_search(pvfm->begin(),pvfm->end(),/* pos+1,*/[&pos](feature_named_min a)->bool{ return a.gend()<pos; });
            // auto first = bvfm_search(pvfm->begin(),pvfm->end(),pos+1,[&pos](int a)->bool{ return a<pos+1; });

//cout << "done checking at : " << (first-pvfm->begin()) << " checking " << first->gend() << "\n";
            // cout << "we start at " << first->gstart() << "\n";

            // was gonna use a function pointer?!? float result = pt2Func(a, b);    // call using 
            auto past_the_end = bvfm_search(pvfm->begin(),pvfm->end(),/*pos+1,*/[&pos](feature_named_min a)->bool{ return a.gstart()<pos+1; });
            // auto past_the_end = bvfm_search(pvfm->begin(),pvfm->end(),pos+1,[&pos](int a)->bool{ return a<pos+1; });

//cout << "done checking at : " << (past_the_end-pvfm->begin()) << " checking " << past_the_end->gstart() << "\n";
            // cout << "we start at " << last->gstart() << "\n";

//cout << "done\n";

//cout << "now we look for overlap with " <<pos<< "\n";

            std::vector<string> list;
            for_each(first,past_the_end,[&list](feature_named_min a) { list.push_back(a.id()); });
            // for(auto it = first ; it != past_the_end ; it++ ) list.push_back(it->id());

            return static_cast<std::vector<string>&&>(list);


/*    if(it==mrna_coords.end()) {
        cout << "\ndoesn't exist\n";
    } else {
        
    }
*/

}


template<typename T> inline void doit(T& it, T& eit) {
    int i=0,last_te=0;
    for ( ; it!=eit; it++,i++) {
    // for ( ; it!=eit; it++,i++) {
         cout <<"i="<<i<< " " << it->gstart() << "\n";

         // if (i==0) 
         if(i==0) {
             // already 0 initialised tstart?!?!
             it->tstart=1;
             // use ternaries e.g. it->tend = strand ? it->gend-it->gstart+1 : ...?!?
             it->tend=it->gend()-it->gstart()+1;
         } else {
             it->tstart=last_te+1;
             it->tend=it->gend()-it->gstart()+it->tstart; 

         }
        
        last_te=it->tend;
    }
}


class rebuildcache { 

  public:

    transcript* pulltranscript(std::string tid);
    rebuildcache(gff_holder& g) : gff(g) {}

    ~rebuildcache() {
//        cout << "\n\ndoing some clean up\n\n";
        for (auto it = tcache.begin() ; it != tcache.end() ; it++ ) 
          delete it->second;
    }

  private:

    std::map<string,transcript*> tcache;
    gff_holder& gff;
    rebuildcache()=delete;


};

transcript* rebuildcache::pulltranscript(std::string tid) {
// std::shared_ptr<transcript> rebuildcache(gff_holder& gff, std::string tid) {

    // should really wrap transcript in something that does automatic allocation and deallocation of transcript?!? - i.e. RAII as with string?!?
    // and just forwards everything to the actual transcript?!?

    //static std::map<string,transcript*> tcache;
    // static std::map<string,std::shared_ptr<Transcript>> tcache;
    // static std::map<string,std::unique_ptr<transcript>> tcache;
    // why the fuck would you use shared - it's just to stop you having to call delete?!? - but since it would be released at termination it's not such a biggie?!?
    
    transcript *tptr;

    if(tcache.count(tid)==0) {
//        cout << "first time access - need to build it up!?!\n";
        /// it's a bit silly to do 2x name lookups - i.e. should possibly use the actual feature_ext for lookup and then just store pointer to them...?!?
//        cout << "SOMETHING IS UP HERE = '" << gff.mrna_coords.find(tid)->second.strand() << "'\n";
        // cout << "SOMETHING IS UP HERE = '" << gff.mrna_coords.find(*it)->second.strand() << "'\n";
        // transcript trans(gff.mrna_coords.find(tid)->second.strand());


        // transcript_holder trans(gff.mrna_coords.find(tid)->second.strand());
        transcript* transptr = new transcript(gff.mrna_coords.find(tid)->second.strand());
        // std::unique_ptr<transcript> transptr(new transcript(gff.mrna_coords.find(tid)->second.strand()));

        /// grab ranges for cds and exons...

        //cout << "There are " << (int)mymm.count(c);
        // cout << " elements with key " << c << ":";
        for (auto eit=gff.exonbymrna_coords.equal_range(tid).first; eit!=gff.exonbymrna_coords.equal_range(tid).second; ++eit) 
        transptr->addexon(feature_converter(eit->second.gstart(),eit->second.gend()));
                    // cout << "exon="<<eit->second.gstart() << "-"<<eit->second.gend()<<"\n";// cout << " " << (*eit).second;
        
        for (auto eit=gff.cdsbymrna_coords.equal_range(tid).first; eit!=gff.cdsbymrna_coords.equal_range(tid).second; ++eit)
        transptr->addcds(feature_converter(eit->second.gstart(),eit->second.gend()));
            // cout << "cds="<<eit->second.gstart() << "-"<<eit->second.gend()<<"\n";// cout << " " << (*eit).second;

//cout << "putting it in!?!\n";
        tcache.insert(std::pair<std::string,transcript*>(tid,transptr));

        tptr=transptr;

    // } else cout << "ALAREADY STORED!!!!!\n";
    
    } else tptr=tcache.find(tid)->second;
//        cout << "\n\nALREADY STORED!!!!!\n\n";

    // tptr=tcache.find(tid)->second;
    return tptr;

}

// void trans2genome(gff_holder &gff);
// void trans2trans(gff_holder &gff);

inline void print_transcripts(gff_holder &gff, std::vector<string> &list_of_transcripts, unsigned int gpos, char mode) {

    for (auto it=list_of_transcripts.begin() ; it!=list_of_transcripts.end() ; it++ ) {

        rebuildcache tc(gff);
        transcript *trans = tc.pulltranscript(*it);

switch (mode) {

//r 2trans
case 't': {

        // int transpos =trans->convertg2t(gpos);
        int tpos =trans->convertg2t(gpos);
        // cout << *it << "\t" << tpos-1 << "\t" << tpos; // << "\n";
        // cout << *it << "\t" << transpos-1 << "\t" << transpos; // << "\n";
// cout << " // G2T : (second of two positions in bed) is " << gpos << " = " <<trans->convertg2t(gpos) << " ";

        // why?!? int tpos = trans->convertg2t(gpos);

        if(tpos!=0) {

            cout << *it << "\t" << tpos-1 << "\t" << tpos; // << "\n";
            

            int pos = trans->convertt2g(tpos);

            /* cout << "// this is really trans2genome reverse process - clearly this followed by genome2trans is trans2trans?!?" // cout << "AND REVERSE " */

            cout << "\t# g2t : reverse : " 
              // wtf?!? << tpos << " = " <<trans->convertt2g(tpos)<<"\n";
              << tpos << " = " << pos << "\n";

        //r clearly not relevant for trans2genome...?!?
        } else if (tpos==0) cout << " - Not exonic.\n";

        else throw runtime_error("what the fuck?!?");

    }
    break;

//r 2cds
case 'c': {

        int cdspos = trans->convertg2cds(gpos);
        cout << *it << "\t" << cdspos << "\t" << cdspos+1;
        cout << "\t# g2c : " << gpos << " = " <<trans->convertg2cds(gpos)<<"\n";
        // cout << "\t# g2c : (first of two position in bed) is " << gpos << " = " <<trans->convertg2cds(gpos)<<"\n";
    }
    break;

//r 2pep
default: { // case 'p': {

        int peppos = trans->convertg2p(gpos);
        cout << *it << "\t" << peppos << "\t" << peppos+1;
        cout << "\t# g2p: " << gpos << " = " <<trans->convertg2p(gpos)<<"\n";
        // cout << "\t# g2p: (first of two positions in bed) is " << gpos << " = " <<trans->convertg2p(gpos)<<"\n";

    }
    break;

    };


}

}

void genome2(gff_holder &gff, char* vcffile, char mode) {

    std::ifstream in(vcffile);
    // std::ifstream in("coordseqer_development/aglaope_2012Sep7.Hmel_1-1.coordseq.test.recode.vcf");
    if(in == 0) throw runtime_error("prob with vcf file");

    VCF vcf;

    COORD_CACHE coordcache(&gff); // this is pretty dirty?!?

//// should really have the genome2... part factored out for genome2 and trans2?!

    while(getvcf(in,&vcf)) { //if(vcf.chr!="HE670865") continue;        

        cout << "# "<< vcf.chr << " : " << vcf.pos << "\n";

        // cout << "\nwe need to grab transcripts for the region " << vcf.chr << " " << vcf.pos << "\n";
        //// atm., transcripts are in gff.mrna_coords - i.e. std::map<std::string,feature_ext>
        /// prolly best to generate sorted vector value with chr key?!? rather than convert to map<start> 
        // just do a sort each time if the key doesn't exist and then use binary search to pull appropriate mRNA?!?
        // grab_transcripts_by_location(gff,vcf.chr,vcf.pos);
        std::vector<string> list_of_transcripts = coordcache.grab_transcripts_by_location(vcf.chr,vcf.pos);
        //cout << "we have " << list_of_transcripts.size() << " to check\n";

        ///f                  ////////////////////
        print_transcripts(gff, list_of_transcripts, vcf.pos, mode);


    }
        /// for vcf->bed need to identify transcripts by genomic location and then do conversion
        // coordcache.grab_transcripts_by_location(vcf.chr,vcf.pos,list_of_transcripts);
        /// for transcript to genomic and transcript to transcript just need to identify transcripts directly by name so no need for previous step...
        /// now that have list of transcripts - i.e. may be multiple transcripts - we're not doing this in gene-centric fashion just do the conversion...?!?
        // if(*it!="HMEL001038-RA") continue;

            //cout << "rebuild "<<*it << " : ";

            ///r except need to reconsistute the transcripts first with a cache - i.e. need to generate transcripts as follows - just grab all exons/cds per
            ///r returned transcript - i.e. above needs to return vector of names if going frm vcf?!? - then need to add ti to cache - thus have conversion object
            ///r to only reconstitute if not having done before?!?... - by having vecotr of transcript names return becomes simple to go from one transcript location to next?!?
            ///r else with bed just go through...

            //// then options for input/output file types, with/without comments giving details?!?
            // trans = tc.pulltranscript(*it);
            //cout << "\nG2T : is " << vcf.pos << " = " <<trans.convertg2t(vcf.pos)<<"\n";;
            //int tpos = trans.convertg2t(vcf.pos);
//             trans.printexons();
    //    cout << "\n\nlet's just put all transcripts into vector and sort with lambda functor with chr comparison\n\n";
    //    cout << "then when you have a specific chr put that subsection into map of vector?!?\n\n";
    // std::map<string,std::vector<feature_ext>> chr_ordered_mrna_cache;
}










void trans2genome(gff_holder &gff,char *bedfile) {

    std::ifstream in(bedfile);
    // std::ifstream in("coordseqer_development/presumably_base_positions_of_variants_in_exon_EXON_N.bed");
    if(in == 0) throw runtime_error("prob with bed file");

    BED bed;
    while(getbed(in,&bed)) { 
        
        if (bed.end-bed.start!=1) 
          throw runtime_error("not currently implemented for ranges?!?");
        // if(auto it = .find() !=?!? vs count,find?!?
        // if((std::map<std::string,feature_ext>::iterator it = gff.mrna_coords.find(bed.trans))==gff.mrna_coords.end()) {

        if (gff.mrna_coords.count(bed.trans)==0)
          throw runtime_error("never heard of "+bed.trans);

        rebuildcache tc(gff);
        transcript *trans = tc.pulltranscript(bed.trans);

        /* should check position is within transcript bounds?!? */
        int gpos = trans->convertt2g(bed.end);
        cout << gff.mrna_coords.find(bed.trans)->second.scfname() << "\t" << gpos-1 << "\t" << gpos;
        cout << " // " << bed.end << " = " << gpos <<"\n";
    }
}

void trans2trans(gff_holder &gff,char *bedfile) {

    std::ifstream in(bedfile);
    // std::ifstream in("coordseqer_development/presumably_base_positions_of_variants_in_exon_EXON_N.bed");
    if(in == 0) throw runtime_error("prob with bed file");

    COORD_CACHE coordcache(&gff); 

    BED bed;
    while(getbed(in,&bed)) { 
        
        cout << "# Running "<< bed.trans << " : " << bed.start << "-" <<  bed.end << "\n";

        if (bed.end-bed.start!=1) 
          throw runtime_error("not currently implemented for ranges?!?");

        if (gff.mrna_coords.count(bed.trans)==0)
          throw runtime_error("never heard of "+bed.trans);

        rebuildcache tc(gff);
        transcript *trans = tc.pulltranscript(bed.trans);

        /* should check position is within transcript bounds?!? */
        int gpos = trans->convertt2g(bed.end);

//////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<string> list_of_transcripts = coordcache.grab_transcripts_by_location(
          gff.mrna_coords.find(bed.trans)->second.scfname(),
          gpos
        );

        print_transcripts(gff, list_of_transcripts, gpos, 't');

        /* for (auto it=list_of_transcripts.begin() ; it!=list_of_transcripts.end() ; it++ ) {

            rebuildcache tc(gff);
            transcript *trans = tc.pulltranscript(*it);
    
            int transpos =trans->convertg2t(gpos);
            cout << *it << "\t" << transpos-1 << "\t" << transpos; // << "\n";
            cout << " // G2T : (second of two positions in bed) is " << gpos << " = " <<trans->convertg2t(gpos) << "\n";

            // genome2trans
            int tpos = trans->convertg2t(gpos);
            if(tpos!=0) {
                int pos = trans->convertt2g(tpos);
                cout << "// this is really trans2genome reverse process - clearly this followed by genome2trans is trans2trans?!?" // cout << "AND REVERSE " 
                << tpos << " = " <<trans->convertt2g(tpos)<<"\n";
            //r clearly not relevant for trans2genome...?!?
            } else if (tpos==0) cout << "Not in an exon - not checking reverse mapping\n";
            else throw runtime_error("what the fuck?!?");

            int cdspos = trans->convertg2cds(gpos);
            cout << *it << "\t" << cdspos << "\t" << cdspos+1;
            cout << " // G2CDS : (first of two position in bed) is " << gpos << " = " <<trans->convertg2cds(gpos)<<"\n";

            int peppos = trans->convertg2p(gpos);
            cout << *it << "\t" << peppos << "\t" << peppos+1;
            cout << " // G2P : (first of two positions in bed) is " << gpos << " = " <<trans->convertg2p(gpos)<<"\n";
            cout << "#\n";

        } */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

    }

}

int main (int argc, char **argv) {

    if (argc <= 3) {
        std::cout << CORRECT_USAGE << std::endl;
        return 0;
    }

    struct stat buf;
    if (stat(*++argv,&buf)!=0) {
    // if (stat(*(argv+2),&buf)!=0) {
        cerr << "\n*** Input gff file " << (char)0x60 << *argv << "' does not exists. ***\n\n" << CORRECT_USAGE << "\n";
        // cerr << "\n*** Input gff file " << (char)0x60 << *(argv+2) << "' does not exists. ***\n\n" << CORRECT_USAGE << "\n";
        return 1;
    }

    cout << "# parsing gff\n";
    clock_t t0, ti;   
    t0 = clock();   

    gff_holder gff;

    parse_gff_for_loading(*argv,&gff);

    ti = clock();   
    float diff = (((float)ti - (float)t0) / 1000000.0F ) * 1000;   
    printf("# parsed gff : %f ms\n",diff);

    // std::string mode(*(av+1));
    // std::cout <<"Mode="<< *++argv << std::endl;
    
    if(strcmp(*++argv,"genome2trans")==0) {

        cerr << "# running genome2trans\n";
        if (stat(*++argv,&buf)!=0) {
            cerr << "\n*** Input vcf file " << (char)0x60 << *argv << "' does not exists. ***\n\n" << CORRECT_USAGE << "\n";
            return 1;
        }
        // cerr << "\n\n*** genome2 ***\n\n";
        // std::ifstream in("coordseqer_development/aglaope_2012Sep7.Hmel_1-1.coordseq.test.recode.vcf");
        genome2(gff,*argv,'t');

    } else if(strcmp(*argv,"genome2cds")==0) {
        cerr << "# running genome2cds\n";
        if (stat(*++argv,&buf)!=0) {
            cerr << "\n*** Input vcf file " << (char)0x60 << *argv << "' does not exists. ***\n\n" << CORRECT_USAGE << "\n";
            return 1;
        }
        genome2(gff,*argv,'c');
    } else if(strcmp(*argv,"genome2pep")==0) {
        cerr << "# running genome2pep\n";
        if (stat(*++argv,&buf)!=0) {
            cerr << "\n*** Input vcf file " << (char)0x60 << *argv << "' does not exists. ***\n\n" << CORRECT_USAGE << "\n";
            return 1;
        }
        genome2(gff,*argv,'p');
    } else if(strcmp(*argv,"trans2genome")==0) {
        cerr << "# running trans2genome\n";
        if (stat(*++argv,&buf)!=0) {
            cerr << "\n*** Input bed file " << (char)0x60 << *argv << "' does not exists. ***\n\n" << CORRECT_USAGE << "\n";
            return 1;
        }
        // std::ifstream in("coordseqer_development/presumably_base_positions_of_variants_in_exon_EXON_N.bed");
        trans2genome(gff,*argv);
    } else if(strcmp(*argv,"trans2trans")==0) {
        cerr << "# running trans2trans\n"; // cout << "\n\n*** trans2trans ***\n\n";
        if (stat(*++argv,&buf)!=0) {
            cerr << "\n*** Input bed file " << (char)0x60 << *argv << "' does not exists. ***\n\n" << CORRECT_USAGE << "\n";
            return 1;
        }
        trans2trans(gff,*argv);
    } else {
        cerr << "\n*** Do not recognise " << (char)0x60 << *argv << "'. ***\n\n" << CORRECT_USAGE << "\n";
        return 1;
    }

    ti = clock();   
    diff = (((float)ti - (float)t0) / 1000000.0F ) * 1000;   
    printf("# done gff : %f ms\n",diff);

    return 0;

/*if (mode=="single") {
} else if (mode=="gff3") {
} else if (mode=="list") {
List(ac-1,av+1);
return 0;*/

/// clearly this will be a mode?!?

// need to access transcripts by position -  either sorted in vector and cache or just in a map?!?

    /// put in mode genome2trans/cds/pep trans2trans, trans2genome...

    /// put in comments option?!?

    /// formats/files?!?
    
    /// option fo rofrce ignoring of gff problems?!?

/* put in the stuff the build up the transcripts!?!?!?
    transcript trans(1);
    trans.addexon(feature_converter(14500,15000));
    trans.addexon(feature_converter(12000,14000));
    trans.addexon(feature_converter(10000,11000));
    trans.addcds(feature_converter(14500,14600));
    trans.addcds(feature_converter(12000,14000));
    trans.addcds(feature_converter(10500,11000));

    int pos =1254;
    trans.printexons();
    cout << "G2T : is " << pos << " = " <<trans.convertt2g(pos)<<"\n";;
    trans.printexons();
    pos = trans.convertt2g(pos);
    cout << "AND REVERSE " << pos << " = " <<trans.convertg2t(pos)<<"\n";;
    trans.printexons();
    cout << "and intron " << pos << " = " <<trans.convertg2t(11500)<<"\n";;
    // for_each(feats.begin(),feats.end(), [](feature_converter& f) { cout << f.gstart <<" "<<f.gend << "\n"; });
        ////check it exists!?!
        /// can you range construct vector from return of 
        //r/ just static_cast<T&&> if created?!?
        // std::for_each(
*/
// cout<<"here="<<sizeof(bool)<<"\nhere="<<sizeof(signed char)<<"\n";

/*
cout <<"\n\n\n";
    transcript trans(1);
    trans.addexon(feature_converter(14500,15000));
    trans.addexon(feature_converter(12000,14000));
    trans.addexon(feature_converter(10000,11000));
    trans.addcds(feature_converter(14500,14600));
    trans.addcds(feature_converter(12000,14000));
    trans.addcds(feature_converter(10500,11000));

    int pos =1254;
//trans.printexons();
    cout << "G2T : is " << pos << " = " <<trans.convertt2g(pos)<<"\n";;
//trans.printexons();
    pos = trans.convertt2g(pos);
    cout << "AND REVERSE " << pos << " = " <<trans.convertg2t(pos)<<"\n";;
//trans.printexons();
    cout << "and intron " << pos << " = " <<trans.convertg2t(11500)<<"\n";;
    // for_each(feats.begin(),feats.end(), [](feature_converter& f) { cout << f.gstart <<" "<<f.gend << "\n"; });

    transcript trans1(0);
    trans1.addexon(feature_converter(14500,15000));
    trans1.addexon(feature_converter(12000,14000));
    trans1.addexon(feature_converter(10000,11000));
    trans1.addcds(feature_converter(14500,14600));
    trans1.addcds(feature_converter(12000,14000));
    trans1.addcds(feature_converter(10500,11000));

    pos =1254;
//trans1.printexons();
    cout << "GET TRANSCRIPT POSITION " << pos << " = " <<trans1.convertg2t(pos)<<"\n";;
    cout << "G2T : is " << pos << " = " <<trans1.convertt2g(pos)<<"\n";;
//trans1.printexons();
    pos = trans1.convertt2g(pos);
    cout << "AND REVERSE " << pos << " = " <<trans1.convertg2t(pos)<<"\n";;
    cout << "and intron " << pos << " = " <<trans1.convertg2t(11500)<<"\n";;
//trans1.printexons();

    // for_each(feats.begin(),feats.end(), [](feature_converter& f) { 

   //  if (1) while ((first!=last)&&(first!=--last)) swap (*first++,*last);
   // this really shouldn't be necessary - need to?!?
//    if(strand) reverse(feats.begin(),feats.end());
    // clearly this is strand dependent - reverse order and use of start/end

    // std::vector<feature_converters>::iterator sit, eit;
    return 0;

// don't want to reverse the elements?!?
// std::vector<feature_converter>::iterator it, eit;
// deque and vector are random access?!?

//////////////////// perhaps use type erasure to allow iterators to be used generaically?!?

*/

}


