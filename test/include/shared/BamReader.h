#include <string>
#include <htslib/sam.h>
#include <iostream>
#include <sstream>

using namespace std;

class BamReader
{
public:
    BamReader(const string& name)
    :bam_name{name}
    {
        sam_file = sam_open(name.c_str(),"r"); // open bam file
        bam_index = sam_index_load(sam_file,name.c_str()); // load index
        bam_header = sam_hdr_read(sam_file); // read header
        aln = bam_init1(); //initialize an alignment
    }
    ~BamReader()
    {
    //    if(iter) sam_itr_destroy(iter);
        if(aln) bam_destroy1(aln);
        bam_hdr_destroy(bam_header);
        sam_close(sam_file);
    }
//    void setRange(const string& range); 
    bool next();


    // 需要被Alignment 访问
    hts_idx_t *bam_index;
    bam_hdr_t *bam_header;
    samFile *sam_file;
    bam1_t *aln = nullptr;
private:
    string bam_name;        
    bool fetch = false;
};

// 如果要多线程，则 iter 必须实现单独的类，不然资源无法合理分配，销毁
class Alignment
{
public:
    Alignment(samFile *sam_file_in,hts_idx_t *bam_index,bam_hdr_t *bam_header,string& ragne_in)
    :range{ragne_in},sam_file{sam_file_in}
    {
        iter = sam_itr_querys(bam_index,bam_header,range.c_str());
        aln = bam_init1(); //initialize an alignment
    }
    ~Alignment()
    {
        if(iter) sam_itr_destroy(iter);
        if(aln) bam_destroy1(aln);
    }
    bool next();
    void close()
    {
        if(iter) sam_itr_destroy(iter);
        if(aln) bam_destroy1(aln);
    }
    bam1_t *aln = nullptr;
private:
    string range; // chr10:112-221
    samFile *sam_file;
    hts_itr_t *iter = nullptr;
};

char fourbits2base(uint8_t val);
string getCigar(const bam1_t *aln);
string getQual(const bam1_t *aln);
string getSeq(const bam1_t *aln);
string getAux(const bam1_t *aln,const char tag[2]);
string getName(const bam1_t *aln);

int getFlag(const bam1_t *aln);

int getPos(const bam1_t *aln);

int getNpos(const bam1_t *aln);

int getMapq(const bam1_t *aln);

int getIsize(const bam1_t *aln);

string getNchrom(const bam1_t* aln,const bam_hdr_t* bam_header);

string getChrom(const bam1_t* aln,const bam_hdr_t* bam_header);