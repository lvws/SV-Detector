#include "shared/BamReader.h"
#include <string>
#include <htslib/sam.h>
#include <iostream>
#include <sstream>

using namespace std;

//Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
char fourbits2base(uint8_t val) {
    switch(val) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            cerr << "ERROR: Wrong base with value "<< (int)val << endl ;
            return 'N';
    }
}

bool BamReader::next()
{
    int tag = sam_read1(sam_file,bam_header,aln);
    if (tag >= 0) return true;
    return false;
}

bool Alignment::next()
{
    int tag = sam_itr_next(sam_file, iter, aln) ;
    if (tag >= 0) return true;
    return false;
}


string getCigar(const bam1_t *aln)
{
    uint32_t *data = (uint32_t *)bam_get_cigar(aln);
    int cigarNum = aln->core.n_cigar;
    stringstream ss;
    for(int i=0; i<cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_opchr(val);
        uint32_t len = bam_cigar_oplen(val);
        ss << len << op;
    }
    return ss.str();
}

// get seq quality
string getQual(const bam1_t *aln) {
    uint8_t *data = bam_get_qual(aln);
    int len = aln->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        s[i] = (char)(data[i] + 33);
    }
    return s;
}

// get seq 序列的信息记录在8bit 的数据结构中，前4bit 是前面的碱基，后4bit 是后面的碱基
string getSeq(const bam1_t *aln) {
    uint8_t *data = bam_get_seq(aln);
    int len = aln->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        char base;
        if(i%2 == 1)
            base = fourbits2base(data[i/2] & 0xF); 
        else
            base = fourbits2base((data[i/2]>>4) & 0xF);
        s[i] = base;
    }
    return s;
}



// 得到辅助字段的信息 tag = "MD" 
string getAux(const bam1_t *aln,const char tag[2])
{
    kstring_t res = KS_INITIALIZE;    // 初始化很重要
    if(bam_aux_get_str(aln,tag,&res) == 1) 
    {
        int len = ks_len(&res);
        char *ks_s = ks_str(&res);
        string s(len, '\0');
        for (int i = 0;i<len;i++ ){
            s[i] = ks_s[i];
        }
        ks_free(&res); // 资源释放很重要
        return s;
    }    
    else 
    {
        cerr << "no tag :" << tag << '\n';
        ks_free(&res);
        return "";
    }
    
}

