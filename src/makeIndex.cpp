#include "shared/reference.h"
#include <unordered_map>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <random>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace std;

// 将碱基记录为两个位的数字，seed 大小为12则 12*2 = 24位，对应16进制值0x00FFFFFF 。
unsigned int trans_kmer(const string& seq){
    unsigned int kmer = 0;
    for(char c:seq){
        kmer = kmer << 2;
        switch (c)
        {
        case 'A':
            kmer += 0;
            break;
        case 'T':
            kmer += 1;
            break;
        case 'C':
            kmer += 2;
            break;
        case 'G':
            kmer += 3;
            break;
        default:
            return 1<<25;
        }
    }
    return kmer;
}

//构建染色体号偏离值对应表。用于以绝对值记录染色体位置。
//@fasta: hg19 基因组文件地址
//return 染色体号与偏离值对应表
map<string,unsigned> make_chrom_offset(const string& fasta){
    ifstream fai{fasta + ".fai"};
    string chrom;
    unsigned length,offset,base_len,byte_len;
    map<string,unsigned> map_chrom_offset;
    while (fai>>chrom>>length>>offset>>base_len>>byte_len)
    {
        map_chrom_offset[chrom] = offset;
    }
    map_chrom_offset["GL000228.1"] = map_chrom_offset["MT"];
    return map_chrom_offset;
}


// 构建索引
unordered_map<unsigned int,vector<unsigned>> make_index(const string& fasta,int seed_len)
{
    unordered_map<unsigned int,vector<unsigned>> map_kmer;
    map<string,uint8_t> chrom_map = {{"1",1},{"2",2},{"3",3},{"4",4},{"5",5},{"6",6},{"7",7},{"8",8},
    {"9",9},{"10",10},{"11",11},{"12",12},{"13",13},{"14",14},{"15",15},{"16",16},{"17",17},{"18",18},
    {"19",19},{"20",20},{"21",21},{"22",22},{"X",23},{"Y",24},{"GL000228.1",25}};
    // map<string,uint8_t> chrom_map = {{"Y",24}};
    fastqReader reference = fastqReader(fasta);
    map<string,unsigned> map_chrom_offset = make_chrom_offset(fasta);
    ifstream fai {fasta+".fai"};
    string chrom;
    unsigned length,offset,lineBase,lineByte;
    unsigned np = 1 << 25;
    while (fai >> chrom >> length >> offset >> lineBase >> lineByte)
    {
        //cout << chrom << '\n';
        if(chrom_map.find(chrom) != chrom_map.end()){
            // uint8_t chrom_n = chrom_map[chrom];
            unsigned n = 1;
            while(n+seed_len-1 <= length){
                string kmer_seq = reference.fetch(chrom,n,n+seed_len-1).seq();
                unsigned int kmer = trans_kmer(kmer_seq);
                unsigned pos = map_chrom_offset[chrom] + n;
                n++;
                if (kmer == np) continue;                
                map_kmer[kmer].push_back(pos);
                 
            }
        }
    }
    return map_kmer;
}

// 将索引序列化
void serialize(const unordered_map<unsigned int,vector<unsigned>>& map_kmer,const string& outfile){
    ofstream of{outfile};
    boost::archive::binary_oarchive oarch(of);
    oarch << map_kmer;
}

int main(int argc,char* argv[])
{
    if(argc < 4){
        cout << argv[0] << "\treference\tseedLen\tindexFileName\n";
        return -1;
    }
    stringstream ss{argv[2]};
    int seed_len;
    ss >> seed_len;
    unordered_map<unsigned int,vector<unsigned>> map_kmer = make_index(argv[1],seed_len);
    cout << "Index Make Complete!\n";
    serialize(map_kmer,argv[3]);
    cout << "FINISHED!\n";
}
