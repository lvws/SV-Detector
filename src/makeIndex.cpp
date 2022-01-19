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
using vp = vector<pair<uint8_t,unsigned int>> ;

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
            return -1;
        }
    }
    return kmer;
}

// 构建索引
unordered_map<unsigned int,vp> make_index(const string& fasta,int seed_len)
{
    unordered_map<unsigned int,vp> map_kmer;
    map<string,uint8_t> chrom_map = {{"1",1},{"2",2},{"3",3},{"4",4},{"5",5},{"6",6},{"7",7},{"8",8},
    {"9",9},{"10",10},{"11",11},{"12",12},{"13",13},{"14",14},{"15",15},{"16",16},{"17",17},{"18",18},
    {"19",19},{"20",20},{"21",21},{"22",22},{"X",23},{"Y",24}};
    fastqReader reference = fastqReader(fasta);
    ifstream fai {fasta+".fai"};
    string chrom;
    int length,offset,lineBase,lineByte;
    while (fai >> chrom >> length >> offset >> lineBase >> lineByte)
    {
        if(chrom_map.find(chrom) != chrom_map.end()){
            uint8_t chrom_n = chrom_map[chrom];
            int n = 1;
            while(n+seed_len-1 <= length){
                string kmer_seq = reference.fetch(chrom,n,n+seed_len-1).seq();
                unsigned int kmer = trans_kmer(kmer_seq);
                pair<uint8_t,unsigned int> pos = make_pair(chrom_n,n);
                n++;
                if (kmer == -1) continue;                
                map_kmer[kmer].push_back(pos);
                 
            }
        }
    }
    return map_kmer;
}

// 将索引序列化
void serialize(const unordered_map<unsigned int,vp>& map_kmer,const string& outfile){
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
    unordered_map<unsigned int,vp> map_kmer = make_index(argv[1],seed_len);
    serialize(map_kmer,argv[3]);
    cout << "FINISHED!\n";
}