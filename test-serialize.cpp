#include <unordered_map>
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
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
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
            return -1;
        }
    }
    return kmer;
}

// 记录序列在基因组中的位置
// struct Position
// {
//     uint8_t chrom;
//     unsigned int pos;
//     // boost serialize 
// private: 
//     friend class boost::serialization::access; 
//     template <typename Archive> 
//     void serialize(Archive &ar, const unsigned int version) { 
//         ar & chrom; 
//         ar & pos; 
//     } 
// };

string random_seq(){
    default_random_engine generator;
    uniform_int_distribution<int> distribution(0,99);
    vector<char> atcg = {'A','T','C','G'};
    string seq;
    for(int i =0 ; i < 12 ;i++){
        seq += atcg[distribution(generator)%4];
    }
    return seq;
}

int main()
{
//    std::unordered_map<int, int> map = {{1,2}, {2,1}};
//    std::stringstream ss;
//    boost::archive::text_oarchive oarch(ss);
//    oarch << map;
//    std::unordered_map<int, int> new_map;
//    boost::archive::text_iarchive iarch(ss);
//    iarch >> new_map;
//    std::cout << (map == new_map) << std::endl;   
    pair<uint8_t,unsigned int> p1 = make_pair(1,111111);
    pair<uint8_t,unsigned int> p2 = make_pair(1,111112);
    pair<uint8_t,unsigned int> p3 = make_pair(1,111113);
    using vp = vector<pair<uint8_t,unsigned int>> ;
    vp v1 = {p1};
    vp v2 = {p1,p2};
    vp v3 = {p1,p2,p3};
    unordered_map<unsigned int,vp> map ;
    map[trans_kmer("GAGTTTTTTTTT")] = v1;
    map[trans_kmer("GAGTTTTTTTTA")] = v2;
    map[trans_kmer("GAGTTTTTTTTG")] = v3;

    default_random_engine generator;
    uniform_int_distribution<int> distribution1(1,200);
    uniform_int_distribution<int> distribution2(1,9999999);
    for (int i = 0;i<100;i++){
        string seq = random_seq();
        unsigned int key = trans_kmer(seq);
        int n = distribution1(generator);
        for (int j = 0;j<n;j++){
            pair<uint8_t,unsigned int> p = make_pair(1,distribution2(generator));
            map[key].push_back(p);
        }
    }
    //stringstream ss;
   { ofstream ss {"b.map"};
// 序列化
    boost::archive::binary_oarchive oarch(ss);
    oarch << map;
   }
// 解序列化
{   ifstream ss {"b.map"};
    unordered_map<unsigned int,vp> new_map;
    boost::archive::binary_iarchive iarch(ss);
    iarch >> new_map;
    cout << (map == new_map) << std::endl; 
}
    ifstream infile {"a.fai"};
    string chrom;
    int length,offset,lineBase,lineByte;
    while (infile >> chrom >> length >> offset >> lineBase >> lineByte)
    {
        cout << "==> " << chrom << ";" << length << ";" << offset << "\n";
    }
    

}