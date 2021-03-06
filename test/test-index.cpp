#include <unordered_map>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <random>
#include <math.h>
#include <algorithm>
#include <iterator>
#include "shared/reference.h"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace std;
using position = pair<uint8_t,unsigned int>; // 基因组位点的表示形式
using vp = vector<position> ; // 基因组位点集合
using mop = map<unsigned int,int>; // 记录read每个kmer对应的reads上的偏离值

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

//切分reads 以kmer 长度，最后的碱基不够可以和倒数第二个序列overlap.
//@read : 待查询的序列
//@kmer_len: 种子长度，应该设定为12
//return: kmer(序列转换成数字)和偏离值的映射表
mop extract_kmer(const string& read,int kmer_len){
    mop kmer_offset;
    string::size_type n = 0;
    auto str_len = read.size();
    while (true)
    {
        string seq = read.substr(n,kmer_len);
        cout << seq << "\n";
        kmer_offset[trans_kmer(seq)] = n ;
        n += kmer_len;
        if(n+kmer_len > str_len){
            n = str_len - kmer_len;
            string seq = read.substr(n,kmer_len);
            cout << seq << "\n";
            kmer_offset[trans_kmer(seq)] = n ;
            break;
        }
    }
    return kmer_offset;
}

//二分法，检查给定的 位点<chrom,pos>， 与一个 位点集合（点位是按顺序存储的） 是否间隔给定的距离dis，允许误差x.
//检索到返回 true ;否则 false.
//@pos： 待检测的位点
//@v_pos ：目标位点集合
//@s: 目标位点集合开始检索位置
//@e: 目标位点集合终止检索位置
//@dis: 待检测位点pos 和目标位点集合中目标位点的距离
//@x: 待检测位点pos 和目标位点集合中目标位点的距离,可以允许的偏差（中间kmer 可能有indel）
//return : true 找到了连锁位点； false 没有找到连锁位点。
bool distance_check(position pos,const vp& v_pos,int s,int e ,int dis,int x){
    if(s > e){
        return false;
    }
    int m = (s + e)/2;
    if(pos.first < v_pos[m].first) e = m - 1;
    else if(pos.first > v_pos[m].first) s = m + 1;
    else if(abs((int)pos.second + dis - (int)v_pos[m].second) <= x ) return true;
    else if(pos.second < v_pos[m].second) e = m - 1;
    else s = m + 1;
    return distance_check(pos,v_pos,s,e,dis,x);
}


// 将kmer 连锁在一起，得到一个 vector 记录可能起始位置，vector 按支持kmer 数排序
// 一半及以上kmer 可以链接才记录成 hotspot。
//@map_kmer_pos: 事先构建好的index；
//@map_kmer_offset: read中kemr 对应的偏离值；
//@v_kmer: 经过，过滤，排序后的kmer集合
//return: 分值按高到低的可能序列比对到的基因组位置
vp make_chain(unordered_map<unsigned int,vp>& map_kmer_pos,mop& map_kmer_offset,vector<unsigned int> v_kmer){
    map<position,int> map_hotspot;
    vp v_hotspot;
    int half = max((int)v_kmer.size()/2,2) ;
    for(vector<uint8_t>::size_type i = 0 ;i < half;i++){
        for(auto p : map_kmer_pos[v_kmer[i]]){
            int offset = map_kmer_offset[v_kmer[i]];
            position newpos = p;
            newpos.second -= offset;
            int n = 1;
            for(auto j = i + 1;j < v_kmer.size();j++){
                if(n + v_kmer.size() - j < half) break;
                int dis = map_kmer_offset[v_kmer[j]] - offset;
                auto vp_j = map_kmer_pos[v_kmer[j]];
                if(distance_check(p,vp_j,0,vp_j.size()-1,dis,2)) n++;
            }
            if ( n >= half && map_hotspot.find(newpos) == map_hotspot.end())
                map_hotspot[newpos] = n;
        }
    }
    for (auto i : map_hotspot) v_hotspot.push_back(i.first);
    sort(v_hotspot.begin(),v_hotspot.end(),[&](position a,position b) {return map_hotspot[a] > map_hotspot[b];});
    cout << "HotSpot:\n";
    for(auto i:v_hotspot) cout <<  "Pos: " << i.second << "\t" << map_hotspot[i] << "\n";
    return v_hotspot;
}


//检索kmer在index中的位置 vector, 过滤超多位点的kmer。将剩余的kmer按位点出现次数，由少到多排序。
//@map_kmer_pos: 事先构建好的index
//@map_kmer_offset: read中kemr 对应的偏离值；
//@threshold: 超过固定数量位点的kmer 就不进行后续的检测；
//return : kmer 集合，按在基因组存在的位点数量由少到多排序
vector<unsigned int> cheap_kmer(unordered_map<unsigned int,vp>& map_kmer_pos,mop& map_kmer_offset,int threshold){
    map<unsigned int,unsigned int> kmer_pos_size;
    vector<unsigned int> v_kmer;
    for(auto i:map_kmer_offset){
        kmer_pos_size[i.first] = map_kmer_pos[i.first].size();
        v_kmer.push_back(i.first);
    }
    sort(v_kmer.begin(),v_kmer.end(),[&](unsigned a,unsigned b){return kmer_pos_size[a] < kmer_pos_size[b];});
    return v_kmer;
}

//检索kmer在index中的位置 vector, 过滤超多位点的kmer。将剩余的kmer按顺序由左至右排序。
//@map_kmer_pos: 事先构建好的index
//@map_kmer_offset: read中kemr 对应的偏离值；
//@threshold: 超过固定数量位点的kmer 就不进行后续的检测；
//return : kmer 集合，按read的 offset 排序，由低到高。
vector<unsigned int> select_kmer(unordered_map<unsigned int,vp>& map_kmer_pos,mop& map_kmer_offset,int threshold){
    map<unsigned int,unsigned int> kmer_pos_size;
    vector<unsigned int> v_kmer;
    for(auto i:map_kmer_offset){
        kmer_pos_size[i.first] = map_kmer_pos[i.first].size();
        v_kmer.push_back(i.first);
    }
    sort(v_kmer.begin(),v_kmer.end(),[&](unsigned a,unsigned b){return map_kmer_offset[a] < map_kmer_offset[b];});
    return v_kmer;
}


//从断点位置开始检索kmer 在基因组中的位置；
//对于n 长度的kmer，断点开始第i 项（ i in [0,n] )；连锁kmer 数量 m，当 m >= n - (i+1) 时终止后续的检索。
//设错误阈值值为e , (可以设 e = n / 3), 非连锁值k, 当 k > e 时，放弃该位点的检索。
//至少要求2 段kmer 可以连锁；
//@map_kmer_pos: 事先构建好的index；
//@map_kmer_offset: read中kemr 对应的偏离值；
//@v_kmer: 经过，过滤，排序后的kmer集合(按离开断裂点的距离增大排序)
//@err: 平均多少个kmer,允许一个错误。
//return: pair数据，记录能找到最靠近断裂点的kmer ; 其在基因组上的位置，以及在read 上的偏离值（由0开始）, 若没能找到连锁，返回值的偏离为 -1.
pair<position,int> global_align(unordered_map<unsigned int,vp>& map_kmer_pos,mop& map_kmer_offset,vector<unsigned int> v_kmer,int err){
    pair<position,int> result = make_pair({0,0},-1);
    int score = 1; //至少要求2 段kmer 可以连锁；
    for(unsigned int i = 0;i < v_kmer.size() - 1;i++){
        auto kmer = v_kmer[i];
        auto offset = map_kmer_offset[kmer];
        auto v_pos = map_kmer_pos[kmer];
        int e = (v_kmer.size() - i) / err;
        for (auto pos : v_pos){
            int m = 1;
            int k = 0;
            for(unsigned int j = i + 1;j < v_kmer.size();j++){
                auto kmer_j = v_kmer[j];
                auto dis = map_kmer_offset[kmer_j] - offset;
                auto v_pos_j = map_kmer_pos[kmer_j];
                if (distance_check(pos,v_pos_j,0,v_pos_j.size()-1,dis,3)) m ++;
                else k++;
                if (k > e) break;
            }
            if(m > score) result = make_pair(pos,offset);
            if ( m >= (v_kmer.size() - i - 1)) return result;
        }        
    }
    return result;
} 

//利用汉明距离算法确定基因的断裂位点；
//若融合基因在下游，序列向左移动匹配找断裂点；(提供反向序列)
//若融合基因在上游，序列向右移动匹配寻找断裂点；（提供正常序列）
//滑动配可以允许最多e个错配点；
//@refer_seq: 基因组截取的断裂位点到融合基因未能完全匹配的序列；
//@query_seq: read 中截取的靠近断裂位点未能和基因组完全匹配的序列；
//@e: 允许的最大错配数量（不允许存在indel）;
//return: 待检测的两段序列可以匹配上的碱基数量；
int break_point_search(const string& refer_seq,const string& query_seq,int e){
    int n = 0; //记录匹配碱基的数量（未达到错配阈值前，错配碱基也会计算进）
    int k = 0; //记录错配数
    for(string::size_type i = 0 ; i < refer_seq.size();i++){
        if(refer_seq[i] == query_seq[i]) n = i + 1;
        else if( k > e) return n;
        else k++;
    }
    return n;
}

//对soft-clip seq 进行重新比对，找到其原本的基因组位置，并找到断裂位点；
//根据 soft-clip seq 在基因组的上下游，以及是否正链匹配（负链匹配），调整不同的检测策略；
//1. soft-clip seq 在基因组的下游、正链匹配： kmer 按离 read offset 距离由小到大排序，检查连锁；确定好基因组位点P 和offset后，进行断点寻找。
//断点寻找，从右至左（使用反向的refer_seq 和 反向的query_seq 寻找匹配碱基数n），最终位点为 P-n;
//2. soft_clip seq 在基因组上游、正链匹配：kmer 按离 read offset 距离由大到小排序，检查连锁；确定好基因组位点P 和offset 后，
//修正位置P+kmer长度、offset + kmer长度，断点寻找，从左至右匹配碱基数为n，最终位点是 P+n;
//3. soft-clip seq 在基因组的下游、负链匹配：取序列的反向互补，后续的计算方针按步骤2进行；
//4. soft_clip seq 在基因组上游、负链匹配：取序列的反向互补，后续的计算方针按步骤1进行；
//@map_kmer_pos: 事先构建好的index；
//@reference: 用于获取基因组序列；
//@soft_clip_seq: 待检索的序列；
//@chrom_map: 染色体号对应表
//@isdown: soft-clip seq 在下游；
//return : 断裂位点，pair的第二位数，-1: 没有找到位点；0: 下游正链；1：下游负链；2:上游正链；3:上游负链
pair<position,int> map_soft_clip_seq(unordered_map<unsigned int,vp>& map_kmer_pos,const fastqReader& reference,const string& soft_clip_seq,map<uint8_t,string>& chrom_map,bool isdown)
{
    pair<position,int> result = make_pair({0,0},-1);    
    mop map_kmer_offset = extract_kmer(soft_clip_seq,12);
    vector<unsigned int> v_kmer = select_kmer(map_kmer_pos,map_kmer_offset,300);
    int error_golbal = 3;
    int error_hanming = 3;
    if(isdown){
        pair<position,int> pos_global = global_align(map_kmer_pos,map_kmer_offset,v_kmer,error_golbal) ;
        if (pos_global.second == 0) {
            result.first = pos_global.first;
            result.second = 1;
            return result;
        }
        string refer_seq = reference.fetch(chrom_map[pos_global.first.first],pos_global.first.second-pos_global.second,pos_global.first.second-1);
        int n = break_point_search()
    }

    
}


int main(int argc,char* argv[])
{
    if (argc < 4){
        cout << argv[0] << "\tindexFile\tfasta\tqureySeq";
        return -1;
    }
    ifstream index{argv[1]};
    fastqReader reference = fastqReader{argv[2]};
    string read = argv[3];
    map<uint8_t,string> chrom_map = {{1,"1"},{2,"2"},{3,"3"},{4,"4"},{5,"5"},{6,"6"},{7,"7"},{8,"8"},
    {9,"9"},{10,"10"},{11,"11"},{12,"12"},{13,"13"},{14,"14"},{15,"15"},{16,"16"},{17,"17"},{18,"18"},
    {19,"19"},{20,"20"},{21,"21"},{22,"22"},{23,"X"},{24,"Y"}};

    unordered_map<unsigned int,vp> map_kmer_pos;
    boost::archive::binary_iarchive iarch(index);
    iarch >> map_kmer_pos;
    
    mop map_kmer_offset = extract_kmer(read,12);
    vector<unsigned int> v_kmer = cheap_kmer(map_kmer_pos,map_kmer_offset,300);
    vp v_hotspot = make_chain(map_kmer_pos,map_kmer_offset,v_kmer);

    cout << "Query read:\n" << read << "\n";
    int n = 0;
    for(auto i : v_hotspot){
        cout <<chrom_map[i.first] <<": " << i.second <<"-" << i.second+read.size()-1 << "\t"
        << reference.fetch(chrom_map[i.first],i.second,i.second+read.size()-1).seq() << "\n";
    }
    // string read1 = "CTGCTGTGTCCCGTGGGGCTCATGGGCCATTTGAGTCCCCTTAGCTGGTGTCTCCCTG";
    // string read2 = "TAGCTCTGGAGGTGTGGGGAGGCTCCCAAGCGGCTTCCTCAGTCACAAATACCTTCAG";

    // mop m1 = extract_kmer(read1,12);
    // for(auto i:m1) cout << i.first << ": " << i.second <<"\n";
    // mop m2 = extract_kmer(read2,12);
    // for(auto i:m2) cout << i.first << ": " << i.second <<"\n";

}