#include "shared/reference.h"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <string>
#include <map>
#include <utility>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_set>

using namespace std;
using position = pair<string,int>; // 基因组位点的表示形式
using mop = map<unsigned int,int>; // 记录read每个kmer对应的reads上的偏离值

// 方便将基因组绝对位置转化为 chrom:pos
vector<string> v_chroms = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
                                "19","20","21","22","X","Y"};

unordered_set<string> s_chroms(v_chroms.begin(),v_chroms.end());

//数据结构，记录kmer连锁的信息
//@pos : unsigned 基因组绝对位置
//@offset: int kmer距离断裂点长度
//@chainNum: int 连锁的kmer 数量
struct Kmer_Chain
{
    unsigned pos; //基因组位置
    int offset; //接近断裂点kmer 距离read
    int chainNum;   //连锁的kmer 数量
};

//构建染色体号偏离值对应表。用于以绝对值记录染色体位置。
//@fasta: hg19 基因组文件地址
//return 染色体号与偏离值对应表
map<string,unsigned int> make_chrom_offset(const string& fasta){
    ifstream fai{fasta + ".fai"};
    string chrom;
    unsigned length,offset,base_len,byte_len;
    map<string,unsigned> map_chrom_offset;
    while (fai>>chrom>>length>>offset>>base_len>>byte_len)
    {
        map_chrom_offset[chrom] = offset;
        // cout << chrom << "\t" << offset << "\n";
    }
    return map_chrom_offset;
}

//二分法寻找染色体绝对位置对应的实际染色体号和实际位置
//@pos: hg19文件碱基绝对位置；
//@map_chrom_offset: 从hg19 索引文件中得到的，各个染色体的初始偏离值；
//@v_chroms： 按染色体号排序的vector;
//@s: 待检索的染色体vector 的起始位置；
//@e: 待检索的染色体vector 的终止位置。
position int2position(unsigned pos,map<string,unsigned>& map_chrom_offset,vector<string>& v_chroms,int s,int e){
    if(s >= e){
        if(pos > map_chrom_offset[v_chroms[s]]) {
            position real_pos = make_pair(v_chroms[s],pos-map_chrom_offset[v_chroms[s]]);
            return real_pos;
        } else {
            position real_pos = make_pair(v_chroms[s-1],pos-map_chrom_offset[v_chroms[s-1]]);
            return real_pos;
        }
    } else {
        int half = (s + e) / 2;
        if (pos > map_chrom_offset[v_chroms[half]]) s = half + 1;
        else e = half - 1;
        // cout << "Search: " << pos << "\t In " << s <<"( " << map_chrom_offset[v_chroms[s]] << " )"<<"-"<<e << "( " <<map_chrom_offset[v_chroms[e]] << " )"<<"\n";
        return int2position(pos,map_chrom_offset,v_chroms,s,e);
    }
}

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
        // cout << seq << "\n";
        kmer_offset[trans_kmer(seq)] = n ;
        n += kmer_len;
        if(n+kmer_len >= str_len){
            n = str_len - kmer_len;
            string seq = read.substr(n,kmer_len);
            // cout << seq << "\n";
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
bool distance_check(int pos,const vector<unsigned>& v_pos,int s,int e ,int dis,int x){
    if(s > e){
        return false;
    }
    int m = (s + e)/2;
    // if(abs(pos + dis - v_pos[m]) <= x ) return true;
    if((pos + dis + x >= v_pos[m]) && (v_pos[m] + x >= pos + dis)) return true;
    else if(pos < v_pos[m]) e = m - 1;
    else s = m + 1;
    return distance_check(pos,v_pos,s,e,dis,x);
}

//检索kmer在index中的位置 vector, 过滤超多位点的kmer。将剩余的kmer按位点出现次数，由少到多排序。
//@map_kmer_pos: 事先构建好的index
//@map_kmer_offset: read中kemr 对应的偏离值；
//@threshold: 超过固定数量位点的kmer 就不进行后续的检测；
//return : kmer 集合，按在基因组存在的位点数量由少到多排序
vector<unsigned int> cheap_kmer(unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,mop& map_kmer_offset,int threshold){
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
vector<unsigned int> select_kmer(unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,mop& map_kmer_offset,int threshold){
    map<unsigned int,unsigned int> kmer_pos_size;
    vector<unsigned int> v_kmer;
    for(auto i:map_kmer_offset){
        if(map_kmer_pos[i.first].size() > threshold) {
            cout << i.first << "\t" << map_kmer_pos[i.first].size() << "\tignore" <<"\n";
            continue;
        }
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
//return: pair数据，记录能找到最靠近断裂点的kmer ; 
// 其在基因组上的位置，以及在read 上的偏离值（由0开始），连锁分值, 若没能找到连锁，返回值的偏离为 -1，连锁分值0.
Kmer_Chain global_align(unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,mop& map_kmer_offset,vector<unsigned int> v_kmer,int err){
    Kmer_Chain result = {0,-1,0};
    int score = 1; //至少要求2 段kmer 可以连锁；
    if(v_kmer.size() < 3) return result;
    cout << "Kmer Chain " << "kemr-vector-len: " << v_kmer.size() << "\n";
    for(unsigned int i = 0;i < v_kmer.size()/2 + 1;i++){
        auto kmer = v_kmer[i];
        auto offset = map_kmer_offset[kmer];
        auto v_pos = map_kmer_pos[kmer];
        int e = (v_kmer.size() - i) / err;

        // int n = 1;
        for (auto pos : v_pos){
            int m = 1;
            int k = 0;
            for(unsigned int j = i + 1;j < v_kmer.size();j++){
                auto kmer_j = v_kmer[j];
                auto dis = map_kmer_offset[kmer_j] - offset;
                auto v_pos_j = map_kmer_pos[kmer_j];
                // cout << "sites: " << v_pos_j.size();
                if (distance_check(pos,v_pos_j,0,(int)(v_pos_j.size())-1,dis,3)) m ++;
                else k++;
                if (k > e) break;
            }
            
            // cout << v_pos.size() << "\t" << "Num." << n << "\t" << i <<"\t"<<m<<"\t"<<e<<"\n";
            // n++;
            if(m > score) {
                result = {pos,offset,m};
                score = m;
                cout << "Kmer vector size: " << v_kmer.size() - i
                    << " Kmer Index: "  << i <<" Kmer Offset: " << offset
                    <<" Pos: " << pos << " Score: " << m << "\n";
            }
            if ( m >= (v_kmer.size() - i - 1) && m > 1) return result;
        }        
    }
    return result;
} 


//反向互补序列
string reverse_complete(const string& inseq){
    string s;
    for(char c : inseq){
        switch (c)
        {
        case 'A':
            s += 'T';
            break;
        case 'T':
            s += 'A';
            break;
        case 'G':
            s += 'C';
            break;
        case 'C':
            s += 'G';
            break;
        default:
            break;
        }
    }
    reverse(s.begin(),s.end());
    return s;
}


//比对到基因组
void mapping(const string& seq,unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,map<string,unsigned> map_chrom_offset){
    mop kmer_f = extract_kmer(seq,12);
    vector<unsigned int> v_kmer = select_kmer(map_kmer_pos,kmer_f,3000);
    Kmer_Chain kcf = global_align(map_kmer_pos,kmer_f,v_kmer,3);
    if (kcf.pos != 0){
        position pf = int2position(kcf.pos,map_chrom_offset,v_chroms,0,23);
        cout << pf.first << ":" << pf.second << "\n";
    } else 
        cout << "cannot mapping!\n";    
}


int main(int argc,char* argv[]){
    if(argc < 4){
        cout << argv[0] << "\tsequence.txt \treference\t12-kmer-index\n";
        return -1;
    }
    
    string seqFile = argv[1];
    string fasta = argv[2];
    map<string,unsigned> map_chrom_offset = make_chrom_offset(fasta);
    fastqReader reference = fastqReader(fasta);
    ifstream kmer_index {argv[3]};
    unordered_map<unsigned int,vector<unsigned>> map_kmer_pos;
    boost::archive::binary_iarchive iarch(kmer_index);
    iarch >> map_kmer_pos;
    kmer_index.close();

    ifstream sf{seqFile};
    string seq;
    while(sf >> seq){
        cout << "Mapping: " << seq << "\n" << "Forward: \n";
        mapping(seq,map_kmer_pos,map_chrom_offset);
        cout << "RP: \n";
        string seq_rp = reverse_complete(seq);
        mapping(seq_rp,map_kmer_pos,map_chrom_offset);
        cout << "\n++++++++++\n";
    }
}


