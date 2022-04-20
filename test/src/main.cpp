#include "shared/reference.h"
#include "shared/fusion.h"
#include "shared/BamReader.h"
// #include "shared/mapping.h"
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
#include <cctype>
#include <algorithm>

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

//数据结构，记录外显子的信息
//@name: 基因名（转录本号+外显子号， SGIP1:NM_001308203）
//@v_pos： 记录每个外显子的 开始位点 结束位点
//@v_extrons:   记录外显子编号顺序；
struct Extron
{
    string name;
    vector<pair<int,int>> v_pos;
    vector<int> v_extrons;
};

//读外显子bed 文件，记录外显子信息
//mge(map_gene_extron)记录基因及其外显子信息； 
//return m_extron 每个染色体对应的所有的转录本名称
map<string,vector<string>> getExtronInfo(const string& bedFile,unordered_map<string,Extron>& mge){
    map<string,vector<string>> m_extron;
    string chrom,gene,name;
    char e;
    int start,end;
    int enm;
    ifstream bed_f{bedFile};
    while(bed_f >> chrom >> start >> end >> gene){
        if(s_chroms.find(chrom) == s_chroms.end()) break;
        vector<string> vs = split(gene,":");
        name = vs[0] + ":" + vs[1];
        m_extron[chrom].push_back(name);
        stringstream ss{vs[2]};
        ss >> e >> enm;
        mge[name].name = name;
        mge[name].v_pos.push_back(make_pair(start+1,end));
        mge[name].v_extrons.push_back(enm);
    }
    return m_extron;
}

//找到基因的阅读框方向
//正向返回true,反之false
//@gene :  基因名
//@mge: 记录基因对应外显子的顺序
bool orf_dr(string& gene,unordered_map<string,Extron>& mge){
    if(gene == "NA") return true;
    vector<string> vs = split(gene,":");
    auto ves = mge[vs[0]+":"+vs[1]];
    if (ves.v_extrons.size() < 2) return true;
    if (ves.v_extrons[0] < ves.v_extrons[1]) return true;
    return false;
}


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
        // cout << "Search: " << pos << "\t In " << s <<"( " << map_chrom_offset[v_chroms[s]] << " )"<<"-"<<e << "( " <<map_chrom_offset[v_chroms[e]] << " )"<< "\t"
        //     << map_chrom_offset[v_chroms[half]] << "\n";
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
        if(map_kmer_pos[i.first].size() > threshold) continue;
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
//针对RNA序列比对进行一些更改； RNA 序列可能存在跨越内含子的问题，所以有时split read 也会出现两端mapping 到不同的外显子上；
//这时可以修正mapping 的细节，当reads 的前几段kmer序列可以连锁， reads 的后几段kmer 也能连锁； 且mapping 距离差别不大时，可以认为是跨越了外显子
//保留最开头的mapping 位置为初始mapping 位点，分值是两个分值之和。
//@map_kmer_pos: 事先构建好的index；
//@map_kmer_offset: read中kemr 对应的偏离值；
//@v_kmer: 经过，过滤，排序后的kmer集合(按离开断裂点的距离增大排序)
//@err: 平均多少个kmer,允许一个错误。
//return: pair数据，记录能找到最靠近断裂点的kmer ; 
// 其在基因组上的位置，以及在read 上的偏离值（由0开始），连锁分值, 若没能找到连锁，返回值的偏离为 -1，连锁分值0.
Kmer_Chain global_align(unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,mop& map_kmer_offset,vector<unsigned int> v_kmer,int err){
    Kmer_Chain result = {0,-1,0};
    int score = 1; //至少要求2 段kmer 可以连锁；
    cout << "Kmer Chain " << "kemr-vector-len: " << v_kmer.size() << "\n";
    if(v_kmer.size() < 2) return result;
    // cout << "Kmer Chain " << "kemr-vector-len: " << v_kmer.size() << "\n";
    for(unsigned int i = 0;i < v_kmer.size()/2 + 1;i++){
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
                // cout << "sites: " << v_pos_j.size();
                if (distance_check(pos,v_pos_j,0,(int)(v_pos_j.size())-1,dis,3)) m ++;
                else k++;
                if (k > e) break;
            }
            if(m > score) {
                result = {pos,offset,m};
                score = m;
                cout << "Kmer vector size: " << v_kmer.size() - i
                    << " Kmer Index: "  << i <<" Kmer Offset: " << offset
                    <<" Pos: " << pos << " Score: " << m << "\n";
            }
            if ( m > (v_kmer.size() - i - 1) && m > 1) return result;
        }        
    }
    return result;

} 

//利用汉明距离算法确定基因的断裂位点；
//若融合基因在下游，序列向左移动匹配找断裂点；(提供反向序列)
//若融合基因在上游，序列向右移动匹配寻找断裂点；（提供正常序列）
//滑动配可以允许最多e个错配点；
//当出现第一个错误时，若后续三个错误在较短的范围内出现，则在第一次出现错误处终止比对。
//@refer_seq: 基因组截取的断裂位点到融合基因未能完全匹配的序列；
//@query_seq: read 中截取的靠近断裂位点未能和基因组完全匹配的序列；
//@e: 允许的最大错配数量（不允许存在indel）;
//return: 待检测的两段序列可以匹配上的碱基数量；
int break_point_search(const string& refer_seq,const string& query_seq,int e){
    int n = 0; //记录匹配碱基的数量（未达到错配阈值前，错配碱基也会计算进）
    int m = 0; //记录完全匹配的碱基数。
    int k = 0; //记录错配数
    for(string::size_type i = 0 ; i < refer_seq.size();i++){
        // cout << refer_seq[i] << "-" << query_seq[i] << "\n";
        if(refer_seq[i] == query_seq[i]) n = i + 1;
        else {
            if(k==0) m = n;
            k++;
            if(k>=e){
                if(n-m < 7) return m;
                else return n;
            }
        }
    }
    return n;
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

//比较两个位点的大小，先比染色体，再比较位点大小
//map<string,int> map_c2i 是对string chrom 转换为 int 用于比较大小；
//break_p1 是 chrom:pos <string>;
//chrom_pos 是 position;
//输入1 大于 输入2  返回 ture else false
bool comparePos(string& break_p1,position chrom_pos){
    vector<string> v_b1 = split(break_p1,":");
    if (map_c2i[v_b1[0]] > map_c2i[chrom_pos.first] || (map_c2i[v_b1[0]] == map_c2i[chrom_pos.first] && string2number<int>(v_b1[1]) > chrom_pos.second))
        return true;
    return false;

}

//移动chrom:pos <string> 的位置
//n 为正向右移动，为负向左移动 
void move_breakp(string& break_p1,int n){
    if(n == 0) return;
    vector<string> v_b1 = split(break_p1,":");
    int pos = string2number<int>(v_b1[1]) + n;
    break_p1 = v_b1[0] + ":" + to_string(pos);
}


//根据断裂点在融合基因的5‘ 端还是 3’ 端，最终调整断裂点的位置
//@break_p1: chrom<string>:pos<int> 记录了第一个断裂点；
//@chain: soft-clip seq 比对到基因组的信息；
//@map_kmer_pos: 事先构建好的index；
//@reference: 用于获取基因组序列；
//@soft_clip_seq: 待检索的序列；
//@match_seq: read soft-clip seq 旁边和参考基因组一至的序列（默认取10bp)
//@chrom_map: 染色体号对应表
//@up: 代表断裂点在基因的上游(基因的5‘端）；
//@map_state  0: 下游正链；1：下游负链；2:上游正链；3:上游负链
position break_point_get(string& break_p1,Kmer_Chain chain,fastqReader& reference,const string& soft_clip_seq,string& match_seq,map<string,unsigned>& map_chrom_offset,vector<string>& v_chroms,int map_state){
    int error_hanming = 3;
    position chrom_pos = int2position(chain.pos,map_chrom_offset,v_chroms,0,v_chroms.size()-1);
    // cout << "Get Pos From: " << chain.pos << " To: " << chrom_pos.first <<":" << chrom_pos.second << "\t"
        // << match_seq << "\t" << soft_clip_seq <<"\n";
    int mlen = (int)match_seq.size();
    if (map_state == 0 || map_state == 3){
        string query_seq = soft_clip_seq.substr(0,chain.offset);
        string refer_seq = reference.fetch(chrom_pos.first,chrom_pos.second-chain.offset,chrom_pos.second-1).seq();
        reverse(refer_seq.begin(),refer_seq.end());
        reverse(query_seq.begin(),query_seq.end());
        int n = break_point_search(refer_seq,query_seq,error_hanming);
        chrom_pos.second -= n;
        // 断点附近可mapping 位置检查是否需要再移动
        if(!comparePos(break_p1,chrom_pos) && n == chain.offset){
            string refer2 = reference.fetch(chrom_pos.first,chrom_pos.second-mlen,chrom_pos.second-1).seq();
            int n = 0;
            string query_seq;
            if(map_state == 0) query_seq = match_seq;
            else query_seq = reverse_complete(match_seq);
            for(int i = mlen-1 ;i >= 0 ; i--){
                if (query_seq[i] == refer2[i]) n++;
                else break;
            }
            if(map_state == 0){
                chrom_pos.second -= n;
                move_breakp(break_p1,-n);
            } else {
                chrom_pos.second -= n;
                move_breakp(break_p1,n);
            }
                       
        }    
    }
    else {
        if (chain.offset == soft_clip_seq.size() - 12) {
            chrom_pos.second += 11;
        }
        else {
            int check_len = soft_clip_seq.size() - chain.offset - 12;
            string query_seq = soft_clip_seq.substr(chain.offset+12,check_len);
            string refer_seq = reference.fetch(chrom_pos.first,chrom_pos.second+12,chrom_pos.second+12+check_len-1).seq();
            int n = break_point_search(refer_seq,query_seq,error_hanming);
            chrom_pos.second += (n + 11);
            if (n != check_len) return chrom_pos;
        }
        if(!comparePos(break_p1,chrom_pos)){
            string refer2 = reference.fetch(chrom_pos.first,chrom_pos.second+1,chrom_pos.second+mlen).seq();
            int n = 0;
            string query_seq;
            if(map_state == 1) query_seq = reverse_complete(match_seq);
            else query_seq = match_seq;
            for(int i = 0;i<mlen;i++){
                if(query_seq[i]==refer2[i]) n++;
                else break;
            }
            if (map_state == 1) {
                chrom_pos.second += n;
                move_breakp(break_p1,-n);
            } else {
                chrom_pos.second += n;
                move_breakp(break_p1,n);
            }
        }      
    }
    return chrom_pos;
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
//@match_seq: read soft-clip seq 旁边和参考基因组一至的序列（默认取10bp)
//@chrom_map: 染色体号对应表
//@isdown: soft-clip seq 在下游；
//@rna: rna 序列的比对和DNA 有区别；
//return : 断裂位点，pair的第二位数， 第一位记录  -1: 没有找到位点；0: 下游正链；1：下游负链；2:上游正链；3:上游负链
//soft seq处于下游时，头部完全匹配另外的位置，往前还可以找到的mapping bases, 此时修改 break_p1值；
pair<position,int> map_soft_clip_seq(string& break_p1,unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,fastqReader& reference,string& soft_clip_seq,string& match_seq,map<string,unsigned>& map_chrom_offset,vector<string>& v_chroms,bool isdown,bool rna)
{
    if(rna) {
        if(soft_clip_seq.size() > 48) soft_clip_seq = soft_clip_seq.substr(0,48);
    }
    pair<position,int> result = make_pair(make_pair("0",0),-1); 
    if(soft_clip_seq.size() < 18) return result; 
    if(!isdown) reverse(soft_clip_seq.begin(),soft_clip_seq.end()); // soft-clip seq 在上游的read ,在pilup 时用反向的read，所以要转换回来。  
    mop map_kmer_offset_forward = extract_kmer(soft_clip_seq,12);
    vector<unsigned int> v_kmer_forward = select_kmer(map_kmer_pos,map_kmer_offset_forward,10000);

    string soft_clip_seq_rp = reverse_complete(soft_clip_seq);
    mop map_kmer_offset_reverse = extract_kmer(soft_clip_seq_rp,12);
    vector<unsigned int> v_kmer_reverse = select_kmer(map_kmer_pos,map_kmer_offset_reverse,10000);

    int error_golbal = 3;

    if(isdown){
        // 下游正链匹配
        Kmer_Chain pos_global_forward = global_align(map_kmer_pos,map_kmer_offset_forward,v_kmer_forward,error_golbal) ;
        // 下游负链匹配
        reverse(v_kmer_reverse.begin(),v_kmer_reverse.end());
        Kmer_Chain pos_global_reverse = global_align(map_kmer_pos,map_kmer_offset_reverse,v_kmer_reverse,error_golbal);
        // cout << "down! \n";
        if(pos_global_forward.chainNum >= pos_global_reverse.chainNum && pos_global_forward.chainNum >=2 ){            
            result.first = break_point_get(break_p1,pos_global_forward,reference,soft_clip_seq,match_seq,map_chrom_offset,v_chroms,0);
            result.second = 0;
        } else if(pos_global_reverse.chainNum >=2)
        {
            result.first = break_point_get(break_p1,pos_global_reverse,reference,soft_clip_seq_rp,match_seq,map_chrom_offset,v_chroms,1);
            result.second = 1;            
        }  
        // cout << "Mapping Finsh!\n";     
    } else {
        reverse(v_kmer_forward.begin(),v_kmer_forward.end());
        Kmer_Chain pos_global_forward = global_align(map_kmer_pos,map_kmer_offset_forward,v_kmer_forward,error_golbal);
        Kmer_Chain pos_global_reverse = global_align(map_kmer_pos,map_kmer_offset_reverse,v_kmer_reverse,error_golbal);

        // cout << "Mapping Finsh!\n";
        if(pos_global_forward.chainNum >= pos_global_reverse.chainNum && pos_global_forward.chainNum >=2){
            result.first = break_point_get(break_p1,pos_global_forward,reference,soft_clip_seq,match_seq,map_chrom_offset,v_chroms,2);
            result.second = 2;
        }else if(pos_global_reverse.chainNum >= 2){
            result.first = break_point_get(break_p1,pos_global_reverse,reference,soft_clip_seq_rp,match_seq,map_chrom_offset,v_chroms,3) ;
            result.second = 3;
        }
    }
    cout << soft_clip_seq << ": " << result.first.first <<  ":" << result.first.second << "\n";
    return result;   
}




//处理 bam1_t *aln 对象（bam 比对结果对象），找到 split-reads 和 discordant reads
//处理含soft-clip 的reads 或 insert size 大于 1000 的reads , mapQ 大于等于 15
//该function 会修改 map_dcp 和 map_split_read 这两个表记录了符合对应条件的reads 信息。
//记录每条序列的气势位点和insertSize;
//@aln bam 的比对记录；
//@bam_header bam的头文件；
//@map_dcp: 散列表用于记录discordant reads;
//@map_split_read: 散列表用于记录split reads;
bool deal_aln(const bam1_t* aln,const bam_hdr_t* bam_header,unordered_map<string,fusionPos>& map_dcp,unordered_map<string,vector<Piled_reads>>& map_split_read,unordered_map<string,vector<string>>& map_alt_split,unordered_map<string,string>& map_transcript){
    string qname = getName(aln);
    string chrom = getChrom(aln,bam_header);
    string nchrom = getNchrom(aln,bam_header);
    string cigar = getCigar(aln);
    string seq = getSeq(aln);
    int pos = getPos(aln);
    int npos = getNpos(aln);
    int flag = getFlag(aln);
    int mapQ = getMapq(aln);
    int isize = getIsize(aln);

    // cout << qname << "\t" << cigar << "\t" << chrom << ":" << pos <<"\t" << nchrom <<":" << npos << "\n";

    if(s_chroms.find(chrom) == s_chroms.end()) return false;
    if(mapQ < 15 || flag & 1024 ) return true;
    else if(parse_split_read(chrom,pos,seq,cigar,qname,map_split_read, map_alt_split,map_transcript)) {
        // cout << "Deal split-reads: " << qname << "\t" << cigar << "\t" << isize << "\t" << chrom << ":" << pos << "\n";
        return true;
    }
    else if(abs(pos-npos) > 10000 || ( chrom != nchrom && nchrom != "NA")) {
        // cout << "Deal diso pairs: " << qname << "\t" << cigar << "\t" << isize << "\t" << chrom << ":" << pos << "\n";
        parse_discordant_read(qname,flag,chrom,pos,nchrom,npos,cigar,map_dcp);
    }

    return true;
}

//对split-read soft clip reads 进行比对找到符合的断点；
//pair<position,int> result = map_soft_clip_seq(map_kmer_pos,reference,read,map_chrom_offset,v_chroms,true)
//void combine_split_reads(string& break_p1,pair<string,int>& break_p2,int num,map<string,fusion>& map_fusion, map<string,int>& map_c2i);
//mapping soft-clip seq 到新的断裂点，并将信息记录到散列表 map_fusion 中；
//@rna: rna 的比对和DNA 不一致
void map_split_reads(string& break_p1,Piled_reads& prd,unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,fastqReader& reference,map<string,unsigned>& map_chrom_offset,vector<string> v_chroms,map<string,fusion>& map_fusion,bool rna){
    cout << "Mapping: " << break_p1 << "\t" << prd.seq << "\n";
    pair<position,int> result = map_soft_clip_seq(break_p1,map_kmer_pos,reference,prd.seq,prd.match_seq,map_chrom_offset,v_chroms,prd.down,rna);
    cout << "first Map " << break_p1 << "\t" << result.first.first << ":" << result.first.second << "\n";
    if(result.second == -1) return;
    combine_split_reads(break_p1,result,prd,map_fusion,map_c2i);
}

//二分法寻找断裂点距离最近的外显子的位置。
//没能在外显子附近找到断裂点，则返回 <NA,0>;
//在外显子左侧找到，返回<gene,-offset>;
//在右侧找到，返回<gene,offset>;
//默认使用的flank=10;
//return string: 基因名（转录本号+外显子号， SGIP1:NM_001308203:e2）； int 与外显子距离; false 表示在外显子外，true 表示在外显子中；
pair<pair<int,int>,bool> find_near_extron(const position& break_p,vector<pair<int,int>>& v_extron,int s,int e,int flank=10){
    if(s > e ){
        return make_pair(make_pair(-1,0),false);
    }
    int n = (s + e) / 2;
    if(break_p.second <= v_extron[n].first ){
        if(v_extron[n].first - break_p.second <= flank) return make_pair(make_pair(n,break_p.second-v_extron[n].first),false);
        else return find_near_extron(break_p,v_extron,s,n-1,flank);
    } else if (break_p.second >= v_extron[n].second){
        if(break_p.second - v_extron[n].second <= flank) return make_pair(make_pair(n,break_p.second - v_extron[n].second),false);
        else return find_near_extron(break_p,v_extron,n+1,e,flank);
    } else {
        int offset = min(break_p.second-v_extron[n].first,v_extron[n].second-break_p.second);
        return make_pair(make_pair(n,offset),true);
    }
}

pair<pair<string,int>,bool> find_extron(const position& break_p,vector<string>& refgenes,unordered_map<string,Extron>& mge){
    pair<pair<string,int>,bool> res = make_pair(make_pair("NA",0),false);
    pair<pair<int,int>,bool> ores_p = make_pair(make_pair(-1,0),false);
    auto &res_p = ores_p.first;
    string gene;
    int extron_num = 0;
    for(auto& ref: refgenes){
        if(mge[ref].v_pos[0].first > break_p.second + 10 || mge[ref].v_pos.back().second < break_p.second - 10 ) continue;
        ores_p = find_near_extron(break_p,mge[ref].v_pos,0,mge[ref].v_pos.size()-1);
        if (res_p.first != -1) {
            gene = ref;
            extron_num = mge[ref].v_extrons[res_p.first];
            break;
        }
    }
    if(res_p.first != -1){
        gene = gene + ":e" + to_string(extron_num) + ":" + to_string(res_p.second);
        res = make_pair(make_pair(gene,res_p.second),ores_p.second);
    }
    return res;
}



//修正基因上下游关系；根据基因的转录方向，融合的位置
//调换了位置返回 true, 否则 false 
bool md_gene_pos(fusion& fu,string& p1_gene,string& p2_gene,unordered_map<string,Extron>& mge){
    bool tag = false;
    char p1dr = '-',p2dr = '-';
    if(orf_dr(p1_gene,mge)) p1dr = '+';
    if(orf_dr(p2_gene,mge)) p2dr = '+';
    if(fu.mpdr.first == '+'){
        if(p1dr == '-' && p2dr == '-'){
            swap(fu.p1,fu.p2);
            swap(fu.p1n,fu.p2n);
            swap(fu.split_seqs.first,fu.split_seqs.second);
            fu.genes = make_pair(p2_gene,p1_gene);
            fu.genedr = make_pair(p2dr,p1dr);
            tag = true;
        } else{
            fu.genes = make_pair(p1_gene,p2_gene);
            fu.genedr = make_pair(p1dr,p2dr);
        }
        
    } else {
        if((p1dr == '-' && p2dr == '+' && fu.rpdr == make_pair('u','u')) ||
            p1dr == '+' && p2dr == '-' && fu.rpdr == make_pair('d','d')){
            swap(fu.p1,fu.p2);
            swap(fu.p1n,fu.p2n);
            swap(fu.split_seqs.first,fu.split_seqs.second);
            fu.genes = make_pair(p2_gene,p1_gene);
            fu.genedr = make_pair(p2dr,p1dr);
            tag = true;
        } else {
            fu.genes = make_pair(p1_gene,p2_gene);
            fu.genedr = make_pair(p1dr,p2dr);
        }
    }
    return tag;
}


//逆向比较两个长度相同序列的碱基，返回相同碱基的数量，以修正断裂点
int cmpSeq(string& seq1,string& seq2){
    int mis = 1;
    reverse(seq1.begin(),seq1.end());
    reverse(seq2.begin(),seq2.end());
    int n = 0,m=0;
    for(unsigned i = 0;i<seq1.size();i++){
        if(seq1[i] == seq2[i]) n++;
        else {
            m++;
            if(m>mis)
                return n;
            else n++;
        }
    }
    return n;
}

//修正RNA的断裂点位置，使得RNA的断裂点尽量处于外显子内；
void md_rna_breakpos(fusion& fu,map<string,vector<string>>& m_extron,fastqReader& reference,unordered_map<string,Extron>& mge){
    auto op1_extron = find_extron(fu.p1,m_extron[fu.p1.first],mge);
    auto op2_extron = find_extron(fu.p2,m_extron[fu.p2.first],mge);
    auto &p1_extron = op1_extron.first;
    auto &p2_extron = op2_extron.first;
    if(md_gene_pos(fu,p1_extron.first,p2_extron.first,mge)) swap(p1_extron,p2_extron);
    if(p1_extron.first == "NA" || p2_extron.first == "NA" || max(p1_extron.second,p2_extron.second) > 10) return;
    string seq1,seq2;
    int len = min(abs(p1_extron.second),abs(p2_extron.second));
    cout << "Move to Extron: " << fu.p1.first << ":" <<fu.p1.second << "-" << fu.p2.first << ":" <<fu.p2.second;
    if(fu.mpdr.first == '+'){
        if(!op1_extron.second){
            if(p1_extron.second < 0){
            // len = -p1_extron.second;
            seq1 = reference.fetch(fu.p1.first,fu.p1.second,fu.p1.second + len -1).seq();
            seq2 = reference.fetch(fu.p2.first,fu.p2.second+1,fu.p2.second + len).seq();
            len = cmpSeq(seq1,seq2);
            fu.p1.second += len;
            fu.p2.second += len;
        
            } else if(p1_extron.second > 0){
                // len = p1_extron.second;
                seq1 = reference.fetch(fu.p1.first,fu.p1.second - len + 1,fu.p1.second).seq();
                seq2 = reference.fetch(fu.p2.first,fu.p2.second - len,fu.p2.second - 1).seq();
                len = cmpSeq(seq1,seq2);
                fu.p1.second -= len;
                fu.p2.second -= len;
            }
        } else if(!op2_extron.second) {
            if (p2_extron.second < 0){
                // len = -p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second,fu.p2.second + len - 1).seq();
                seq2 = reference.fetch(fu.p1.first,fu.p1.second+1,fu.p1.second + len).seq();
                len = cmpSeq(seq1,seq2);
                fu.p1.second += len;
                fu.p2.second += len;
            
            } else if(p2_extron.second > 0){
                // len = p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second - len + 1,fu.p2.second).seq();
                seq2 = reference.fetch(fu.p1.first,fu.p1.second - len,fu.p1.second - 1).seq();
                len = cmpSeq(seq1,seq2);
                fu.p1.second -= len;
                fu.p2.second -= len;
                
            }

        }
        
    } else {
        if(!op1_extron.second){
            if(p1_extron.second < 0){
                // len = -p1_extron.second;
                seq1 = reference.fetch(fu.p1.first,fu.p1.second,fu.p1.second + len -1).seq();
                seq2 = reference.fetch(fu.p2.first,fu.p2.second - len,fu.p2.second - 1).seq();
                seq1 = reverse_complete(seq1);
                len = cmpSeq(seq1,seq2);
                fu.p1.second += len;
                fu.p2.second -= len;
                
            }else if(p1_extron.second > 0){
                // len = p1_extron.second;
                seq1 = reference.fetch(fu.p1.first,fu.p1.second - len + 1,fu.p1.second).seq();
                seq2 = reference.fetch(fu.p2.first,fu.p2.second+1,fu.p2.second + len).seq();
                seq1 = reverse_complete(seq1);
                len = cmpSeq(seq1,seq2);
                fu.p1.second -= len;
                fu.p2.second += len;
        
            }
        } else if (!op2_extron.second){
            if(p2_extron.second < 0){
                // len = -p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second,fu.p2.second + len - 1).seq();
                seq2 = reference.fetch(fu.p1.first,fu.p1.second - len,fu.p1.second - 1).seq();
                seq1 = reverse_complete(seq1);
                len = cmpSeq(seq1,seq2);
                fu.p2.second += len;
                fu.p1.second -= len;
                
            } else if(p2_extron.second > 0){
                // len = p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second - len + 1,fu.p2.second).seq();
                seq2 = reference.fetch(fu.p1.first,fu.p1.second+1,fu.p1.second + len).seq();
                seq1 = reverse_complete(seq1);
                len = cmpSeq(seq1,seq2);
                fu.p2.second -= len;
                fu.p1.second += len;
                
            }                                    
        }
    }
    cout << " Move to Extron: " << fu.p1.first << ":" <<fu.p1.second << "-" << fu.p2.first << ":" <<fu.p2.second;
    cout << "\t" << p1_extron.second << ";"<<p2_extron.second<<"\t" << seq1 << "\t" << seq2 << "\n";
}



//vector 去重
template<class T>
void sort_unique(vector<T>& vp){
    sort(vp.begin(),vp.end());
    auto it = unique(vp.begin(),vp.end());
    vp.resize(distance(vp.begin(),it));
    return;
}

//根据cigar 找不同的模板数量
int templateNum(vector<string>& vcigars){
    vector<int> mv;
    vector<int> sv;
    int n;
    char c;
    for(auto& cigar: vcigars){
        stringstream ss{cigar};
        while(ss>>n>>c){
            if(c=='M') mv.push_back(n);
            else if(c=='S') sv.push_back(n);
        }
    }
    sort_unique<int>(mv);
    sort_unique<int>(sv);
    return min(mv.size(),sv.size());
}


//整合相邻的融合突变，将证据加和，输出结果
//@vmf: 记录了相邻的融合突变集合
//@dis : 合并的距离
void combine_rna_fus(vector<pair<fusion,bool>>& vmf,int dis){
    if(vmf.size() == 0) return;
    for(unsigned i = 0; i < vmf.size();i++){
        if(vmf[i].second) continue;
        cout << "Combine: " << vmf[i].first.p1.first << ":" << vmf[i].first.p1.second << "-" << vmf[i].first.p2.first << ":"
            <<vmf[i].first.p2.second << "\t";
        for(unsigned j = i+ 1;j < vmf.size();j++){
            if(vmf[i].first.p2.first == vmf[j].first.p2.first && abs(vmf[i].first.p2.second - vmf[j].first.p2.second) <= dis && abs(vmf[i].first.p1.second - vmf[j].first.p1.second) <= dis){
                vmf[i].first.p1n += vmf[j].first.p1n;
                vmf[i].first.p2n += vmf[j].first.p2n;
                vmf[i].first.dcp += vmf[j].first.dcp;
                for(auto s:vmf[j].first.vp) vmf[i].first.vp.push_back(s);
                for(auto s:vmf[j].first.p1_cigars) vmf[i].first.p1_cigars.push_back(s);
                for(auto s:vmf[j].first.p2_cigars) vmf[i].first.p2_cigars.push_back(s);
                if(vmf[i].first.split_seqs.first == "") vmf[i].first.split_seqs.first = vmf[j].first.split_seqs.first;
                if(vmf[i].first.split_seqs.second == "") vmf[i].first.split_seqs.second = vmf[j].first.split_seqs.second;
                vmf[j].second = true;
                cout << vmf[j].first.p1.first << ":" << vmf[j].first.p1.second << "-" << vmf[j].first.p2.first << ":"
                    <<vmf[j].first.p2.second << "\t";
            }
        }
        cout << "\n";
        auto& f = vmf[i].first;
        // cout<< f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" << f.p1n << ":"<<f.p2n<<":" << f.dcp  << "\n";
        if (min(f.p1n,f.p2n) >= 1)
        {
            if( (f.p1.first != f.p2.first || abs(f.p1.second - f.p2.second) > 5000) ){
                sort_unique<string>(f.vp);
                int p1c = templateNum(f.p1_cigars);
                int p2c = templateNum(f.p2_cigars);
                cerr << f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" 
                    << f.mpdr.first << ";" << f.mpdr.second << "\t" 
                    << f.genedr.first << ";" << f.genedr.second << "\t"
                    << f.rpdr.first << ';' << f.rpdr.second << "\t"
                    << f.genes.first << ";" << f.genes.second << "\t"
                    << f.p1n << ":"<<f.p2n<<":" << f.dcp  <<  "\t" 
                    << p1c << ":" << p2c << ":" << f.vp.size() << "\t"
                    << f.split_seqs.first << ":" << f.split_seqs.second << "\n";
            }
        }
    }
}

//对融合进行排序
bool cmpFu(pair<string,fusion>& aa,pair<string,fusion>& bb){
    auto &a = aa.second;
    auto &b = bb.second;
    if(a.p1.first > b.p1.first) return true;
    else if(a.p1.first < b.p1.first) return false;
    else if(a.p1.second > b.p1.second) return true;
    else return false;
}

//RNA Mapping 时会存在soft-clip 断裂点会根据外显子位置修正，导致我们的断裂点修正存在偏差；
//针对这种情况，可以合并 双断裂点在 d (10) 范围内的，融合变异；
//@map_fusion 记录了所有检测到的潜在融合变异，且按顺序输出；
//@dis 定义距离在 10 范围内的突变可以整合；
//输出符合需求的融合突变
void rna_print_sv(map<string,fusion>& map_fusion,int dis=10){
    vector<pair<string,fusion>> vfs(map_fusion.begin(),map_fusion.end()); 
    sort(vfs.begin(),vfs.end(),cmpFu);
    vector<pair<fusion,bool>> vmf;
    fusion sf = vfs.begin()->second;
    for(auto& mf : vfs){
        auto &f = mf.second;
        cout<< "未合并前： " << mf.first << "\t" << f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" << f.mpdr.first << "\t" << f.mpdr.second << "\t"  << f.p1n << ":"<<f.p2n<<":" << f.dcp  << "\n";
        if(mf.second.p1.first == sf.p1.first && abs(mf.second.p1.second - sf.p1.second) <= dis){
            vmf.push_back(make_pair(mf.second,false));
            sf = mf.second;
        } else {
            combine_rna_fus(vmf,dis);
            sf = mf.second;
            vmf = {make_pair(sf,false)};
        }
    }
    combine_rna_fus(vmf,dis);
}

//输出非整合附近断裂位点融合结果，适用于DNA数据；
//@map_fusion 记录了所有检测到的潜在融合变异，且按顺序输出；
void dna_print_sv(map<string,fusion>& map_fusion){
    for(auto mf:map_fusion){
        auto &f = mf.second;
        cout << f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" << f.mpdr.first << "\t" << f.mpdr.second << "\t" << f.p1n << ":"<<f.p2n<<":" << f.dcp  << "\n";
        if (f.p1n + f.p2n >= 6 && min(f.p1n,f.p2n) >= 2)
        {
            if(f.p1.first != f.p2.first || abs(f.p1.second - f.p2.second) > 5000) 
                cerr << f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" << f.mpdr.first << "\t" << f.mpdr.second << "\t" << f.p1n << ":"<<f.p2n<<":" << f.dcp  << "\n";
        }            
    }
}

//获得特定转录本的剪切位点信息
unordered_map<string,string> getTranscript(string infile){
    ifstream ifs{infile};
    unordered_map<string,string> map_transcript;
    string chrom;
    int start,end;
    string extron;
    while(ifs >> chrom >> start >> end >> extron){
        map_transcript[chrom + ":" + to_string(start)] = extron;
        map_transcript[chrom + ":" + to_string(end)] = extron;
    }
    return map_transcript;
}

int main(int argc,char* argv[]){
    if(argc < 7){
        cout << argv[0] << "\tbam\treference\t12-kmer-index\trefseq_trans.exon.bed\ttranscript-altsplit\t0:DNA(1:RNA)\n";
        return -1;
    }
    
    string bamFile = argv[1];
    string fasta = argv[2];
    string kmerFile = argv[3];
    string extronFile = argv[4];
    string transFile = argv[5];
    string mode = argv[6];
    bool rna = false;
    if (mode == "0") rna = false;
    else if(mode == "1") rna = true;
    else {
        cout << "Mode need be RNA or DNA!";
        return -1;
    }

    map<string,unsigned> map_chrom_offset = make_chrom_offset(fasta);
    fastqReader reference = fastqReader(fasta);
    BamReader brd = BamReader(bamFile);
    unordered_map<string,vector<string>> map_alt_split;
    unordered_map<string,fusionPos> map_dcp;
    unordered_map<string,vector<Piled_reads>> map_split_read;
    map<string,fusion> map_fusion;
    unordered_map<string,Extron> mge;

    //获得外显子、基因信息
    map<string,vector<string>> m_extron = getExtronInfo(extronFile,mge);
    // 获得要进行可变剪切检测的相关转录本信息
    unordered_map<string,string> map_transcript = getTranscript(transFile);

    cout << "开始获取融合支持reads ...\n";
    // 寻找要处理的reads
    while(brd.next()){
        if (!deal_aln(brd.aln,brd.bam_header,map_dcp,map_split_read,map_alt_split,map_transcript)) break;       
    }
    cout << "融合支持reads获取结束！\n";

    ifstream kmer_index {kmerFile};
    unordered_map<unsigned int,vector<unsigned>> map_kmer_pos;
    boost::archive::binary_iarchive iarch(kmer_index);
    iarch >> map_kmer_pos;
    kmer_index.close();

    //mapping split-reads 把相同断点的 split-reads 整合到一起。
    cout << "把相同断点的 split-reads 整合到一起。\n";
    for(auto& sr:map_split_read){
        string break_p1 = sr.first;
        cout << "Piled reads: " << break_p1 << "\n";
 
        for(auto& fu:sr.second){
            string break_p1 = sr.first;
            cout << "FF: " << break_p1 << "\t" <<fu.seq << "\t" << fu.count << "\t" << fu.down << "\n";
            map_split_reads(break_p1,fu,map_kmer_pos,reference,map_chrom_offset,v_chroms,map_fusion,rna);
        }
    }

    //将discordant reads 的支持信息也增加到split-reads 信息中
    cout << "将discordant reads 的支持信息也增加到split-reads 信息中\n";
    unordered_map<string,vector<string>> map_fus_dcps = combine_discordant_reads(map_dcp);
    combine_fusion(map_fusion,map_fus_dcps,map_dcp);

    //修正基因上下游关系
    for(auto& mfu:map_fusion){
        auto& fu = mfu.second;
        md_rna_breakpos(fu,m_extron,reference,mge);
    }

    //展示查找到的融合信息
    cerr << "break_pos1\tbreak_pos2\tsoft_map_dr(p1);soft_map_dr(p2)\tgene_p1_dr;gene_p2_dr\tsupport_p1_dr;support_p2_dr\t"
        << "gene_p1;gene_p2\tsplit_p1:split_p2:discordant reads\tp1_cigars;p2_cigars;total_reads\n"; 
    if(!rna) dna_print_sv(map_fusion);
    else rna_print_sv(map_fusion);

    cerr << "\nAlterNative SplitInfos\n";
    for(auto& a : map_alt_split){
        cerr << a.first << "\t" << a.second.size() <<"\n";
    }

}


