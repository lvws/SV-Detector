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
// 基因组位点的表示形式
using position = pair<string,int>; 
// 记录read每个kmer对应的reads上的偏离值
using mop = map<unsigned int,int>; 
// 记录断点Mapping 的位置，以及mapping 的形式；
using mstat = pair<position,int> ; 

// 方便将基因组绝对位置转化为 chrom:pos
vector<string> v_chroms = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
                                "19","20","21","22","X","Y","GL000228.1"};

// map<string,int> map_c2i = {{"1",1},{"2",2},{"3",3},{"4",4},{"5",5},{"6",6},{"7",7},{"8",8},
//     {"9",9},{"10",10},{"11",11},{"12",12},{"13",13},{"14",14},{"15",15},{"16",16},{"17",17},{"18",18},
//     {"19",19},{"20",20},{"21",21},{"22",22},{"X",23},{"Y",24}};

unordered_set<string> s_chroms(v_chroms.begin(),v_chroms.end());

//数据结构，记录kmer连锁的信息
//@pos : unsigned 基因组绝对位置
//@offset: int kmer距离断裂点长度
//@chainNum: int 连锁的kmer 数量
//@sa: bool 是否直接来自SA tag
struct Kmer_Chain
{
    unsigned pos; //基因组位置
    int offset; //接近断裂点kmer 距离read
    int chainNum;   //连锁的kmer 数量
    bool sa = false; //默认不是
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

//vector 去重
template<class T>
void sort_unique(vector<T>& vp){
    sort(vp.begin(),vp.end());
    auto it = unique(vp.begin(),vp.end());
    vp.resize(distance(vp.begin(),it));
    return;
}

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
    map_chrom_offset["GL000228.1"] = map_chrom_offset["MT"]; // 增加这个contig
    return map_chrom_offset;
}

//将染色体位点转换为绝对位置。
//@chrom : 染色体号
//@pos: 位置
unsigned position2int(string chrom,int pos,map<string,unsigned>& map_chrom_offset){
    if (map_chrom_offset.find(chrom) == map_chrom_offset.end())
        return 0;
    return map_chrom_offset[chrom] + pos;
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
            //cout << pos << " " << v_chroms[s] << " " << map_chrom_offset[v_chroms[s]] << "\n";
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


//找某一断点处的，断点类型计数 （10bp 上下）
//@pos 要计算的融合断点
//@map_breaks 记录每个断点的不同类型的reads的数量（primaray alignment）
//@flank 允许的偏移空间
//@return 一串数字，记录断点附近的不同潜在融合reads的数量
vector<int> break_point_info(pair<string,int>& pos,map<pair<string,int>,vector<int>>& map_breaks,unordered_map<string,int>& map_spreads_count,int flank=10){
    vector<int> vp_counts;
    for(int s=pos.second-flank;s <= pos.second+flank;s++){
        pair<string,int> check_pos = make_pair(pos.first,s);
        if (map_breaks.find(check_pos) != map_breaks.end()){
            // cout << check_pos.first << ':' << check_pos.second << '\n'; 
            // sort_unique<string>(map_breaks[check_pos]);
            for(auto i:map_breaks[check_pos]) {
                vp_counts.push_back(i);
            }
        }
    }
    return vp_counts;
}

//计算该断点处的， 断点reads总数,本融合断点数，断点reads类型总数、本断点数量排名
//@vp_counts 上一步计算出的附近断裂位点信息；
//@pu_count 本断裂点支持reads数；
//@m 本断裂点合并的reads 类型数；
//@return  按顺序记录 【断点reads总数，本融合断点数，断点reads类型总数，本断点数量排名】
vector<int> break_point_info2(vector<int>& vp_counts,int pu_count,int m){
    int total_count = 0;
    int total_type = vp_counts.size() - m ;
    int n = 1;
    for(int i :vp_counts){
        total_count += i;
        if (i > pu_count) n++;
    }
    vector<int> vb_infos = {total_count,pu_count,total_type,n};
    return vb_infos;
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
    else if(pos + dis < v_pos[m]) e = m - 1;
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


//有时存在soft-clip 较长（RNA 只选择了前48bp），但是所选择的序列在基因组中特异性太差，不能有效形成kmer-chain 被过滤。 
//所以考虑使用全长的soft-clip，但是经过（特异性筛检后）最多保留前4个kmer（防止RNA 中某些外显子过短引起的错误）
//检索kmer在index中的位置 vector, 过滤超多位点的kmer。将剩余的kmer按顺序由左至右排序。
//@map_kmer_pos: 事先构建好的index
//@map_kmer_offset: read中kemr 对应的偏离值；
//@threshold: 超过固定数量位点的kmer 就不进行后续的检测；
//@isdown: 表示soft 在reads 下游，取kmer 时要顺向取;
//return : kmer 集合，按read的 offset 排序，由低到高。
vector<unsigned int> select_kmer(unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,mop& map_kmer_offset,int threshold,bool isdown,bool rna,bool positive){
    //map<unsigned int,unsigned int> kmer_pos_size;
    vector<pair<unsigned,unsigned>> v_kmer_size;
    vector<unsigned int> v_kmer;
    for(auto i:map_kmer_offset){
        if(map_kmer_pos[i.first].size() > threshold) {
            // cout << map_kmer_pos[i.first].size() << "\t" ;
            //kmer_pos_size[i.first] = map_kmer_pos[i.first].size();
            if(map_kmer_pos[i.first].size() < 10000)
                v_kmer_size.push_back(make_pair(i.first,map_kmer_pos[i.first].size()));
            continue; 
        }
        v_kmer.push_back(i.first);
    }

    //对于某些kmer 唯一性差，补充到三个kmer(起码有一个Kmer 是比较唯一的)
    // int n = 3 - v_kmer.size();
    // if(n > 0 && n < 3){
    //     sort(v_kmer_size.begin(),v_kmer_size.end(),[&](pair<unsigned,unsigned> a,pair<unsigned,unsigned> b){return a.second < b.second;}); //优先取唯一性好的kmer;
    //     for(auto kps :v_kmer_size){
    //         v_kmer.push_back(kps.first);
    //         n = n - 1;
    //         if (n == 0) break;
    //     }
    // }
    //RNA 存在外显子短的情况，为避免跨越外显子，使用靠近断点的、较短的序列作为kmer
    if ((!isdown && positive) || (isdown && !positive)) sort(v_kmer.begin(),v_kmer.end(),[&](unsigned a,unsigned b){return map_kmer_offset[a] > map_kmer_offset[b];});
    else sort(v_kmer.begin(),v_kmer.end(),[&](unsigned a,unsigned b){return map_kmer_offset[a] < map_kmer_offset[b];});
    if(!rna) return v_kmer;
    vector<unsigned int> v2 = {v_kmer.begin(),v_kmer.begin()+min(4,(int)v_kmer.size())};
    return v2;
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
//@map_state: 预计的mapping 方式，用来选择SA mapping.
//@least: 记录最少要连锁的长度。
//@sct: 记录该断点的支持reads 数；uniq count 少于一条不Map（有可能0），total count 少于三条不进行SA 检索；
//return: pair数据，记录能找到最靠近断裂点的kmer ; 
// 其在基因组上的位置，以及在read 上的偏离值（由0开始），连锁分值, 若没能找到连锁，返回值的偏离为 -1，连锁分值0.
vector<Kmer_Chain> global_align(unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,mop& map_kmer_offset,vector<unsigned int> v_kmer,
int err,vector<pair<unsigned,bool>> v_sa,int sct,int usct, int map_state,int least = 2){
    vector<Kmer_Chain> results,results2;
    Kmer_Chain result = {0,-1,0};
    results2.push_back(result);
    int score = max(2,(int)v_kmer.size()/2) ; //至少要求 一半的 kmer 可以连锁；
    int bias = 50; // supplementary 位点和重新mapping 的位点间允许的偏差；
    least = max(2,least); // 至少保证2 长度；
    //cout << "Kmer Chain " << "kemr-vector-len: " << v_kmer.size() << "\n";

    if(sct < 3)
        v_sa = {};
    else if (usct < 1)
        least = 99;

    //cout << "Kmer Chain " << "kemr-vector-len: " << v_kmer.size() << "\tleast:"<<least<< "\n";
    if (v_kmer.size() >= least) {
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
                    //cout << "sites: " << v_pos_j.size();
                    if (distance_check(pos,v_pos_j,0,(int)(v_pos_j.size())-1,dis,3)) m ++;
                    else k++;
                    //cout << "start offset: " << offset << '\t' << "link offset: " << map_kmer_offset[kmer_j] << '\t' << "m: " << m << '\n';
                    if (k > e) break;
                }
                //cout << m << "\t" << score << "\n";
                if(m >= score) {
                    result = {pos,offset,m};
                    results.push_back(result);
                    score = m;
                    //cout << "Kmer vector size: " << v_kmer.size() - i
                    //    << " Kmer Index: "  << i <<" Kmer Offset: " << offset
                    //    <<" Pos: " << pos << " Score: " << m << "\n";
                }


            }   
            if(results.size() != 0){
                if(results[-1].chainNum > (v_kmer.size() - i - 1) && results[-1].chainNum > 1)
                    break;
            }     
        }
    }

    
    vector<Kmer_Chain> v_tmp;
    map<pair<unsigned,bool>,bool> map_sa;
    bool sa_tag = false;
    // 优先采用SA tag 位置，防止重复碱基造成的mapping 位置不一致
    for(auto r:results){
        if(r.chainNum == score){
            // v_tmp.push_back(r);
            sa_tag = false;
            for(auto sa_pos : v_sa){
                //cout << r.pos << "\t" << sa_pos << "\n";
                if(sa_pos.first <= r.pos + bias && sa_pos.first + bias >= r.pos){
                    // map_sa[sa_pos] = true;
                    //v_tmp.push_back({sa_pos,0,score,true});
                    sa_tag = true;
                    break;
                } 
            }
            if(!sa_tag) v_tmp.push_back(r);
        }
    }

    // 太多的匹配位置只留SA的mapping。
    if(v_tmp.size() > 100) {
        v_tmp = {};
        for (auto& sa:map_sa) sa.second = false; 
    }

    for (auto sa_pos:v_sa){
        // MS
        if(sa_pos.second == false ){
            if(map_state == 0 || map_state == 3) continue;
        }
        // SM
        else {
            if(map_state == 1 || map_state == 2) continue;
        }
        // 只记录没有map到的SA 位置， （防止存在map 间隔造成的 不匹配）
        if(map_sa.find(sa_pos) == map_sa.end() ) v_tmp.push_back({sa_pos.first,0,score,true}); 
    }

    // if(v_tmp.size() > 100) return results2;
    if(v_tmp.size() == 0)
        v_tmp.push_back(result);
    
    return v_tmp;

} 

//动态规划算法计算两个序列之前的相似度，可解决错配、delins的影响
//默认设定match 1 ; mismatch -1; delins -2; 允许错配空间 5 bp;
//计算得到该匹配的最高得分。
int SW_map(string& seq1,string& seq2, int match = 1, int mismatch = -1, int delins = -2){
    map<pair<int,int>,pair<int,int>> map_path;
    map<pair<int,int>,int> map_score;
    int up = 0,left = 0 ,dd = 0,now = 0;
    int dis = 7;
    for (int i = 0 ; i < (int)seq1.size();i++){
        for(int j = max(i-dis,0); j < min((int)seq2.size(),i+dis);j++){
            // 检测当位点的值
            if(seq1[i] == seq2[j]) 
                now = match;
            else 
                now = mismatch;
            auto p_now = make_pair(i,j);
            auto p_dd = make_pair(i-1,j-1);
            auto p_up = make_pair(i-1,j);
            auto p_left = make_pair(i,j-1);
            // up socre
            if(i-1 >= 0) 
                up = map_score[p_up] + delins ;
            else 
                up = 0;
            if(j-1 >= 0)
                left = map_score[p_left] + delins ;
            else
                left = 0;
            if(i-1>=0 && j-1 >= 0)
                dd = map_score[p_dd] + now ;
            else
                dd = 0;
            
            if( dd > 0 && dd >= up && dd >= left){
                map_path[p_now] = p_dd;
                map_score[p_now] = dd;
            } 
            else if(up > 0 && up > dd && up >= left){
                map_path[p_now] = p_up;
                map_score[p_now] = up;
            }
            else if(left > 0 && left > dd && left > up){
                map_path[p_now] = p_left;
                map_score[p_now] = left;
            }
            else{
                map_path[p_now] = make_pair(-1,-1);
                map_score[p_now] = max(now,0);
            }
        }
    }
    vector<pair<pair<int,int>,int>> v_score;
    for (auto& i : map_score) v_score.push_back(i);
    sort(v_score.begin(),v_score.end(),[&](pair<pair<int,int>,int>& a,pair<pair<int,int>,int>& b){return a.second > b.second;});
    // 这注释的代码是，为了可视化比对效果；
    // auto chain = v_score.begin()->first;
    // int s = v_score.begin()->second;
    // int pre_i = chain.first ;
    // int pre_j = chain.second ;
    // for(;s > 0;){
    //     auto next_chain = map_path[chain];
    //     if (pre_i == next_chain.first ){
    //         cout << '*' << '-' << seq2[pre_j--] << '\t' << s  << '\n';
    //     }
    //     else if (pre_j == next_chain.second){
    //             cout << seq1[pre_i--] << '-' << "*\t" << s  << '\n';
    //     }
    //     else cout << seq1[chain.first] << '-' << seq2[chain.second] <<  "\t" <<  s <<  '\n';
    //     pre_i = next_chain.first;
    //     pre_j = next_chain.second;
    //     chain = next_chain;
    //     if (chain.first == -1) break;
    //     s = map_score[chain];
    // }
    return v_score.begin()->second;

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
            if(k==0 || n - m >= 3) m = n; //距离第一次出现错配，有多个正确匹配，可以刷新m值。
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
string move_breakp(string& break_p1,int n){
    string new_break_p1 = break_p1;
    if(n == 0) return new_break_p1;
    vector<string> v_b1 = split(break_p1,":");
    int pos = string2number<int>(v_b1[1]) + n;
    new_break_p1 = v_b1[0] + ":" + to_string(pos);
    return new_break_p1;
}


//检查两条序列是否一致（可能存在delins ; mismatch)


//根据断裂点在融合基因的5‘ 端还是 3’ 端，最终调整断裂点的位置
//在输出reMap的断点2的位置同时会生成新的断点1的信息
//@break_p1: chrom<string>:pos<int> 记录了第一个断裂点；
//@chain: soft-clip seq 比对到基因组的信息；
//@map_kmer_pos: 事先构建好的index；
//@reference: 用于获取基因组序列；
//@soft_clip_seq: 待检索的序列；
//@match_seq: read soft-clip seq 旁边和参考基因组一至的序列（默认取10bp)
//@chrom_map: 染色体号对应表
//@up: 代表断裂点在基因的上游(基因的5‘端）；
//@map_state  0: 下游正链；1：下游负链；2:上游正链；3:上游负链
// **** 代表soft seq ; ======代表 match seq;
// =======******   --->   *****======   map_state: 0 ; seq 正向map;
// =======******   --->   =======****   map_state: 1 ; seq 逆向map;
// *****========   --->   =======****   map_state: 2 ; seq 正向map;
// *****========   --->   *****======   map_state: 3 ; seq 逆向map;

pair<string,mstat> break_point_get(string& break_p1,Kmer_Chain chain,fastqReader& reference,const string& soft_clip_seq,string& match_seq,map<string,unsigned>& map_chrom_offset,vector<string>& v_chroms,int map_state){
    int error_hanming = 3;
    string new_break_p1 = break_p1;
    position chrom_pos = int2position(chain.pos,map_chrom_offset,v_chroms,0,v_chroms.size()-1);
    cout << "Get Pos From: " << chain.pos << " To: " << chrom_pos.first <<":" << chrom_pos.second << "\t"
         << match_seq << "\t" << soft_clip_seq << "\tMap Sate:" << map_state <<"\n";
    // 需要对SA的候选断点进行验证，防止某些位置多个断点，造成假的mapping 点。
    if (chain.sa) {
        int test_len = 15;
        string query_seq,refer_seq;
        if(map_state == 0 || map_state == 3){
            refer_seq = reference.fetch(chrom_pos.first,chrom_pos.second,chrom_pos.second+test_len-1);
            //query_seq = soft_clip_seq.substr(0,test_len);
        } else {
            refer_seq = reference.fetch(chrom_pos.first,chrom_pos.second-test_len+1,chrom_pos.second);
            //query_seq = soft_clip_seq.substr(soft_clip_seq.size()-test_len,test_len);
        }
        // 为了解决可能存在的融合片段间有插入序列，对soft-clip seq 进行滑窗检测（幅度 7 bp)
        int n = 0 ;// 计算分值
        int shift = 0; // 滑动偏移值 
        while (soft_clip_seq.size() >= shift + test_len ){
            if(map_state == 0 || map_state == 3) query_seq = soft_clip_seq.substr(min(shift,(int)soft_clip_seq.size()-test_len),test_len);
            else query_seq = soft_clip_seq.substr(max((int)soft_clip_seq.size()-test_len-shift,0),test_len);

            // cout << "SA refer: \t" << refer_seq << "\nquery_seq:\t" << query_seq << '\n';
            n = SW_map(refer_seq,query_seq);
            if(n >= 8) break;
            shift += 7;
        }
        // if(map_state == 0 || map_state == 1) query_seq = soft_clip_seq.substr(0,test_len);
        // else query_seq = soft_clip_seq.substr(soft_clip_seq.size()-test_len,test_len);

        //cout << "Score: " << n << "\t" << "Shift: " << shift <<'\n';
        cout << "SA refer: \t" << refer_seq << "\nquery_seq:\t" << query_seq << '\n';
        // int n = SW_map(refer_seq,query_seq);
        if (n >= 8)
            return make_pair(break_p1,make_pair(chrom_pos,map_state));
        else
            return make_pair(break_p1,make_pair(make_pair("",0),map_state)); // 不合格的SA不计入后续的融合计算
    }
    int mlen = (int)match_seq.size();
    if (map_state == 0 || map_state == 3){
        string query_seq = soft_clip_seq.substr(0,chain.offset);
        string refer_seq = reference.fetch(chrom_pos.first,chrom_pos.second-chain.offset,chrom_pos.second-1);
        reverse(refer_seq.begin(),refer_seq.end());
        reverse(query_seq.begin(),query_seq.end());
        int n = break_point_search(refer_seq,query_seq,error_hanming);
        chrom_pos.second -= n;
        // 断点附近可mapping 位置检查是否需要再移动
        if(!comparePos(break_p1,chrom_pos) && n == chain.offset){
            string refer2 = reference.fetch(chrom_pos.first,chrom_pos.second-mlen,chrom_pos.second-1);
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
                new_break_p1 = move_breakp(break_p1,-n);
            } else {
                chrom_pos.second -= n;
                new_break_p1 = move_breakp(break_p1,n);
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
            string refer_seq = reference.fetch(chrom_pos.first,chrom_pos.second+12,chrom_pos.second+12+check_len-1);
            int n = break_point_search(refer_seq,query_seq,error_hanming);
            chrom_pos.second += (n + 11);
            if (n != check_len) return make_pair(new_break_p1,make_pair(chrom_pos,map_state));
        }
        if(!comparePos(break_p1,chrom_pos)){
            string refer2 = reference.fetch(chrom_pos.first,chrom_pos.second+1,chrom_pos.second+mlen);
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
                new_break_p1 = move_breakp(break_p1,-n);
            } else {
                chrom_pos.second += n;
                new_break_p1 = move_breakp(break_p1,n);
            }
        }      
    }
    auto fusion = make_pair(new_break_p1,make_pair(chrom_pos,map_state));
    //cout << "原始断点移动为： "<< break_p1 <<"\t新的点位为： " << fusion.first << '\n';
    return fusion;
}

//对soft-clip seq 进行重新比对，找到其原本的基因组位置，并找到断裂位点（同时修正断点位置）；
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
vector<pair<string,mstat>> map_soft_clip_seq(string& break_p1,unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,
fastqReader& reference,Piled_reads& prd,map<string,unsigned>& map_chrom_offset,
vector<string>& v_chroms,map<string,vector<pair<unsigned,bool>>>& map_sa ,map<string,vector<int>>& map_breaks_pre,bool rna)
{
    string& soft_clip_seq = prd.seq;
    string& match_seq = prd.match_seq;
    bool isdown = prd.down;
    vector<pair<string,mstat>> results;
    if(soft_clip_seq.size() < 18) {
        //results.push_back(result)
        return results;
    }

    //只对满足条件的断点,进行后续的断点信息检索
    //1. ucount >= 10; 单独满足这就OK
    //2. check_num 至少2 ；
    //3. check_state 排名前三；
    //4. check_num 占总 break_points 10% 以上
    if(prd.ucount >= 10) ;
    else {
        //【断点reads总数，本融合断点数，断点reads类型总数，本断点数量排名】
        auto vinfos = break_point_info2(map_breaks_pre[break_p1],prd.ucount,0); 
        //cout <<prd.count << '\t' << vinfos[1] << '\t' << vinfos[3] << '\t' << vinfos[0] << '\n';
        //cout << (vinfos[1]+1.0)/(vinfos[0]+1.0) << '\n';
        if( prd.count >= 2 && vinfos[3] <= 3 && (vinfos[1]+1.0)/(vinfos[0]+1.0) >= 0.1) ; 
        else return results;
    }


    if(!isdown) reverse(soft_clip_seq.begin(),soft_clip_seq.end()); // soft-clip seq 在上游的read ,在pilup 时用反向的read，所以要转换回来。  
    mop map_kmer_offset_forward = extract_kmer(soft_clip_seq,12);
    vector<unsigned int> v_kmer_forward = select_kmer(map_kmer_pos,map_kmer_offset_forward,5000,isdown,rna,true);

    // 看看是选择了哪几个kmer
    //cout << "Forward: ";
    //for(auto i :v_kmer_forward) 
        //cout << map_kmer_offset_forward[i] << ":" << soft_clip_seq.substr(map_kmer_offset_forward[i],12) <<
        //"(" << map_kmer_pos[i].size() << ")" << "\t";
    //cout << '\n';
    string soft_clip_seq_rp = reverse_complete(soft_clip_seq);
    mop map_kmer_offset_reverse = extract_kmer(soft_clip_seq_rp,12);
    vector<unsigned int> v_kmer_reverse = select_kmer(map_kmer_pos,map_kmer_offset_reverse,5000,isdown,rna,false);
    //cout << "Reverse: ";
    //for(auto i :v_kmer_reverse) 
        //cout << map_kmer_offset_reverse[i] << ":" << soft_clip_seq_rp.substr(map_kmer_offset_reverse[i],12) << 
        //"(" << map_kmer_pos[i].size() << ")" << "\t";
    //cout << '\n';

    int error_golbal = 3;

    //获得可能的supplementary 位点,半径10bp范围；
    vector<pair<unsigned,bool>> v_sa1,v_sa;
    vector<string> v_break_p1 = split(break_p1,":");
    int break_p1_pos = string2number<int>(v_break_p1[1]);
    string sa_pos1 = v_break_p1[0] + ":" + to_string(break_p1_pos/10 - 1);
    string sa_pos2 = v_break_p1[0] + ":" + to_string(break_p1_pos/10 );
    string sa_pos3 = v_break_p1[0] + ":" + to_string(break_p1_pos/10 + 1);
    vector<string> sa_poses = {sa_pos1,sa_pos2,sa_pos3};
    //cout << break_p1 << '\t' << sa_pos1 << '\t' << sa_pos2 << '\t' << sa_pos3 << "\t";
    for(auto pos : sa_poses){
        if (map_sa.find(pos) != map_sa.end()){
            for(auto i:map_sa[pos]){
              v_sa1.push_back(i);  
              //cout << i << '\t';
            } 
        }
    }
    //cout << '\n';
    sort_unique<pair<unsigned,bool>>(v_sa1);
    unsigned sa_tmp = 0;
    for(auto p_sa:v_sa1){
        if(p_sa.first + 2 > sa_tmp && sa_tmp + 2 > p_sa.first) ;
        else{
            v_sa.push_back(p_sa);
            sa_tmp = p_sa.first;
        }
    }

    // 这些参数用来方便后续断点的移动检索；
    // vector<Kmer_Chain> *select_chains;
    vector<Kmer_Chain> tmp_chains;
    // string *select_soft_seq;
    int map_state = -1;

    vector<Kmer_Chain> pos_global_forwards;
    vector<Kmer_Chain> pos_global_reverses;
    if(isdown){
        // 下游正链匹配
        pos_global_forwards = global_align(map_kmer_pos,map_kmer_offset_forward,v_kmer_forward,error_golbal,v_sa,prd.count,prd.ucount,0) ;
        auto& pos_global_forward = pos_global_forwards[0];
        // 下游负链匹配
        pos_global_reverses = global_align(map_kmer_pos,map_kmer_offset_reverse,v_kmer_reverse,error_golbal,v_sa,prd.count,prd.ucount,1,pos_global_forward.chainNum);
        auto& pos_global_reverse = pos_global_reverses[0];

        map_state = 10; 
        // cout << "down! \n";
        // 保留SA的mapping ,防止因为snp 导致的，遗漏；
        if(pos_global_forward.chainNum > pos_global_reverse.chainNum && pos_global_forward.chainNum >=2 ){   
            for(auto i=pos_global_reverses.size()-1;;i--){
                if (pos_global_reverses[i].sa) tmp_chains.push_back(pos_global_reverses[i]);
                else break;
                if (i==0) break;
            }
            pos_global_reverses = tmp_chains;       
        } 
        else if(pos_global_forward.chainNum < pos_global_reverse.chainNum && pos_global_reverse.chainNum >=2)
        {
            for(auto i=pos_global_forwards.size()-1;;i--){
                if (pos_global_forwards[i].sa) tmp_chains.push_back(pos_global_forwards[i]);
                else break;
                if (i==0) break;
            }
            pos_global_forwards = tmp_chains;
                                 
        } 
        else if(pos_global_forward.chainNum < 2 ){
            map_state = -1;
        }

        // cout << "Mapping Finsh!\n";     
    } else {
        pos_global_forwards = global_align(map_kmer_pos,map_kmer_offset_forward,v_kmer_forward,error_golbal,v_sa,prd.count,prd.ucount,2);
        auto& pos_global_forward = pos_global_forwards[0];
        pos_global_reverses = global_align(map_kmer_pos,map_kmer_offset_reverse,v_kmer_reverse,error_golbal,v_sa,prd.count,prd.ucount,3,pos_global_forward.chainNum);
        auto& pos_global_reverse = pos_global_reverses[0];
        
        map_state = 32; 
        // cout << "Mapping Finsh!\n";
        if(pos_global_forward.chainNum > pos_global_reverse.chainNum && pos_global_forward.chainNum >=2){
            for(auto i=pos_global_reverses.size()-1;;i--){
                if (pos_global_reverses[i].sa) tmp_chains.push_back(pos_global_reverses[i]);
                else break;
                if (i==0) break;
            }
            pos_global_reverses = tmp_chains;         
        }
        else if(pos_global_forward.chainNum < pos_global_reverse.chainNum && pos_global_reverse.chainNum >= 2){
            for(auto i=pos_global_forwards.size()-1;;i--){
                if (pos_global_forwards[i].sa) tmp_chains.push_back(pos_global_forwards[i]);
                else break;
                if (i==0) break;
            }
            pos_global_forwards = tmp_chains;
        }
        else if(pos_global_reverse.chainNum < 2){
            map_state = -1;
        }
    }

    // 设定断点滑动检索使用的参数；
    if(map_state == -1) return results;
    // if (map_state == 0 || map_state == 2){
    //     select_chains = &pos_global_forwards;
    //     select_soft_seq = &soft_clip_seq;
    // } else {
    //     select_chains = &pos_global_reverses;
    //     select_soft_seq = &soft_clip_seq_rp;
    // }

    // 滑动检索断点信息；
    if(map_state == 10){
        for(auto sc:pos_global_forwards){
            auto result = break_point_get(break_p1,sc,reference,soft_clip_seq,match_seq,map_chrom_offset,v_chroms,0);
            results.push_back(result);
        }
        for(auto sc:pos_global_reverses){
            auto result = break_point_get(break_p1,sc,reference,soft_clip_seq_rp,match_seq,map_chrom_offset,v_chroms,1);
            results.push_back(result);
        }
    }
    else if(map_state == 32){
        for(auto sc:pos_global_forwards){
            auto result = break_point_get(break_p1,sc,reference,soft_clip_seq,match_seq,map_chrom_offset,v_chroms,2);
            results.push_back(result);
        }
        for(auto sc:pos_global_reverses){
            auto result = break_point_get(break_p1,sc,reference,soft_clip_seq_rp,match_seq,map_chrom_offset,v_chroms,3);
            results.push_back(result);
        }
    }
    
    return results;   
}




//处理 bam1_t *aln 对象（bam 比对结果对象），找到 split-reads 和 discordant reads
//处理含soft-clip 的reads 或 insert size 大于 1000 的reads , mapQ 大于等于 15
//该function 会修改 map_dcp 和 map_split_read 这两个表记录了符合对应条件的reads 信息。
//记录每条序列的气势位点和insertSize;
//得到supplementary alignment 位置（ 以防多位点匹配时遗漏，SA:Z:7,55863708,-,78M73S,255,0;）
//@aln bam 的比对记录；
//@bam_header bam的头文件；
//@map_dcp: 散列表用于记录discordant reads;
//@map_split_read: 散列表用于记录split reads;
//@map_sa: 记录first alignment --> supplementary alignment的映射
bool deal_aln(const bam1_t* aln,const bam_hdr_t* bam_header,unordered_map<string,fusionPos>& map_dcp,unordered_map<string,vector<Piled_reads>>& map_split_read,unordered_map<string,vector<string>>& map_alt_split,unordered_map<string,string>& map_transcript,
map<string,unsigned>& map_chrom_offset,map<string,vector<pair<unsigned,bool>>>& map_sa){
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
    pair<unsigned,string> p_sa = make_pair(0,"");
    string sa = getAux(aln,"SA");
    if(sa!="") {
        vector<string> v_sa = split(sa,",");
        vector<string> v_sa_1 = split(v_sa[0],":");
        unsigned sa_pos = position2int(v_sa_1[2],string2number<int>(v_sa[1]),map_chrom_offset);
        p_sa = make_pair(sa_pos,v_sa[3]);
        //cout << "Get SA: " << pos << "\t" << sa << "\t" <<v_sa_1[2] <<":"<< v_sa[1] << "\n";
    }
    bool ptag = true;
    if(flag & 256 || flag & 2048) ptag = false;
    // cout << qname << "\t" << cigar << "\t" << chrom << ":" << pos <<"\t" << nchrom <<":" << npos << "\n";

    if(s_chroms.find(chrom) == s_chroms.end()) return false;
    if(mapQ < 0 || flag & 1024 )  return true;
    else if(parse_split_read(chrom,pos,seq,cigar,qname,map_split_read, map_alt_split,map_transcript,p_sa,map_sa,ptag)) {
        return true;
    }
    else if(abs(pos-npos) > 10000 || ( chrom != nchrom && nchrom != "NA")) {
        parse_discordant_read(qname,flag,chrom,pos,nchrom,npos,cigar,map_dcp);
    }

    return true;
}

//对split-read soft clip reads 进行比对找到符合的断点；
//pair<position,int> result = map_soft_clip_seq(map_kmer_pos,reference,read,map_chrom_offset,v_chroms,true)
//void combine_split_reads(string& break_p1,pair<string,int>& break_p2,int num,map<string,fusion>& map_fusion, map<string,int>& map_c2i);
//mapping soft-clip seq 到新的断裂点，并将信息记录到散列表 map_fusion 中；
//@rna: rna 的比对和DNA 不一致
void map_split_reads(string& break_p1,Piled_reads& prd,unordered_map<unsigned int,vector<unsigned>>& map_kmer_pos,fastqReader& reference,map<string,unsigned>& map_chrom_offset,vector<string> v_chroms,map<string,fusion>& map_fusion,
map<string,vector<pair<unsigned,bool>>>& map_sa ,map<string,vector<int>>& map_breaks_pre,bool rna){
    vector<pair<string,mstat>> results = map_soft_clip_seq(break_p1,map_kmer_pos,reference,prd,map_chrom_offset,v_chroms,map_sa,map_breaks_pre,rna);
    if(results.size() == 0) return;
    for (auto res: results){
        auto& result = res.second;
        if (result.first.second == 0) continue;
        cout << "first Map " << res.first << "\t" << result.first.first << ":" << result.first.second << "\n";
        combine_split_reads(res.first,result,prd,map_fusion,map_c2i);  
    }
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
            seq1 = reference.fetch(fu.p1.first,fu.p1.second,fu.p1.second + len -1);
            seq2 = reference.fetch(fu.p2.first,fu.p2.second+1,fu.p2.second + len);
            len = cmpSeq(seq1,seq2);
            fu.p1.second += len;
            fu.p2.second += len;
        
            } else if(p1_extron.second > 0){
                // len = p1_extron.second;
                seq1 = reference.fetch(fu.p1.first,fu.p1.second - len + 1,fu.p1.second);
                seq2 = reference.fetch(fu.p2.first,fu.p2.second - len,fu.p2.second - 1);
                len = cmpSeq(seq1,seq2);
                fu.p1.second -= len;
                fu.p2.second -= len;
            }
        } else if(!op2_extron.second) {
            if (p2_extron.second < 0){
                // len = -p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second,fu.p2.second + len - 1);
                seq2 = reference.fetch(fu.p1.first,fu.p1.second+1,fu.p1.second + len);
                len = cmpSeq(seq1,seq2);
                fu.p1.second += len;
                fu.p2.second += len;
            
            } else if(p2_extron.second > 0){
                // len = p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second - len + 1,fu.p2.second);
                seq2 = reference.fetch(fu.p1.first,fu.p1.second - len,fu.p1.second - 1);
                len = cmpSeq(seq1,seq2);
                fu.p1.second -= len;
                fu.p2.second -= len;
                
            }

        }
        
    } else {
        if(!op1_extron.second){
            if(p1_extron.second < 0){
                // len = -p1_extron.second;
                seq1 = reference.fetch(fu.p1.first,fu.p1.second,fu.p1.second + len -1);
                seq2 = reference.fetch(fu.p2.first,fu.p2.second - len,fu.p2.second - 1);
                seq1 = reverse_complete(seq1);
                len = cmpSeq(seq1,seq2);
                fu.p1.second += len;
                fu.p2.second -= len;
                
            }else if(p1_extron.second > 0){
                // len = p1_extron.second;
                seq1 = reference.fetch(fu.p1.first,fu.p1.second - len + 1,fu.p1.second);
                seq2 = reference.fetch(fu.p2.first,fu.p2.second+1,fu.p2.second + len);
                seq1 = reverse_complete(seq1);
                len = cmpSeq(seq1,seq2);
                fu.p1.second -= len;
                fu.p2.second += len;
        
            }
        } else if (!op2_extron.second){
            if(p2_extron.second < 0){
                // len = -p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second,fu.p2.second + len - 1);
                seq2 = reference.fetch(fu.p1.first,fu.p1.second - len,fu.p1.second - 1);
                seq1 = reverse_complete(seq1);
                len = cmpSeq(seq1,seq2);
                fu.p2.second += len;
                fu.p1.second -= len;
                
            } else if(p2_extron.second > 0){
                // len = p2_extron.second;
                seq1 = reference.fetch(fu.p2.first,fu.p2.second - len + 1,fu.p2.second);
                seq2 = reference.fetch(fu.p1.first,fu.p1.second+1,fu.p1.second + len);
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
void combine_rna_fus(vector<pair<fusion,bool>>& vmf,map<pair<string,int>,vector<int>>& map_breaks,unordered_map<string,int>& map_spreads_count,int dis){
    if(vmf.size() == 0) return;
    for(unsigned i = 0; i < vmf.size();i++){
        if(vmf[i].second) continue;
        // cout << "Combine: " << vmf[i].first.p1.first << ":" << vmf[i].first.p1.second << "-" << vmf[i].first.p2.first << ":"
        //     <<vmf[i].first.p2.second << "\t";
        int s_max = vmf[i].first.p1n + vmf[i].first.p2n;
        auto p1 = vmf[i].first.p1;
        auto p2 = vmf[i].first.p2;
        int p1sn = 0; //用于记录合并的不同点位的p1 点位的次数
        int p2sn = 0; //用于记录合并的不同点位的p2 点位的次数
        for(unsigned j = i+ 1;j < vmf.size();j++){
            if(vmf[i].first.p2.first == vmf[j].first.p2.first && abs(vmf[i].first.p2.second - vmf[j].first.p2.second) <= dis && abs(vmf[i].first.p1.second - vmf[j].first.p1.second) <= dis){
                if(vmf[i].first.p2.second - vmf[j].first.p2.second != 0) p2sn++;
                if(vmf[i].first.p1.second - vmf[j].first.p1.second !=0 ) p1sn++;

                vmf[i].first.p1n += vmf[j].first.p1n;
                vmf[i].first.p2n += vmf[j].first.p2n;
                vmf[i].first.dcp += vmf[j].first.dcp;
                for(auto s:vmf[j].first.vp1) vmf[i].first.vp1.push_back(s);
                for(auto s:vmf[j].first.vp2) vmf[i].first.vp2.push_back(s);
                for(auto s:vmf[j].first.p1_cigars) vmf[i].first.p1_cigars.push_back(s);
                for(auto s:vmf[j].first.p2_cigars) vmf[i].first.p2_cigars.push_back(s);
                for(auto s:vmf[j].first.vdcp) vmf[i].first.vdcp.push_back(s);
                if(vmf[i].first.split_seqs.first == "") vmf[i].first.split_seqs.first = vmf[j].first.split_seqs.first;
                if(vmf[i].first.split_seqs.second == "") vmf[i].first.split_seqs.second = vmf[j].first.split_seqs.second;
                vmf[j].second = true;
                // cout << vmf[j].first.p1.first << ":" << vmf[j].first.p1.second << "-" << vmf[j].first.p2.first << ":"
                //     <<vmf[j].first.p2.second << "\t";
                int s_tmp = vmf[j].first.p1n + vmf[j].first.p2n;
                if(s_tmp > s_max) {
                    p1 = vmf[j].first.p1; //记录是附近最多split-reads 支持的点位；
                    p2 = vmf[j].first.p2;
                    s_max = s_tmp;
                }
            }
        }
        // cout << "\n";
    
        auto& f = vmf[i].first;
        f.p1 = p1;
        f.p2 = p2;
        cout << f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" << f.p1n << ":"<<f.p2n<<":" << f.dcp  << "\n";
        auto vp_counts1 = break_point_info(p1,map_breaks,map_spreads_count);
        auto vb_infos_p1 = break_point_info2(vp_counts1,f.p1n,p1sn);
        auto vp_counts2 = break_point_info(p2,map_breaks,map_spreads_count);
        auto vb_infos_p2 = break_point_info2(vp_counts2,f.p2n,p2sn);

        // cout<< f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" << f.p1n << ":"<<f.p2n<<":" << f.dcp  << "\n";
        if (min(f.p1_cigars.size(),f.p2_cigars.size()) >= 1 || f.dcp >= 2)
        {
            if( true ){
                //sort_unique<string>(f.vp);
                sort_unique<string>(f.vp1);
                sort_unique<string>(f.vp2);
                sort_unique<string>(f.vdcp);

                int p1c = templateNum(f.p1_cigars);
                int p2c = templateNum(f.p2_cigars);
                cerr << f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" 
                    << f.mpdr.first << ";" << f.mpdr.second << "\t" 
                    << f.genedr.first << ";" << f.genedr.second << "\t"
                    << f.rpdr.first << ';' << f.rpdr.second << "\t"
                    << f.genes.first << ";" << f.genes.second << "\t"
                    << f.p1n << ":"<<f.p2n<<":" << f.vdcp.size()  <<  "\t" 
                    << p1c << ":" << p2c << "\t"
                    << vb_infos_p1[0] << ':' << vb_infos_p1[1] << ':' << vb_infos_p1[2] << ':' << vb_infos_p1[3] << '\t'
                    << vb_infos_p2[0] << ':' << vb_infos_p2[1] << ':' << vb_infos_p2[2] << ':' << vb_infos_p2[3] << '\t'
                    << f.split_seqs.first << ":" << f.split_seqs.second << "\t";
                for (auto& i :f.vp1) cerr << i << ";";
                cerr << "\t";
                for(auto& i :f.vp2) cerr << i << ";";
                cerr << "\t";
                for(auto& i:f.vdcp) cerr << i << ";";
                cerr << "\n";
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
void rna_print_sv(map<string,fusion>& map_fusion,map<pair<string,int>,vector<int>>& map_breaks,unordered_map<string,int>& map_spreads_count,int dis=10){
    vector<pair<string,fusion>> vfs(map_fusion.begin(),map_fusion.end()); 
    if(vfs.size() == 0) return;
    sort(vfs.begin(),vfs.end(),cmpFu);
    vector<pair<fusion,bool>> vmf;
    fusion sf = vfs.begin()->second;
    for(auto& mf : vfs){
        auto &f = mf.second;
        cout<< "未合并前： " << mf.first << "\t" << f.p1.first <<":" <<f.p1.second << "\t" << f.p2.first<<":"<<f.p2.second << "\t" << f.mpdr.first << "\t" << f.mpdr.second << "\t"  << f.p1n << ":"<<f.p2n<<":" << f.dcp  << "\n";
        // 对于只有单边支持（支持度高），以及 pair_end  支持信息的融合，回顾性检测另一断点支持度信息；

        if(mf.second.p1.first == sf.p1.first && abs(mf.second.p1.second - sf.p1.second) <= dis){
            vmf.push_back(make_pair(mf.second,false));
            sf = mf.second;
        } else {
            combine_rna_fus(vmf,map_breaks,map_spreads_count,dis);
            sf = mf.second;
            vmf = {make_pair(sf,false)};
        }
    }
    combine_rna_fus(vmf,map_breaks,map_spreads_count,dis);
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
        cout << argv[0] << "\tbam\treference\t12-kmer-index\trefseq_trans.exon.bed\ttranscript-altsplit\t0:DNA(1:RNA)\tChr.start\tChr.end\n";
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
    map<string,vector<pair<unsigned,bool>>> map_sa; // 记录supplementary aignment的mapping位置,其中bool 记录soft-clip 的位置，false 为MS， true 为 SM.
    map<pair<string,int>,vector<int>> map_breaks; // 记录断点处的split reads 数量及种类（记录断点复杂度，primary Alignment）
    //获得外显子、基因信息
    map<string,vector<string>> m_extron = getExtronInfo(extronFile,mge);
    // 获得要进行可变剪切检测的相关转录本信息
    unordered_map<string,string> map_transcript = getTranscript(transFile);
    //记录split reads 的数量， key 是最长soft clip seq ，value 是符合的数量
    unordered_map<string,int> map_spreads_count;


    cout << "开始获取融合支持reads ...\n";
    // 寻找要处理的reads 
    if(argc == 7){
        while(brd.next()){
            if (!deal_aln(brd.aln,brd.bam_header,map_dcp,map_split_read,map_alt_split,map_transcript,map_chrom_offset,map_sa)) break;       
        }
    } else if(argc == 9){
        vector<string> v_pos1 = split(argv[7],":");
        vector<string> v_pos2 = split(argv[8],":");
        string regin1 = v_pos1[0] + ":" + to_string(string2number<int>(v_pos1[1])-10000) + "-" + to_string(string2number<int>(v_pos1[1])+10000);
        string regin2 = v_pos2[0] + ":" + to_string(string2number<int>(v_pos2[1])-10000) + "-" + to_string(string2number<int>(v_pos2[1])+10000);

        Alignment align1 {brd.sam_file,brd.bam_index,brd.bam_header,regin1};
        while(align1.next()){
            if (!deal_aln(align1.aln,brd.bam_header,map_dcp,map_split_read,map_alt_split,map_transcript,map_chrom_offset,map_sa)) break;       
        }

        Alignment align2 {brd.sam_file,brd.bam_index,brd.bam_header,regin2};
        while(align2.next()){
            if (!deal_aln(align2.aln,brd.bam_header,map_dcp,map_split_read,map_alt_split,map_transcript,map_chrom_offset,map_sa)) break;       
        }
    } else {
        cout << "input args numbers wrong(7 or 9)";
        return -1;
    }
    cout << "融合支持reads获取结束！\n";

    ifstream kmer_index {kmerFile};
    unordered_map<unsigned int,vector<unsigned>> map_kmer_pos;
    boost::archive::binary_iarchive iarch(kmer_index);
    iarch >> map_kmer_pos;
    kmer_index.close();

    // map_sa 的值uniq
    for(auto& sa:map_sa) sort_unique<pair<unsigned,bool>>(sa.second);

    // 统计断点的花reads信息
    map<string,vector<int>> map_breaks_pre;
    for (auto&sr :map_split_read){
        string break_p1 = sr.first;
        for(auto& fu:sr.second){
            map_breaks_pre[break_p1].push_back(fu.ucount);
        }
    }

    //mapping split-reads 把相同断点的 split-reads 整合到一起。
    cout << "把相同断点的 split-reads 整合到一起。\n";
    for(auto& sr:map_split_read){
        string break_p1 = sr.first;
        //cout << "Piled reads: " << break_p1 << "\n";
 
        for(auto& fu:sr.second){
            string break_p1 = sr.first;
            //cout << "FF: " << break_p1 << "\t" <<fu.seq << "\t" << fu.count << "\t" << fu.down << "\n";
            map_split_reads(break_p1,fu,map_kmer_pos,reference,map_chrom_offset,v_chroms,map_fusion,map_sa,map_breaks_pre,rna);
        }
    }

    //将discordant reads 的支持信息也增加到split-reads 信息中
    cout << "将discordant reads 的支持信息也增加到split-reads 信息中\n";
    unordered_map<string,vector<string>> map_fus_dcps = combine_discordant_reads(map_dcp);
    combine_fusion(map_fusion,map_fus_dcps,map_dcp);

    //修正基因上下游关系 并 计算断点的复杂度
    for(auto& mfu:map_fusion){
        auto& fu = mfu.second;
        md_rna_breakpos(fu,m_extron,reference,mge);
        if(fu.p1n != 0) map_breaks[fu.p1].push_back(fu.p1n);
        if(fu.p2n !=0) map_breaks[fu.p2].push_back(fu.p2n);
        // map_spreads_count[fu.split_seqs.first] = fu.p1n;
        // map_spreads_count[fu.split_seqs.second] = fu.p2n;
    }

    //展示查找到的融合信息
    cerr << "break_pos1\tbreak_pos2\tsoft_map_dr(p1);soft_map_dr(p2)\tgene_p1_dr;gene_p2_dr\tsupport_p1_dr;support_p2_dr\t"
        << "gene_p1;gene_p2\tsplit_p1:split_p2:discordant reads\tp1_cigars;p2_cigars;total_reads\n"; 
    if(false) dna_print_sv(map_fusion);
    else rna_print_sv(map_fusion,map_breaks,map_spreads_count);

    cerr << "\nAlterNative SplitInfos\n";
    for(auto& a : map_alt_split){
        cerr << a.first << "\t" << a.second.size() <<"\n";
    }

}


