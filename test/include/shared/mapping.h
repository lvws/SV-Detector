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
//#include "shared/reference.h"

using namespace std;
using position = pair<string,int>; // 基因组位点的表示形式
using mop = map<unsigned int,int>; // 记录read每个kmer对应的reads上的偏离值

// 方便将基因组绝对位置转化为 chrom:pos
vector<string> v_chroms = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
                                "19","20","21","22","X","Y"};

//构建染色体号偏离值对应表。用于以绝对值记录染色体位置。
//@fasta: hg19 基因组文件地址
//return 染色体号与偏离值对应表
map<string,int> make_chrom_offset(const string& fasta);

//二分法寻找染色体绝对位置对应的实际染色体号和实际位置
//@pos: hg19文件碱基绝对位置；
//@map_chrom_offset: 从hg19 索引文件中得到的，各个染色体的初始偏离值；
//@v_chroms： 按染色体号排序的vector;
//@s: 待检索的染色体vector 的起始位置；
//@e: 待检索的染色体vector 的终止位置。
position int2position(int pos,map<string,int>& map_chrom_offset,vector<string>& v_chroms,unsigned s,unsigned e);

// 将碱基记录为两个位的数字，seed 大小为12则 12*2 = 24位，对应16进制值0x00FFFFFF 。
unsigned int trans_kmer(const string& seq);

//切分reads 以kmer 长度，最后的碱基不够可以和倒数第二个序列overlap.
//@read : 待查询的序列
//@kmer_len: 种子长度，应该设定为12
//return: kmer(序列转换成数字)和偏离值的映射表
mop extract_kmer(const string& read,int kmer_len);

//二分法，检查给定的 位点<chrom,pos>， 与一个 位点集合（点位是按顺序存储的） 是否间隔给定的距离dis，允许误差x.
//检索到返回 true ;否则 false.
//@pos： 待检测的位点
//@v_pos ：目标位点集合
//@s: 目标位点集合开始检索位置
//@e: 目标位点集合终止检索位置
//@dis: 待检测位点pos 和目标位点集合中目标位点的距离
//@x: 待检测位点pos 和目标位点集合中目标位点的距离,可以允许的偏差（中间kmer 可能有indel）
//return : true 找到了连锁位点； false 没有找到连锁位点。
bool distance_check(int pos,const vector<int>& v_pos,int s,int e ,int dis,int x);

//检索kmer在index中的位置 vector, 过滤超多位点的kmer。将剩余的kmer按位点出现次数，由少到多排序。
//@map_kmer_pos: 事先构建好的index
//@map_kmer_offset: read中kemr 对应的偏离值；
//@threshold: 超过固定数量位点的kmer 就不进行后续的检测；
//return : kmer 集合，按在基因组存在的位点数量由少到多排序
vector<unsigned int> cheap_kmer(unordered_map<unsigned int,vector<int>>& map_kmer_pos,mop& map_kmer_offset,int threshold);

//检索kmer在index中的位置 vector, 过滤超多位点的kmer。将剩余的kmer按顺序由左至右排序。
//@map_kmer_pos: 事先构建好的index
//@map_kmer_offset: read中kemr 对应的偏离值；
//@threshold: 超过固定数量位点的kmer 就不进行后续的检测；
//return : kmer 集合，按read的 offset 排序，由低到高。
vector<unsigned int> select_kmer(unordered_map<unsigned int,vector<int>>& map_kmer_pos,mop& map_kmer_offset,int threshold);

//数据结构，记录kmer连锁的信息
struct Kmer_Chain
{
    int pos; //基因组位置
    int offset; //接近断裂点kmer 距离read
    int chainNum;   //连锁的kmer 数量
};


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
Kmer_Chain global_align(unordered_map<unsigned int,vector<int>>& map_kmer_pos,mop& map_kmer_offset,vector<unsigned int> v_kmer,int err);

//利用汉明距离算法确定基因的断裂位点；
//若融合基因在下游，序列向左移动匹配找断裂点；(提供反向序列)
//若融合基因在上游，序列向右移动匹配寻找断裂点；（提供正常序列）
//滑动配可以允许最多e个错配点；
//@refer_seq: 基因组截取的断裂位点到融合基因未能完全匹配的序列；
//@query_seq: read 中截取的靠近断裂位点未能和基因组完全匹配的序列；
//@e: 允许的最大错配数量（不允许存在indel）;
//return: 待检测的两段序列可以匹配上的碱基数量；
int break_point_search(const string& refer_seq,const string& query_seq,int e);

//反向互补序列
string reverse_complete(const string& inseq);


//根据断裂点在融合基因的5‘ 端还是 3’ 端，最终调整断裂点的位置
//@map_kmer_pos: 事先构建好的index；
//@reference: 用于获取基因组序列；
//@soft_clip_seq: 待检索的序列；
//@chrom_map: 染色体号对应表
//@up: 代表断裂点在基因的上游(基因的5‘端）；
position break_point_get(Kmer_Chain chain,fastqReader& reference,const string& soft_clip_seq,map<string,int>& map_chrom_offset,vector<string>& v_chroms,bool up);

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
pair<position,int> map_soft_clip_seq(unordered_map<unsigned int,vector<int>>& map_kmer_pos,fastqReader& reference,const string& soft_clip_seq,map<string,int>& map_chrom_offset,vector<string>& v_chroms,bool isdown);