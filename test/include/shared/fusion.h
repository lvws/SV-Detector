#include<string>
#include<vector>
#include<utility>
#include <algorithm>
#include <math.h>
#include <map>
#include <unordered_map>
#include <regex>
#include <sstream>

using namespace std;

//pileup reads struct
//记录堆叠的split reads,最长序列和堆叠次数。
//@seq : 记录最长的soft-clip seq;
//@match_seq: read soft-clip seq 旁边和参考基因组一至的序列（默认取10bp)
//@count: 记录堆叠的次数；
//@qnames: 记录堆叠的read的名称；
//@cigars: 记录堆叠的cigar 名称
//@down： true 表示soft-clip seq 处于下游，反之～。
struct Piled_reads
{
    string seq;
    string match_seq;
    int count=0;
    vector<string> qnames;
    vector<string> cigars;
    bool down = false;
};

using mps = unordered_map<string,vector<Piled_reads>>; // 记录break point 对应的 soft reads 堆叠，和部分 mapping seq.

//map chrom to int
map<string,int> map_c2i = {{"1",1},{"2",2},{"3",3},{"4",4},{"5",5},{"6",6},{"7",7},{"8",8},
    {"9",9},{"10",10},{"11",11},{"12",12},{"13",13},{"14",14},{"15",15},{"16",16},{"17",17},{"18",18},
    {"19",19},{"20",20},{"21",21},{"22",22},{"X",23},{"Y",24}};


template <class T>
T string2number(const string& s);

//汉明距离比较两个序列，
//在预定错误阈值下，通过比对返回较长的序列，否则是 空字符串；
//默认堆叠read 允许错误为2
string hanming(const string& s1,const string& s2,int td=2);

//堆叠新加入的soft-clip seq
//堆叠使用的方法是，汉明距离，容错默认是3个碱基
//保留最长的seq.
//@v_ss: 记录着该位点堆叠的序列和数量；
//@s: 需要进行比较的序列；
//@m: reads 上和 参考基因组match 的序列（10bp）
//@qname: reads 的 name;
//@down: soft-clip是否在断裂点下游;
//@p_sa: 记录supplementary alignment 的绝对位置和 cigar 
//@map_sa: 记录first alignment --> supplementary alignment的映射
//return 更新堆叠序列集合。
// void pileup(vector<Piled_reads>& v_ss,string& s,string& m,string qname,bool down);

//处理split-read, 对一断点的soft-clip seq 进行堆叠；
//该function 会修改 map_split_read 的内容，
//map_split_read : key 为 chrom<string>:pos<int> ; value 为 Piled_reads 向量。
bool parse_split_read(string& chrom,int pos,string& seq,string& cigar,string& qname,unordered_map<string,vector<Piled_reads>>& map_split_read,unordered_map<string,vector<string>>& map_alt_split,unordered_map<string,string>& map_transcript,
pair<unsigned,string>& p_sa,map<string,vector<unsigned>>& map_sa);

//融合断点结构,记录两个断点的位置。
//自动会把较低的染色体位点记录为p1, 较高的位点记录为p2；
//记录每个点位左右偏移cigar记录距离的点位，便于后续和split reads 支持的融合进行合并；
// p1pos 记录低位置，p2pos 记录高位置
struct fusionPos{
    pair<string,int> p1;
    pair<string,int> p2;
    pair<int,int> p1pos;
    pair<int,int> p2pos;
};

//discordant read 记录到字典中（序列id -> 比对位置）
//方便后续使用该字典计算支持reads
void parse_discordant_read(string& qname,int flag,string& rname,int pos,string& nchrom,int npos,string& cigar,unordered_map<string,fusionPos>& map_dcp);

//将融合点位进行整合方便后续的检索比较。
//pos 除去 1000 记录这个区间的reads。
//减少后续的计算支持该断点融合的，discordant read 的数量。
//@map_dcp:  以discordant read 的 qname 为键，融合断点为值的字典；
//return map_fus_dcps: 汇总1000bp 范围内的discordant reads ,简化后续计算
unordered_map<string,vector<string>> combine_discordant_reads(unordered_map<string,fusionPos>& map_dcp);




//记录read 的mapping 情况
//@slen : 记录soft-clip seq 长度；
//@offset: 记录偏离长度（主要用于soft-clip 在下游的情况）；
//@alsp: 记录可变剪切的偏移位点；
struct altSplit
{
    int slen = 0 ,offset = 0;
    vector<pair<int,int>> altsp;
};



//记录融合信息的结构
//@p1: 记录上游断点；
//@p2: 记录下游断点；
//@check_points 用于在discordant reads 集合中检索合适的支持(按染色体位置大小顺序构造融合名称)；
//@p1n: 位于上游断点的split reads 数量；
//@p2n: 位于下游断点的split reads 数量;
//@rpdr: 记录match seq  read 处与上游（‘u’）,下游('d');
//@mpdr: 记录上下游Soft-Clip read 的 map 方向 ‘+’ 正向 ‘-’ 反向；
//@vp: 断点的split reads name 集合；
//@dcp: 支持断点的discordant reads数。
//@genes: 记录融合的上下游基因。
//@genedr: 记录基因的转录方向, ‘+’ 正向 ‘-’ 反向, '0' 未找到基因。
//@split_seqs: 记录soft-clip seqs。
//@p1_cigars: 支持p1 断点的reads 的cigar 集合；
//@p2_cigars: 支持p2 断点的reads 的cigar 集合；
struct fusion
{
    pair<string,int> p1;
    pair<string,int> p2;
    pair<char,char> rpdr = make_pair('0','0');
    pair<char,char> mpdr;
    pair<char,char> genedr = make_pair('0','0');
    pair<string,string> split_seqs;
    string check_point;
    int p1n = 0;
    int p2n = 0;
    vector<string> vp;
    vector<string> p1_cigars;
    vector<string> p2_cigars;
    int dcp = 0;
    pair<string,string> genes;
};

//整合split read 的break points
//使用 断点 上游-下游 作为Key 记录 fusion value 。
//对每一个split read，进行处理。融合结果更新到map_fusion 表中。
//@break_p1 : 以 chrom:pos 记录的断点（检索bam 时得到的第一个断裂点）；
//@break_p2: 以<string(chrom),int(pos)> 记录的断点，是soft-clip reads 重新比对到的位置；
//@num: 记录了这个断点-soft-clip reads 支持的数量；
//@map_fusion:  以break_point1-break_point2 为 key, fusion 为值的 散列表；
//@map_c2i: 染色体号和数字对应表，用于比较断裂点先后；
void combine_split_reads(string& break_p1,pair<pair<string,int>,int>& map_p2,Piled_reads& prd,map<string,fusion>& map_fusion, map<string,int>& map_c2i);

//根据split read 断点信息，寻找支持的 discordant reads;
//split read 的 fusion 信息记录在 散列表 map_fusion 中；
//储备设定在断裂点附近。目前粗略计算在他们附近100bp,都作为支持项；
void combine_fusion(map<string,fusion>& map_fusion,unordered_map<string,vector<string>> map_fus_dcps,unordered_map<string,fusionPos>& map_dcp);
