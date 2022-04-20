#include<string>
#include<vector>
#include<utility>
#include <algorithm>
#include <math.h>
#include <map>
#include <unordered_map>
#include <regex>
#include <sstream>
#include "shared/reference.h"
#include "shared/fusion.h"

using namespace std;

//string to number
template <class T>
T string2number(const string& s){
    stringstream ss{s};
    T n;
    ss >> n;
    return n;
}

//汉明距离比较两个序列，
//在预定错误阈值下，通过比对返回较长的序列，否则是 空字符串；
string hanming(const string& s1,const string& s2,int td){
    int len = s1.size() < s2.size()? s1.size(): s2.size();
    len = min(len,24);
    int err = 0;
    for (int i = 0;i < len;i++){
        if(s1[i] != s2[i]) {
            err++;
            if (err > td) return "";
        }
    }
    return s1.size() > s2.size()? s1 : s2;
    // return "";
}

//堆叠新加入的soft-clip seq
//堆叠使用的方法是，汉明距离，容错默认是3个碱基
//保留最长的seq.
//@v_ss: 记录着该位点堆叠的序列和数量；
//@s: 需要进行比较的序列；
//@m: reads 上和 参考基因组match 的序列（10bp）
//@qname: reads 的 name;
//@down: soft-clip是否在断裂点下游;
//return 更新堆叠序列集合。
void pileup(vector<Piled_reads>& v_ss,string& s,string& m,string& qname,string& cigar,bool down){
    bool flag = false;
    for(auto& p :v_ss){
        if (p.down != down) continue;        
        string hs = hanming(p.seq,s);  
        if(hs != "") {
            p.seq = hs;
            p.count ++;
            p.qnames.push_back(qname);
            p.cigars.push_back(cigar);
            flag = true ;
        }
    }
    if (flag) return;
    v_ss.push_back({s,m,1,{qname},{cigar},down});
    return ;
}

//寻找cigar 中 S 的 长度，以及， soft-clip 在下游时，断点距离比对起点的偏移量(MND)。
//找到reads 中的可变剪切的位点
altSplit getCigarInfo(string& cigar){
    altSplit asp;
    int a = 0;
    char c ;

    stringstream ss{cigar};
    while(ss >> a >>c){
        switch (c)
        {
        case 'S':
            asp.slen = a;
            break;
        case 'M':
        case 'D':
            asp.offset += a;
            break;
        case 'N':
            asp.altsp.push_back(make_pair(asp.offset,asp.offset+a));
            asp.offset += a;
            break;
        default:
            break;
        }
    }
    return asp;
}


//处理split-read, 对一断点的soft-clip seq 进行堆叠；
//split-read,最少12bp
//@chrom: read map 的染色体号；
//@pos: read map的位置；
//@seq: read 的序列；
//@cigar: read 的cigar值；
//@qname: read 的名字；
//@map_split_read: 记录每个断点堆叠信息的映射表；
//@return:  检查read 是否含soft-clip , 是返回 true ,若满足soft-clip 长度大于等于12 则进一步堆叠信息； 
bool parse_split_read(string& chrom,int pos,string& seq,string& cigar,string& qname,mps& map_split_read,unordered_map<string,vector<string>>& map_alt_split,unordered_map<string,string>& map_transcript ){
    smatch m;
    string split_read,match_seq;
    bool down = false;
    auto asp = getCigarInfo(cigar);

// 找可变剪切
    cout << "可变剪切： " << chrom << ":" << pos << "\t" << cigar<<"\t";
    for (auto as : asp.altsp){
        string p1 = chrom + ':' + to_string(as.first + pos);
        string p2 = chrom + ':' + to_string(as.second + pos);
        cout << p1 << "-" << p2 << "; ";
        if (map_transcript.find(p1) != map_transcript.end() && map_transcript.find(p2) != map_transcript.end()){
            map_alt_split[map_transcript[p1] + '-' + map_transcript[p2]].push_back(qname);
        }
    }
    cout << "\n";

    if(cigar[cigar.size()-1] == 'S'){
        // if(!isdigit(cigar[cigar.size()-3])) return false;
        down = true;
        // cout << "Treat Cigar: " << cigar << "\tsoft: " << rs.first << "\toffset: " <<  rs.second << "\n";
        pos += asp.offset ;
        int n = asp.slen;
        if(n<12) return false;
        split_read = seq.substr(seq.size()-n,n);
        if (seq.size()-10 < n) match_seq = seq.substr(0,seq.size()-n);
        else match_seq = seq.substr(seq.size()-n-10,10);
        
    } else {
        int n = asp.slen;
        if (n >= 12 ){
            split_read = seq.substr(0,n);
            reverse(split_read.begin(),split_read.end());
            pos++;
            match_seq = seq.substr(n,10);
        } else return false;
    }
    
    string p = chrom + ":" + to_string(pos); 
    for (char c:split_read) {
        if (c == 'N') return true;
    }

    pileup(map_split_read[p],split_read,match_seq,qname,cigar,down);
    return true;
}

//discordant read 记录到字典中（序列id -> 比对位置）
//方便后续使用该字典计算支持reads
//不记录 dup secondAln subAln 的 reads
void parse_discordant_read(string& qname,int flag,string& rname,int pos,string& nchrom,int npos,string& cigar,unordered_map<string,fusionPos>& map_dcp){
    if(flag & 3328 || map_c2i.find(nchrom) == map_c2i.end()) return;
    auto rs = getCigarInfo(cigar);
    int pos1 = pos - rs.offset;
    int pos2 = pos + rs.offset;

    if(map_c2i[rname] < map_c2i[nchrom] || (map_c2i[rname] == map_c2i[nchrom] && pos < npos)){
        map_dcp[qname].p1 = make_pair(rname,pos);
        map_dcp[qname].p2 = make_pair(nchrom,npos);
        map_dcp[qname].p1pos = make_pair(pos1,pos2);
    } else {
        map_dcp[qname].p2 = make_pair(rname,pos);
        map_dcp[qname].p1 = make_pair(nchrom,npos);
        map_dcp[qname].p2pos = make_pair(pos1,pos2);
    }
    
}

//将融合点位进行整合方便后续的检索比较。
//pos 除去 1000 记录这个区间的reads。
//减少后续的计算支持该断点融合的，discordant read 的数量。
//@map_dcp:  以discordant read 的 qname 为键，融合断点为值的字典；
//return map_fus_dcps: 汇总1000bp 范围内的discordant reads ,简化后续计算
unordered_map<string,vector<string>> combine_discordant_reads(unordered_map<string,fusionPos>& map_dcp){
    unordered_map<string,vector<string>> map_fus_dcps;
    for(auto x:map_dcp){
        int b1_pos = x.second.p1.second / 100000;
        int b2_pos = x.second.p2.second / 100000;
        string fus = x.second.p1.first + ":" + to_string(b1_pos) + "-" + x.second.p2.first + ":" + to_string(b2_pos);
        map_fus_dcps[fus].push_back(x.first);
    }
    return map_fus_dcps;
}

//整合split read 的break points
//使用 断点 上游-下游 作为Key 记录 fusion value 。
//对每一个split read，进行处理。融合结果更新到map_fusion 表中。
//对于soft-clip seq正向mapping 的 split-read, 设left_side(p1) 点位是，mapping 端在左侧的点位，反之为 right_side(p2);
//对于soft-clip seq反向互补mapping 的 split-read, 设染色体位点低的为 left-side, 染色体高的为 right-side；
//@break_p1 : 以 chrom:pos 记录的断点（检索bam 时得到的第一个断裂点）；
//@map_p2: 以<string(chrom),int(pos)> 记录的断点，是soft-clip reads 重新比对到的位置,以及比对的方式（0，1，2，3）；
// ==========xxxxxxxx   +(0);
// ==========xxxxxxxx   -(1);
// xxxxxxxxxx========   +(2);
// xxxxxxxxxx========   -(3);
//@num: 记录了这个断点-soft-clip reads 支持的数量；
//@map_fusion:  以break_point1-break_point2 为 key, fusion 为值的 散列表；
//@map_c2i: 染色体号和数字对应表，用于比较断裂点先后；
void combine_split_reads(string& break_p1,pair<pair<string,int>,int>& map_p2,Piled_reads& prd,map<string,fusion>& map_fusion, map<string,int>& map_c2i){
    pair<string,int> p1,p2,left_side,right_side;
    pair<char,char> mpdr;
    auto &break_p2 = map_p2.first;
    vector<string> v_b1 = split(break_p1,":");
    string chrom_b1 = v_b1[0];
    int p1n = 0;
    int p2n = 0;
    int pos_b1 = string2number<int>(v_b1[1]);
    pair<char,char> rpdr;
    if (map_c2i[chrom_b1] > map_c2i[break_p2.first] || (map_c2i[chrom_b1] == map_c2i[break_p2.first] && pos_b1 > break_p2.second)){
        p1 = break_p2;
        p2 = make_pair(chrom_b1,pos_b1);
        p2n = prd.count;
    } else {
        p1 = make_pair(chrom_b1,pos_b1);
        p2 = break_p2;
        p1n = prd.count;
    }

    mpdr = make_pair('+','+');
    rpdr = make_pair('u','d');
    if (map_p2.second == 1 || map_p2.second == 3){
        left_side = p1;
        right_side = p2;
        mpdr = make_pair('-','-');
        if(map_p2.second == 1) rpdr = make_pair('u','u');
        else rpdr = make_pair('d','d');
    } else if (map_p2.second == 0){
        left_side = make_pair(chrom_b1,pos_b1);
        right_side = break_p2;
        p1n = prd.count;
        p2n = 0;
    } else {
        left_side = break_p2;
        right_side = make_pair(chrom_b1,pos_b1);
        p2n = prd.count;
        p1n = 0;
    }

    // cout << break_p1 <<"\t" << break_p2.first <<":" << break_p2.second << "\n";
    string fusion_n = left_side.first + ":" + to_string(left_side.second) + "-" + right_side.first + ":" + to_string(right_side.second);
    fusion* fs =  &map_fusion[fusion_n];
    
    fs->p1 = left_side;
    fs->p2 = right_side;
    int trans_b1 = p1.second / 100000;
    int trans_b2 = p2.second / 100000;
    fs->check_point = p1.first + ":" + to_string(trans_b1) + "-" + p2.first + ":" + to_string(trans_b2);
    fs->p1n = fs->p1n + p1n;
    fs->p2n = fs->p2n + p2n;
    fs->mpdr = mpdr;
    fs->rpdr = rpdr;
    auto &vp = fs->vp ;
    for(auto s:prd.qnames) vp.push_back(s);
    // cigars 记录到对应的集合
    if(right_side == break_p2) {
        for(auto s:prd.cigars) fs->p1_cigars.push_back(s);
        fs->split_seqs.first = prd.seq;
    } else{
        for(auto s:prd.cigars) fs->p2_cigars.push_back(s);
        fs->split_seqs.second = prd.seq;
    }
}

//根据split read 断点信息，寻找支持的 discordant reads;
//split read 的 fusion 信息记录在 散列表 map_fusion 中；
//储备设定在断裂点附近。目前粗略计算在他们附近100bp,都作为支持项；
void combine_fusion(map<string,fusion>& map_fusion,unordered_map<string,vector<string>> map_fus_dcps,unordered_map<string,fusionPos>& map_dcp){
    int dis = 2000;
    for(auto& fs:map_fusion){
        string check_point = fs.second.check_point;
        // cout << "Check Point: " << check_point << "\t" << fs.second.p1.first << ":" << fs.second.p1.second << "-"
        //     << fs.second.p2.first << ":" << fs.second.p2.second << "\n";
        if(map_fus_dcps.find(check_point) != map_fus_dcps.end()){
            int n = 0;
            for(string& qname:map_fus_dcps[check_point]){
                fusionPos tmp_p = map_dcp[qname];
                // cout << qname << "\t" << tmp_p.p1.first << ":" << tmp_p.p1pos.first << "-" << tmp_p.p1pos.second << "\t"
                //     << tmp_p.p2.first << ":" << tmp_p.p2pos.first << "-" << tmp_p.p2pos.second ;
                auto &sp = fs.second;
                if(tmp_p.p1.first == sp.p1.first && tmp_p.p2.first == sp.p2.first){
                    if( (abs(tmp_p.p1pos.first - sp.p1.second) < dis || abs(tmp_p.p1pos.second - sp.p1.second) < dis ) && (abs(tmp_p.p2pos.first - sp.p2.second) < dis || abs(tmp_p.p2pos.second - sp.p2.second) < dis))
                    {
                        n++;
                        cout << "\tOK";
                        continue;
                    }
                }
                if (tmp_p.p2.first == sp.p1.first && tmp_p.p1.first == sp.p2.first){
                    if( (abs(tmp_p.p2pos.first - sp.p1.second) < dis || abs(tmp_p.p2pos.second - sp.p1.second) < dis ) && (abs(tmp_p.p1pos.first - sp.p2.second) < dis || abs(tmp_p.p1pos.second - sp.p2.second) < dis))
                    {
                        n++;
                        cout <<"\tOK";
                    }
                }
                cout << "\n";
            }
            fs.second.dcp = n;
        }
    }
}

int main(int argc,char* argv[]){
    cout << "OK\n";
    return 0;
}
