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
//@ptag: 是否为primary alignment
//return 更新堆叠序列集合。
void pileup(vector<Piled_reads>& v_ss,string& s,string& m,string& qname,string& cigar,bool down,bool ptag){
    bool flag = false;
    for(auto& p :v_ss){
        if (p.down != down) continue;        
        string hs = hanming(p.seq,s);  
        if(hs != "") {
            p.seq = hs;
            p.count ++;
            if (ptag) p.ucount ++; 
            p.qnames.push_back(qname);
            p.cigars.push_back(cigar);
            flag = true ;
        }
    }
    if (flag) return;
    int uc = 0;
    if (ptag) uc = 1;
    v_ss.push_back({s,m,1,uc,{qname},{cigar},down});
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
            asp.slen.push_back(a);
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

 
bool parse_split_read(string& chrom,int pos,string& seq,string& cigar,string& qname,mps& map_split_read,unordered_map<string,vector<string>>& map_alt_split,unordered_map<string,string>& map_transcript,
pair<unsigned,string>& p_sa,map<string,vector<pair<unsigned,bool>>>& map_sa,bool ptag ){
    smatch m;
    string split_read,match_seq;
    bool down = false;
    auto asp = getCigarInfo(cigar);
    if (asp.slen.size() == 0) return false; // 不存在soft-clip reads
    // 找可变剪切
    //cout << "可变剪切： " << chrom << ":" << pos << "\t" << cigar<<"\t";
    for (auto as : asp.altsp){
        string p1 = chrom + ':' + to_string(as.first + pos);
        string p2 = chrom + ':' + to_string(as.second + pos);
        //cout << p1 << "-" << p2 << "; ";
        if (map_transcript.find(p1) != map_transcript.end() && map_transcript.find(p2) != map_transcript.end()){
            map_alt_split[map_transcript[p1] + '-' + map_transcript[p2]].push_back(qname);
        }
    }
    //cout << "\n";
   
    // 检查soft 若有两个选择最长的, 若只有一个，则检查时在开端还是结尾
    bool stag = false; // soft 位于开头的标记
    int n = asp.slen[0];
    if ( asp.slen.size() > 1){
        if (asp.slen[0] >= asp.slen[1]) {
            stag = true;
        } else n = asp.slen[1];
    }  else if(cigar[cigar.size()-1] != 'S') stag = true;

    //处理reads信息
    if(!stag){
        // if(!isdigit(cigar[cigar.size()-3])) return false;
        down = true;
        // cout << "Treat Cigar: " << cigar << "\tsoft: " << rs.first << "\toffset: " <<  rs.second << "\n";
        pos += asp.offset ;
        if(n<12) return false;
        split_read = seq.substr(seq.size()-n,n);
        if (seq.size()-10 < n) match_seq = seq.substr(0,seq.size()-n);
        else match_seq = seq.substr(seq.size()-n-10,10);
    }                
    else {
        if (n >= 12 ){
            split_read = seq.substr(0,n);
            reverse(split_read.begin(),split_read.end());
            pos++;
            match_seq = seq.substr(n,10);
        } else return false;
    }

     //处理SA 信息
    if(p_sa.first !=0){
        string f_sa = chrom + ":" + to_string(pos/10);
        if(p_sa.second.back() == 'S'){
            auto asp_sa = getCigarInfo(p_sa.second);
            map_sa[f_sa].push_back(make_pair(p_sa.first+asp_sa.offset-1,false));
        }
        else map_sa[f_sa].push_back(make_pair(p_sa.first,true));
        //cout << "SA: " << f_sa << ": " << p_sa.first << "\t" << p_sa.second << "\n"; 
    }
    
    string p = chrom + ":" + to_string(pos); 
    // cout << "Count Split Reads: " << p << '\t' << qname << "\n";
    for (char c:split_read) {
        if (c == 'N') return true;
    }

    pileup(map_split_read[p],split_read,match_seq,qname,cigar,down,ptag);
    return true;
}

//discordant read 记录到字典中（序列id -> 比对位置）
//方便后续使用该字典计算支持reads
//不记录 dup secondAln subAln 的 reads
void parse_discordant_read(string& qname,int flag,string& rname,int pos,string& nchrom,int npos,string& cigar,unordered_map<string,fusionPos>& map_dcp){
    if(flag & 3328 || map_c2i.find(nchrom) == map_c2i.end()) return;
    auto rs = getCigarInfo(cigar);
    int pos1 = pos ;
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
        p2n = prd.ucount;
    } else {
        p1 = make_pair(chrom_b1,pos_b1);
        p2 = break_p2;
        p1n = prd.ucount;
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
        p1n = prd.ucount;
        p2n = 0;
    } else {
        left_side = break_p2;
        right_side = make_pair(chrom_b1,pos_b1);
        p2n = prd.ucount;
        p1n = 0;
    }

    //cout << break_p1 <<"\t" << break_p2.first <<":" << break_p2.second << "\n";
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
    
    if (left_side == break_p2){
        for(auto s:prd.qnames) fs->vp2.push_back(s);
    } else for(auto s:prd.qnames) fs->vp1.push_back(s);

    // cigars 记录到对应的集合
    if(right_side == break_p2) {
        for(auto s:prd.cigars) fs->p1_cigars.push_back(s);
        fs->split_seqs.first = prd.seq;
    } else{
        for(auto s:prd.cigars) fs->p2_cigars.push_back(s);
        fs->split_seqs.second = prd.seq;
    }
    //cout << "Combine split reads: " << fusion_n << "\t" << fs->p1n << ':' << fs->p2n <<'\t'<< fs->vp1.size() << ':' << fs->vp2.size() << '\n';
}

//根据split read 断点信息，寻找支持的 discordant reads;g
//储备设定在断裂点附近。目前粗略计算在他们附近100bp,都作为支持项；
void combine_fusion(map<string,fusion>& map_fusion,unordered_map<string,vector<string>> map_fus_dcps,unordered_map<string,fusionPos>& map_dcp){
    int dis1 = 3000; //远距离
    int dis2 = 20; // 近距离，防止融合位点偏移 
    for(auto& fs:map_fusion){
        string check_point = fs.second.check_point;
        // cout << "Check Point: " << check_point << "\t" << fs.second.p1.first << ":" << fs.second.p1.second << "-"
        //     << fs.second.p2.first << ":" << fs.second.p2.second << "\n";
        if(map_fus_dcps.find(check_point) != map_fus_dcps.end()){
            int n = 0;
            for(string& qname:map_fus_dcps[check_point]){
                fusionPos tmp_p = map_dcp[qname];
                // cout << qname << "\t" << tmp_p.p1.first << ":" << tmp_p.p1pos.first << "-" << tmp_p.p1pos.second << "\t"
                //      << tmp_p.p2.first << ":" << tmp_p.p2pos.first << "-" << tmp_p.p2pos.second ;
                auto &sp = fs.second;
                int check_num = 0;
                if(tmp_p.p1.first == sp.p1.first && tmp_p.p2.first == sp.p2.first){
                    // mapping seq 在下游，soft seq 在上游
                    if(sp.rpdr.first == 'd' && abs(tmp_p.p1pos.first - sp.p1.second) < dis1 && (tmp_p.p1pos.first + dis2) > sp.p1.second) check_num++;
                    // mapping seq 在上游，soft seq 在下游
                    else if (sp.rpdr.first == 'u' && abs(tmp_p.p1pos.second - sp.p1.second) < dis1 && (tmp_p.p1pos.second - dis2) < sp.p1.second ) check_num ++ ;

                    if(sp.rpdr.second == 'd' && abs(tmp_p.p2pos.first - sp.p2.second) < dis1 && (tmp_p.p2pos.first + dis2) > sp.p2.second ) check_num ++ ;
                    else if(sp.rpdr.second == 'u' && abs(tmp_p.p2pos.second - sp.p2.second) < dis1 && (tmp_p.p2pos.second - dis2) < sp.p2.second) check_num ++;

                    if(check_num == 2){
                        n++;
                        fs.second.vdcp.push_back(qname);
                        //cout << "\tOK!\n";
                        continue;
                    }
                    check_num = 0; // 防止第二次计算不符合，但积累了数值；
                }
                if (tmp_p.p2.first == sp.p1.first && tmp_p.p1.first == sp.p2.first){
                    // mapping seq 在下游，soft seq 在上游
                    if(sp.rpdr.first == 'd' && abs(tmp_p.p2pos.first - sp.p1.second) < dis1 && (tmp_p.p2pos.first + dis2) > sp.p1.second) check_num++;
                    // mapping seq 在上游，soft seq 在下游
                    else if (sp.rpdr.first == 'u' && abs(tmp_p.p2pos.second - sp.p1.second) < dis1 && (tmp_p.p2pos.second - dis2) < sp.p1.second ) check_num ++ ;

                    if(sp.rpdr.second == 'd' && abs(tmp_p.p1pos.first - sp.p2.second) < dis1 && (tmp_p.p1pos.first + dis2) > sp.p2.second ) check_num ++ ;
                    else if(sp.rpdr.second == 'u' && abs(tmp_p.p1pos.second - sp.p2.second) < dis1 && (tmp_p.p1pos.second - dis2) < sp.p2.second) check_num ++;

                    if(check_num == 2){
                        n++;
                        fs.second.vdcp.push_back(qname);
                        //cout << "\tOK!\n";
                        continue;
                    }
                }
                //cout << "\tfail!\n";
            }
            fs.second.dcp = n;
        }
    }
}

