#include <map>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;

int SW_map(string& seq1,string& seq2, int match = 1, int mismatch = -1, int delins = -2,int dis=5){
    map<pair<int,int>,pair<int,int>> map_path;
    map<pair<int,int>,int> map_score;
    int up = 0,left = 0 ,dd = 0,now = 0;
    dis = 7;
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
    auto chain = v_score.begin()->first;
    int s = v_score.begin()->second;
    int pre_i = chain.first ;
    int pre_j = chain.second ;
    for(;s > 0;){
        auto next_chain = map_path[chain];
        if (pre_i == next_chain.first ){
            cout << '*' << '-' << seq2[pre_j--] << '\t' << s  << '\n';
        }
        else if (pre_j == next_chain.second){
                cout << seq1[pre_i--] << '-' << "*\t" << s  << '\n';
        }
        else cout << seq1[chain.first] << '-' << seq2[chain.second] <<  "\t" <<  s <<  '\n';
        pre_i = next_chain.first;
        pre_j = next_chain.second;
        chain = next_chain;
        if (chain.first == -1) break;
        s = map_score[chain];
    }
    return v_score.begin()->second;

}

int main(int argc, char* argv[]){
    if(argc < 3){
        cout << argv[0] << "\tseq1\tseq2" << "\n";
        return -1;
    }
    string seq1 = argv[1];
    string seq2 = argv[2];
    int s = SW_map(seq1,seq2);
    cout << "Match Socre: " << s << '\n';

}