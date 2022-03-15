#include <sstream>
#include <iostream>
#include <string>
#include <map>
#include <regex>
#include <utility>
#include <vector>
#include <algorithm>

using namespace std;

int main(){
	vector<pair<int,int>> vps ;
	vps.push_back(make_pair(1,2));
	vps.push_back(make_pair(1,1));
	vps.push_back(make_pair(2,1));
	vps.push_back(make_pair(2,1));
	vps.push_back(make_pair(3,1));
	vps.push_back(make_pair(3,2));
	vps.push_back(make_pair(2,2));
	sort(vps.begin(),vps.end());
	for(auto p:vps) cout << p.first <<':' <<  p.second << "\t";
	auto it = unique(vps.begin(),vps.end());
	vps.resize(distance(vps.begin(),it));
	cout << "\n";
	for(auto p:vps) cout << p.first <<':' <<  p.second << "\t";
}
