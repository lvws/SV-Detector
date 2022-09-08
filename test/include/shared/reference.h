#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include<vector>
#include<sstream>
#include<typeinfo>

using namespace std;

class fastqReader
{
public:
    fastqReader(string fasta_name)
        :fasta{fasta_name}{
            index_dic = get_index(fasta_name);
            fasta.seekg(0,fasta.end);
            length = fasta.tellg();
            fasta.seekg(0,fasta.beg);
        }
    ~fastqReader(){fasta.close();}
    string fetch(string chrom,int start,int end);
private:
    unsigned length;
    ifstream fasta;
    map<string,vector<unsigned>> index_dic;
    map<string,vector<unsigned>> static get_index(string& fasta_name);
};

vector<string> split (string s, string delimiter);
