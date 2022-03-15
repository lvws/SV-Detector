#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include<vector>
#include<sstream>
#include<typeinfo>

using namespace std;

class ChunkParse
{
public:
    ChunkParse(char* chunk_in,int len):
        chunk{chunk_in},length{len}{}
    ~ChunkParse(){ delete[] chunk;}
    string seq();

private:
    char* chunk;
    int length;
};

class fastqReader
{
public:
    fastqReader(string fasta_name)
        :fasta{fasta_name}{
            index_dic = get_index(fasta_name);
        }
    ~fastqReader(){fasta.close();}
    ChunkParse fetch(string chrom,int start,int end);
private:
    ifstream fasta;
    map<string,vector<unsigned>> index_dic;
    map<string,vector<unsigned>> static get_index(string& fasta_name);
};

vector<string> split (string s, string delimiter);
