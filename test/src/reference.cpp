#include "shared/reference.h"

using namespace std;

map<string,vector<unsigned>> fastqReader::get_index(string& fasta_name){
        ifstream fai {fasta_name + ".fai",ifstream::binary};
        string chrom;
        unsigned seq_len,offset,line_bases,line_bytes;
        map<string,vector<unsigned>> index_dic;
        for(;fai >> chrom >> seq_len >> offset >> line_bases >> line_bytes;)
            index_dic[chrom] = vector<unsigned>{seq_len,offset,line_bases,line_bytes};
        fai.close();
        return index_dic;
    }

/*
   def fech(self,chrom,start,end):
        start -= 1
        end -= 1 # 索引是以0位开始计算碱基位，即 pos 是0 时对应第一个碱基，但是IGV 中，以及使用习惯上，人们以 1 为第一个碱基。
        assert chrom in self.index_dic, "chrom is not in fai index!\n"
        seq_len,offset,line_base,line_byte = self.index_dic[chrom]
        assert start >= 0 ,"start need bigger than 0 \n"
        assert start <= end , "start is bigger than end!\n"
        assert seq_len >= end, "end is bigger than seq length!\n"
        chunk_size = end - start + 1
        lines = int(start/line_base)
        pos = start % line_base
        shift = offset + lines*line_byte + pos  
        chunk_size  = chunk_size + int((chunk_size+pos)/line_byte) # 补充pos位，以计算跨越的行，加上每一行的LF
        self.file.seek(shift)
        chunk = self.file.read(chunk_size)
        return chunkParser(chunk)
    
*/
ChunkParse fastqReader::fetch(string chrom,int start,int end)
{
    unsigned seq_len,offset,line_bases,line_bytes;
    if(end < start){
        char* buffer = new char[0];
        return ChunkParse(buffer,0);        
    }
    vector<unsigned> v = index_dic[chrom];   
    seq_len = v[0];
    offset = v[1];
    line_bases = v[2];
    line_bytes = v[3];
    start--;
    end --;
    int chunk_size = end - start + 1;
    int lines = start/line_bases;
    int pos = start % line_bases;
    unsigned shift = offset + lines*line_bytes + pos;
//    cout << start << '/' << line_bases << '=' << lines<<'\n';
//    cout << start << '%' << line_bases << '=' << pos << '\n';
//    cout << offset << '+' << lines << '*' << line_bytes << '+' << pos << '=' << shift << '\n';
    chunk_size = chunk_size + (chunk_size+pos-1)/line_bases;
    fasta.seekg(shift);
    char* buffer = new char[chunk_size];
    fasta.read(buffer,chunk_size);
    return ChunkParse(buffer,chunk_size);
}

string ChunkParse::seq()
{
    string s;
    for(int i=0;i<length;i++)
    {
        if(chunk[i] != '\n') s += chunk[i];
    }
    return s;
}

// for string delimiter
vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}
