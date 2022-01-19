#include <iostream>
#include <string>
#include <stdint.h>
#include<vector>
#include <utility>


using namespace std;

struct postion
{
    uint8_t chrom;
    unsigned int pos;
};


int main(int argc,char *argv[])
{
    postion pos1{1,111111};
    int dint = 12;
    int8_t dint_8 = 12;
    char c = 'a';
    string s = "1:111111";
    int16_t dint16 = 12;
    pair<int8_t,int> p = make_pair(1,11111111);
    unsigned int ui = 1 ;
    ui = (ui & 0xFFFFFFFF) << 28;
    cout <<hex << ui << "\n";
    ui = ui + 123456789;
    unsigned int ch = (ui & 0xFF000000) >> 28;
    unsigned int pos = ui & 0x0FFFFFFF;
    cout << "Size of:\n"
    << "default int: " << sizeof(dint) << "\n"
    << "int 8: " << sizeof(dint_8) << "\n"
    <<"int 16: " << sizeof(dint16) << "\n"
    << "char: " << sizeof(c) << "\n"
    << "string: " << sizeof(s) << "\n"
    <<"pair: " << sizeof(p) << "\n"
    << "UI: " << ui << " Size:" << sizeof(ui) << "\n"
    << "Chrom: " << ch << " Pos: " << pos << "\n"
    <<"Struct pos: " << pos1.chrom << ":" << pos1.pos << " Size: " << sizeof(pos1)<<"\n";


}