/*
	compare field 1 and field 3 char by char
	Bo W Han 07/06/2014
*/
#include <cstdlib> // for exit
#include <iostream>
#include <fstream>
using namespace std;
const int str_len = 13;
int main(int argc, char** argv)
{
	if(argc==1) {
		cerr << "this program compare field 1 and 3 char by char and count"
			<< " write positions with identify char to an integer and record"
			<< " it in a new field\n";
		cerr << "usage: " << argv[0] << " input.txt" << endl;
		exit(0);
	}
	ifstream in(argv[1]);
	if(!in) {
		cerr << "cannot open file " << argv[1] << endl; exit(1);
	}
	unsigned int bits[str_len];
	for(int i=0;i<str_len;++i) {
		bits[i] = 1<< (str_len-i-1);
	}
	char buffer[1024];
	char* ptr1;
	char* ptr2;
	unsigned int bit;
	while(in.getline(buffer, 1024, '\n')) {
		ptr1 = ptr2 = buffer;
		bit = 0;
		while(*ptr2++!='\t');
		while(*ptr2++!='\t');
		for(int i=0;i<str_len;++i){
			if(*ptr1++==*ptr2++) {
				bit |= bits[i];
			}
		}
		cout << buffer << '\t' << bit << '\n';
	}
}
