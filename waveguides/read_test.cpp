#define real_t   double //float
#include <fstream>
#include <iostream>
using namespace std;
int main()
{
 real_t density;
 ifstream skvrna("/r/mhs_data.out", ios::in | ios::binary);
 skvrna.seekg (0, ios::beg);
 skvrna.read (reinterpret_cast<char *>(&density), sizeof (real_t));
 std::cout<<density;
 skvrna.read (reinterpret_cast<char *>(&density), sizeof (real_t));
 std::cout<<density;
 return 0;   
}