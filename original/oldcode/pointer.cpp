#include<iostream>

using namespace std;

int main() {

	double v[3] = {1,2,3};
	double* p = v;
	cout<<*(p++)<<endl;
	p=v;
	cout<<*p++<<endl;
	
	return 0;
	}
