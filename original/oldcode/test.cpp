#include<vector>
#include<iostream>
#include<algorithm>

using namespace std;

int main(){
	vector<double> one(4,4);
	vector<double> two(4,5);
	vector<double> three = one*two;
	
	copy(three.begin(), three.end(), ostream_iterator<double>(cout, "\n") )l
	
	return 0;
	
	}
	
