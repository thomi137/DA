#include<limits>
#include<iostream>

using namespace std;

int main(){
cout << "Numeric epsilon on system" << endl;
cout << "double: "<< numeric_limits<double>::epsilon()<<endl;
cout << "int: "<< numeric_limits<int>::epsilon()<<endl;
cout << "float: "<< numeric_limits<float>::epsilon()<<endl;
}
