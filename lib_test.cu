#include <iostream>

using namespace std;

void printTest() {
	#ifdef ONLYCPU
		cout << "Test 1" << endl;
	#else
		cout << "Test 2" << endl;
	#endif
}