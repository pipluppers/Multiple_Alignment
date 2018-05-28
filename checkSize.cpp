#include <iostream>
#include <fstream>

using namespace std;

int main() {

	string input, str = "";
	cout << "Enter name of file: ";
	cin >> input;
	cout << endl;

	ifstream fin(input.c_str());
	while (fin >> input) {
		str += input;
	}
	cout << "Size of string: " << str.size() << endl;
	return 0;
}
