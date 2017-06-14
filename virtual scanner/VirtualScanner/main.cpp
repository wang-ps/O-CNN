#include <iostream>
#include <string>
#include <time.h>
#include <io.h>
#include "VirtualScanner.h"
using namespace std;

void get_all_filenames(vector<string>& _all_filenames, string _filename);
int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		cout << "Usage: VirtualScanner.exe <file name> "
			"[view_num] [flags]" << endl;
		return 0;
	}
	string filename(argv[1]);

	int view_num = 6; // scanning view number
	if (argc >= 3) view_num = atoi(argv[2]);
	
	bool flag = false; // output normal flipping flag
	if (argc >= 4) flag = atoi(argv[3]);

	vector<string> all_files;
	get_all_filenames(all_files, filename);

	//#pragma omp parallel for
	for (int i = 0; i < all_files.size(); i++)
	{
		clock_t t1 = clock();
		VirtualScanner scanner;
		scanner.scanning(all_files[i], view_num, flag);
		clock_t t2 = clock();

		string messg = all_files[i].substr(all_files[i].rfind('\\') + 1) +
			" done! Time: " + to_string(t2 - t1) + "\n";
		cout << messg;
	}

	return 0;
}

void get_all_filenames(vector<string>& _all_filenames, string _filename)
{
	// reset data
	_all_filenames.clear();

	// find
	size_t p0 = _filename.rfind('\\') + 1;
	size_t p1 = _filename.rfind('.');

	// file path
	string file_path(_filename, 0, p0);

	// get the regular expression
	_filename.replace(p0, p1 - p0, "*");

	// find all the file
	_finddata_t c_file;
	intptr_t hFile = _findfirst(_filename.c_str(), &c_file);
	do
	{
		if (hFile == -1) break;
		_all_filenames.push_back(file_path + string(c_file.name));
	} while (_findnext(hFile, &c_file) == 0);
	_findclose(hFile);
}