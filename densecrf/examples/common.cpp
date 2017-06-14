/*
    Copyright (c) 2013, Philipp Krähenbühl
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the Stanford University nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY Philipp Krähenbühl ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL Philipp Krähenbühl BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "common.h"
#include <fstream>
#include <iostream>
#include <io.h>
using namespace std;


// Store the colors we read, so that we can write them again.
int nColors = 0;
int colors[255];
int getColor( const unsigned char * c ){
	return c[0] + 256*c[1] + 256*256*c[2];
}
void putColor( unsigned char * c, int cc ){
	c[0] = cc&0xff; c[1] = (cc>>8)&0xff; c[2] = (cc>>16)&0xff;
}
// Produce a color image from a bunch of labels
unsigned char * colorize( const VectorXs & labeling, int W, int H ){
	unsigned char * r = new unsigned char[ W*H*3 ];
	for( int k=0; k<W*H; k++ ){
		int c = colors[ labeling[k] ];
		putColor( r+3*k, c );
	}
	//printf("%d %d %d \n",r[0],r[1],r[2]);
	return r;
}
// Read the labeling from a file
VectorXs getLabeling( const unsigned char * im, int N, int M ){
	VectorXs res(N);
	//printf("%d %d %d \n",im[0],im[1],im[2]);
	for( int k=0; k<N; k++ ){
		// Map the color to a label
		int c = getColor( im + 3*k );
		int i;
		for( i=0;i<nColors && c!=colors[i]; i++ );
		if (c && i==nColors){
			if (i<M)
				colors[nColors++] = c;
			else
				c=0;
		}
		res[k] = c?i:-1;
	}
	return res;
} 

void load_pointcloud(MatrixXf& pt, MatrixXf& normal,
	VectorXs& label, const string& filename)
{
	ifstream infile(filename, ios::binary);
	if (!infile) cout << "Open file error: " << filename << endl;

	int n;
	infile.read((char*)(&n), sizeof(int));

	pt.resize(3, n);
	infile.read((char*)pt.data(), sizeof(float) * 3 * n);

	normal.resize(3, n);
	infile.read((char*)normal.data(), sizeof(float) * 3 * n);

	label.resize(n);
	vector<int> buffer(n);
	infile.read((char*)buffer.data(), sizeof(int)*n);
	for (int i = 0; i < n; ++i)
	{
		label[i] = buffer[i];
	}

	infile.close();
}

void load_prob(MatrixXf& prob, const string& filename)
{
	ifstream infile(filename, ios::binary);
	if (!infile) cout << "Open file error: " << filename << endl;

	int c, n;
	infile.read((char*)(&c), sizeof(int));
	infile.read((char*)(&n), sizeof(int));

	prob.resize(c, n);
	infile.read((char*)prob.data(), sizeof(int)*prob.size());

	infile.close();
}

void save_label(const VectorXs& map, const string& filename)
{
	ofstream outfile(filename);
	for (int i = 0; i < map.size(); ++i)
	{
		outfile << map[i] << endl;
	}
	outfile.close();
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
		_all_filenames.push_back(file_path + std::string(c_file.name));
	} while (_findnext(hFile, &c_file) == 0);
	_findclose(hFile);
}
