#include <fstream>
#include <iostream>
#include "VirtualScanner.h"

const Float EPS = 1.0e-20;

bool read_obj(std::string filename, MatrixXF& V, MatrixXi& F)
{
	std::ifstream infile(filename, std::ifstream::binary);
	if (!infile)
	{
		std::cout << "Open OBJ file error!" << std::endl;
		return false;
	}

	// get length of file
	infile.seekg(0, infile.end);
	int len = infile.tellg();
	infile.seekg(0, infile.beg);

	// load the file into memory
	char* buffer = new char[len + 1];
	infile.read(buffer, len);
	buffer[len] = 0;
	infile.close();

	// parse buffer data
	std::vector<char*> pVline, pFline;
	char* pch = strtok(buffer, "\n");
	while (pch != nullptr)
	{
		if (pch[0] == 'v' && pch[1] == ' ')
		{
			pVline.push_back(pch + 2);
		}
		else if (pch[0] == 'f' && pch[1] == ' ')
		{
			pFline.push_back(pch + 2);
		}

		pch = strtok(nullptr, "\n");
	}

	// load V
	V.resize(3, pVline.size());
	#pragma omp parallel for
	for (int i = 0; i < pVline.size(); i++)
	{
		char* p = strtok(pVline[i], " ");
		for (int j = 0; j < 3; j++)
		{
			V(j, i) = atof(p);
			p = strtok(nullptr, " ");
		}
	}

	// load F
	F.resize(3, pFline.size());
	#pragma omp parallel for
	for (int i = 0; i < pFline.size(); i++)
	{
		char* p = strtok(pFline[i], " ");
		for (int j = 0; j < 3; j++)
		{
			F(j, i) = atoi(p) - 1;
			p = strtok(nullptr, " ");
		}
	}

	// release
	delete[] buffer;
	return true;
}

bool read_off(std::string filename, MatrixXF& V, MatrixXi& F)
{
	std::ifstream infile(filename, std::ios::binary);
	if (!infile)
	{
		std::cout << "Open " + filename + " error!" << std::endl;
		return false;
	}

	// eat head
	std::string head;
	infile >> head;
	if (!(head == "OFF"))
	{
		std::cout << filename + " is not an OFF file!" << std::endl;
		return false;
	}

	// face/vertex number
	int nv, nf, ne;
	infile >> nv >> nf >> ne;

	// get length of file
	int p1 = infile.tellg();
	infile.seekg(0, infile.end);
	int p2 = infile.tellg();
	infile.seekg(p1, infile.beg);
	int len = p2 - p1;

	// load the file into memory
	char* buffer = new char[len + 1];
	infile.read(buffer, len);
	buffer[len] = 0;

	// close file
	infile.close();

	// parse buffer data
	std::vector<char*> pV;
	pV.reserve(3 * nv);
	char* pch = strtok(buffer, " \r\n");
	pV.push_back(pch);
	for (int i = 1; i < 3 * nv; i++)
	{
		pch = strtok(nullptr, " \r\n");
		pV.push_back(pch);
	}
	std::vector<char*> pF;
	pF.reserve(3 * nf);
	for (int i = 0; i < nf; i++)
	{
		// eat the first data
		pch = strtok(nullptr, " \r\n");
		for (int j = 0; j < 3; j++)
		{
			pch = strtok(nullptr, " \r\n");
			pF.push_back(pch);
		}
	}

	// load vertex
	V.resize(3, nv);
	Float* p = V.data();
	#pragma omp parallel for
	for (int i = 0; i < 3 * nv; i++)
	{
		*(p + i) = atof(pV[i]);
	}

	// load face
	F.resize(3, nf);
	int* q = F.data();
	#pragma omp parallel for
	for (int i = 0; i < 3 * nf; i++)
	{
		*(q + i) = atoi(pF[i]);
	}

	//release
	delete[] buffer;
	return true;
}

bool read_mesh(std::string filename, MatrixXF& V, MatrixXi& F)
{
	size_t found = filename.rfind('.');
	if (found != std::string::npos)
	{
		std::string suffix(filename, found + 1);
		std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);

		bool succ = false;
		if (suffix == "obj")
		{
			succ = read_obj(filename, V, F);
		}
		else if (suffix == "off")
		{
			succ = read_off(filename, V, F);
		}
		else
		{
			cout << "Error : Unsupported file formate!" << endl;
		}
		return succ;
	}
	return false;
}

void compute_face_normal(MatrixXF& Nf, VectorXF& f_areas,
	const MatrixXF& V, const MatrixXi& F)
{
	size_t nf = F.cols();
	Nf.resize(3, nf);
	f_areas.resize(nf);

	#pragma omp parallel for
	for (int i = 0; i < nf; i++)
	{
		Vector3F p01 = V.col(F(1, i)) - V.col(F(0, i));
		Vector3F p02 = V.col(F(2, i)) - V.col(F(0, i));
		Vector3F n = p01.cross(p02);

		Float d = n.norm();
		if (d < EPS) d = EPS;
		n /= d;

		f_areas[i] = d*0.5;
		Nf.col(i) = n;
	}
}

void delete_faces(MatrixXF& V, MatrixXi& F, const VectorXb& valid_f)
{
	int nf = F.cols();
	int nv = V.cols();

	// valid F idx
	int nvf = 0;
	VectorXi fi = VectorXi::Constant(F.cols(), -1);
	for (int i = 0; i < F.cols(); i++)
	{
		if (valid_f[i] == 1)
		{
			fi[i] = nvf++;
		}
	}

	// valid V idx
	VectorXb valid_v = VectorXb::Constant(nv, false);
	for (int i = 0; i < nf; ++i)
	{
		if (valid_f[i])
		{
			for (int j = 0; j < 3; ++j)
			{
				valid_v[F(j, i)] = true;
			}
		}
	}
	int nvv = 0;
	VectorXi vi = VectorXi::Constant(V.cols(), -1);
	for (int i = 0; i < V.cols(); i++)
	{
		if (valid_v[i])
		{
			vi[i] = nvv++;
		}
	}

	// valid F
	for (int i = 0; i < nf; i++)
	{
		if (valid_f[i])
		{
			F.col(fi[i]) << vi[F(0, i)], vi[F(1, i)], vi[F(2, i)];
		}
	}
	F.conservativeResize(3, nvf);

	// valid V
	for (int i = 0; i < nv; i++)
	{
		if (valid_v[i])
		{
			V.col(vi[i]) = V.col(i);
		}
	}
	V.conservativeResize(3, nvv);
}

void remove_zero_faces(MatrixXF& V, MatrixXi& F)
{
	MatrixXF N;
	VectorXF f_area;
	compute_face_normal(N, f_area, V, F);

	int k = 0;
	int nf = F.cols();
	VectorXb valid_f = VectorXb::Constant(nf, true);
	for (int j = 0; j < nf; ++j)
	{
		if (f_area[j] < 1.0e-15)
		{
			valid_f[j] = false;
			k++;
		}
	}
	if (k > 0)
	{
		delete_faces(V, F, valid_f);
	}
}

void normalize_mesh(MatrixXF& V)
{
	int size = V.size();
	auto buffer = V.data();
	float maxExtent = 0.0;
	for (int j = 0; j < size; ++j)
	{
		float vertexExtent = fabs(buffer[j]);
		if (vertexExtent > maxExtent)
		{
			maxExtent = vertexExtent;
		}
	}
	for (int j = 0; j < size; ++j)
	{
		buffer[j] /= maxExtent;
	}
}

bool VirtualScanner::scanning(string filename, int view_num, bool flags, bool normalize)
{
	if (view_num < 1) view_num = 1;
	if (view_num > total_view_num_) view_num = total_view_num_;

	// load mesh
	bool succ = read_mesh(filename, V_, F_);
	if (!succ) return false;
	
	if (normalize)
	{
	    normalize_mesh(V_);
	}

	// remove zero-area faces and unreferenced vertices,
	// otherwise, CGAL will collapse
	remove_zero_faces(V_, F_);

	// calc normal
	VectorXF f_area;
	compute_face_normal(N_, f_area, V_, F_);
	
	// calc views
	bbMin_ = V_.rowwise().minCoeff();
	bbMax_ = V_.rowwise().maxCoeff();
	calc_views();

	// data type conversion
	vector<Point> vecP; vecP.reserve(V_.cols());
	for (int i = 0; i < V_.cols(); i++)
	{
		vecP.push_back(Point(V_(0, i), V_(1, i), V_(2, i)));
	}
	vector<Triangle> vecF; vecF.reserve(F_.cols());
	for (int i = 0; i < F_.cols(); i++)
	{
		vecF.push_back(Triangle(vecP[F_(0, i)], vecP[F_(1, i)], vecP[F_(2, i)]));
	}

	// constructs AABB tree
	Tree tree(vecF.begin(), vecF.end());
	auto start_iteration = vecF.begin();

	// scanning
	int W = 2 * resolution_ + 1;
	int Num = view_num*W*W;
	MatrixXF buffer_pt(3, Num), buffer_n(3, Num);
	VectorXi buffer_flag = VectorXi::Constant(Num, 0);
	#pragma omp parallel for
	for (int i = 0; i < Num; ++i)
	{
		int y = i % W - resolution_;
		int j = i / W;
		int x = j % W - resolution_;
		int v = j / W;

		Vector3F pt = view_center_[v] + 
			Float(x)*dx_[v] + Float(y)*dy_[v];
		
		Point p0(pt[0], pt[1], pt[2]);
		Point p1(pt[0] - view_dir_[v][0], 
			pt[1] - view_dir_[v][1], pt[2] - view_dir_[v][2]);
		Ray ray_query(p0, p1); // ray with source p0 and passing through p1
		std::list<Ray_intersection> intersections;
		tree.all_intersections(ray_query, std::back_inserter(intersections));
		//Ray_intersection intersection = tree.first_intersection(ray_query);

		if (!intersections.empty())
		{
			Vector3F rst(-1.0e30, -1.0e30, -1.0e30);
			int id = -1;
			Float distance = 1.0e30;
			for (auto& intersection : intersections )
			{
				// gets intersection object
				Point* pPt = boost::get<Point>(&(intersection->first));
				if (!pPt) continue;

				Vector3F tmp(pPt->x(), pPt->y(), pPt->z());
				Float dis = (tmp - pt).squaredNorm();
				if (dis < distance)
				{
					distance = dis;
					rst = tmp;
					id = std::distance(start_iteration, intersection->second);
				}
			}
			if (id < 0 || rst[0]<bbMin_[0] - 1 || rst[0]>bbMax_[0] + 1 ||
				rst[1]<bbMin_[1] - 1 || rst[1]>bbMax_[1] + 1 ||
				rst[2]<bbMin_[2] - 1 || rst[2]>bbMax_[2] + 1) continue;

			// normal
			Vector3F normal = N_.col(id);
			int flag = 1;
			if (normal.dot(view_dir_[v]) < 0)
			{
				normal = -normal;
				flag = -1;
			}

			// save to buffer
			buffer_pt.col(i) = rst;
			buffer_n.col(i) = normal;
			buffer_flag(i) = flag;			
		}
	}

	// get points and normals
	pts_.clear();
	normals_.clear();
	for (int i = 0; i < Num; ++i)
	{
		if (buffer_flag(i) != 0)
		{
			for (int j = 0; j < 3; ++j)
			{
				pts_.push_back(buffer_pt(j, i));
				normals_.push_back(buffer_n(j, i));
			}
			flags_.push_back(buffer_flag(i));
		}
	}
	
	// save
	string filename_pc = filename.substr(0, filename.rfind('.'));
	//save_ply(filename_pc + "_pc.ply");
	save_binary(filename_pc + ".points", flags);
	return true;
}

void VirtualScanner::calc_views()
{
	view_center_.resize(total_view_num_);
	view_dir_.resize(total_view_num_);
	dx_.resize(total_view_num_);
	dy_.resize(total_view_num_);
	Vector3F center = (bbMin_ + bbMax_) * 0.5;
	Float len = (bbMax_ - bbMin_).norm() * 0.5;

	// view directions
	view_dir_[0] << 0.0, 0.0, 1.0;
	view_dir_[1] << 1.0, 0.0, 0.0;
	view_dir_[2] << -1.0, 0.0, 0.0;
	view_dir_[3] << 0.0, 1.0, 0.0;
	view_dir_[4] << 0.0, -1.0, 0.0;
	view_dir_[5] << 0.0, 0.0, -1.0;
	view_dir_[6] = (Vector3F(bbMax_[0], bbMax_[1], bbMax_[2]) - center).normalized();
	view_dir_[7] = (Vector3F(bbMax_[0], bbMin_[1], bbMax_[2]) - center).normalized();
	view_dir_[8] = (Vector3F(bbMin_[0], bbMax_[1], bbMax_[2]) - center).normalized();
	view_dir_[9] = (Vector3F(bbMin_[0], bbMin_[1], bbMax_[2]) - center).normalized();
	view_dir_[10] = (Vector3F(bbMax_[0], bbMax_[1], bbMin_[2]) - center).normalized();
	view_dir_[11] = (Vector3F(bbMax_[0], bbMin_[1], bbMin_[2]) - center).normalized();
	view_dir_[12] = (Vector3F(bbMin_[0], bbMax_[1], bbMin_[2]) - center).normalized();
	view_dir_[13] = (Vector3F(bbMin_[0], bbMin_[1], bbMin_[2]) - center).normalized();

	// view centers
	for (int i = 0; i < total_view_num_; i++)
	{
		view_center_[i] = view_dir_[i] * len + center;
	}

	// dx & dy
	Float d = (bbMax_ - bbMin_).maxCoeff() / (2 * resolution_ + 1);
	dx_[0] << d, 0, 0;
	dy_[0] << 0, d, 0;
	dx_[5] << d, 0, 0;
	dy_[5] << 0, d, 0;
	for (int i = 1; i < total_view_num_; ++i)
	{
		if (i == 5)continue;
		Vector3F Z(0.0, 0.0, 1.0);
		Float dot01 = Z.dot(view_dir_[i]);
		if (dot01 < -1)dot01 = -1;
		if (dot01 > 1)dot01 = 1;
		Float angle = acos(dot01);
		Vector3F axis = Z.cross(view_dir_[i]).normalized();

		Matrix3F Rot; Rot = Eigen::AngleAxis<Float>(angle, axis);

		dx_[i] = Rot.col(0)*d;
		dy_[i] = Rot.col(1)*d;
	}
}

bool VirtualScanner::save_ply(string filename)
{
	ofstream outfile(filename, ios::binary);
	if (!outfile)
	{
		cout << "Open " << filename << "error!" << endl;
		return false;
	}

	// write header
	int n = pts_.size() / 3;
	outfile << "ply" << endl
		<< "format ascii 1.0" << endl
		<< "element vertex " << n << endl
		<< "property float x" << endl
		<< "property float y" << endl
		<< "property float z" << endl
		<< "property float nx" << endl
		<< "property float ny" << endl
		<< "property float nz" << endl
		<< "element face 0" << endl
		<< "property list uchar int vertex_indices" << endl
		<< "end_header" << endl;

	// wirte contents
	const int len = 128;
	char* buffer = new char[n*len];
	char* pstr = buffer;
	#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		sprintf(pstr + i*len,
			"%.6f %.6f %.6f %.6f %.6f %.6f\n",
			pts_[3 * i], pts_[3 * i + 1], pts_[3 * i + 2],
			normals_[3 * i], normals_[3 * i + 1], normals_[3 * i + 2]);
	}
	
	int k = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = len*i; j<len*(i+1); ++j)
		{
			if (buffer[j] == 0) break;
			buffer[k++] = buffer[j];
		}
	}
	outfile.write(buffer, k);
	
	outfile.close();
	delete[] buffer;
	return true;
}

bool VirtualScanner::save_binary(string filename, bool flags)
{
	ofstream outfile(filename, ios::binary);
	if (!outfile)
	{
		cout << "Opening file error!" << endl;
		return false;
	}
	
	int n = pts_.size();

	int num = n / 3;  // point number
	outfile.write((char*)(&num), sizeof(int));
	
 	//int pt_n = 1 | 2; // has_point | has_normal
 	//outfile.write((char*)(&pt_n), sizeof(int));
	
	outfile.write((char*)pts_.data(), sizeof(float)*n);
	outfile.write((char*)normals_.data(), sizeof(float)*n);	
	
	outfile.close();

	if (flags)
	{
		string filename_flag = filename;
		filename_flag.replace(filename.rfind('.') + 1, string::npos, "flags");
		outfile.open(filename_flag, ios::binary);
		outfile.write((char*)flags_.data(), sizeof(int)*flags_.size());
		outfile.close();
	}

	return true;
}
