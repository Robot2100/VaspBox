#pragma once
#include "stdafx.h"
using namespace std;
flo _quad(flo);
int _sign(flo);

#ifndef POINT_H
struct Point
{
	flo x;
	flo y;
	flo z;
	flo r() const
	{
		return sqrt(_quad(x) + _quad(y) + _quad(z));
	}
	Point() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
	Point(flo _x, flo _y, flo _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}
	Point operator-(Point & a)
	{
		Point r;
		r.x = -a.x + x;
		r.y = -a.y + y;
		r.z = -a.z + z;
		return r;
	}
	Point operator+(Point & a)
	{
		Point r(*this);
		r.x = a.x + x;
		r.y = a.y + y;
		r.z = a.z + z;
		return r;
	}
	Point operator*(const Point & a)
	{
		Point r(*this);
		r.x *= a.x;
		r.y *= a.y;
		r.z *= a.z;
		return r;
	}
	Point operator*(const flo & a)
	{
		Point r(*this);
		r.x *= a;
		r.y *= a;
		r.z *= a;
		return r;
	}
	Point operator*(flo && a)
	{
		Point r(*this);
		r.x *= a;
		r.y *= a;
		r.z *= a;
		return r;
	}
	Point & operator+=(const Point & a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	}
	Point & operator-=(const Point & a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;
		return *this;
	}
	bool operator!=(const Point & a)
	{
		return (a.x != x) || (a.y != y) || (a.z != z);
	}
	void normalize()
	{
		flo R = r();
		x /= R;
		y /= R;
		z /= R;
	}
	static Point corrf(Point & p1, Point & p2) {
		return Point(((p1.x)*(p2.x)), ((p1.y)*(p2.y)), ((p1.z)*(p2.z)));
	}
};
#endif
#ifndef CELL_MATRIX_H
const flo RadtoGrad = flo(90 / 1.57079632679489661923);
const flo GradtoRad = flo(1.57079632679489661923 / 90);
struct Matrix
{             // 0  1  2  3  4  5
	flo U[6]; // x2 y2 z2 xy xz yz
	Matrix()
	{
		for (int i = 0; i<3; i++)
			U[i] = 1;
		for (int i = 3; i<6; i++)
			U[i] = 0;

	}
	void RotateZ(flo angle)
	{
		flo Cos = cos(angle);
		flo Sin = sin(angle);

		flo temp[6];

		temp[0] = U[0] * _quad(Cos) + U[1] * _quad(Sin) - U[3] * Cos*Sin;
		temp[1] = U[0] * _quad(Sin) + U[1] * _quad(Cos) + U[3] * Cos*Sin;
		temp[2] = U[2];
		temp[3] = U[3] * (_quad(Cos) - _quad(Sin)) + 2 * Cos*Sin*(U[0] - U[1]);
		temp[4] = U[4] * Cos - U[5] * Sin;
		temp[5] = U[5] * Cos + U[4] * Sin;

		for (int i = 0; i<6; i++)
			U[i] = temp[i];
	}
	void RotateX(flo angle)
	{
		flo Cos = cos(angle);
		flo Sin = sin(angle);

		flo temp[6];

		temp[0] = U[0];
		temp[1] = U[1] * _quad(Cos) + U[2] * _quad(Sin) - U[0] * Cos*Sin;
		temp[2] = U[1] * _quad(Sin) + U[2] * _quad(Cos) + U[0] * Cos*Sin;
		temp[3] = U[3] * Cos - U[4] * Sin;
		temp[4] = U[3] * Cos + U[4] * Sin;
		temp[5] = U[5] * (_quad(Cos) - _quad(Sin)) + 2 * Cos*Sin*(U[1] - U[2]);

		for (int i = 0; i<6; i++)
			U[i] = temp[i];
	}
	Matrix Mirror() const
	{
		Matrix mat;
		mat.U[0] = 1 / U[0];
		mat.U[1] = 1 / U[1];
		mat.U[2] = 1 / U[2];
		mat.U[3] = -U[3] / (U[0] * U[1]);
		mat.U[4] = (U[3] * U[5] - U[4] * U[1]) / (U[0] * U[1] * U[2]);
		mat.U[5] = -U[5] / (U[2] * U[1]);
		return mat;
	}
	Point Transform(Point a) const
	{
		Point ret;
		ret.x = U[0] * a.x + U[3] * a.y + U[4] * a.z;
		ret.y = U[1] * a.y + U[5] * a.z;
		ret.z = U[2] * a.z;
		return ret;
	}
};
struct Cell
{
	flo a, b, c;
	flo alpha, beta, gamma;
	Cell(Matrix Mat) {
		for (int i = 0; i< 6; i++)
			if (Mat.U[i] == 0) Mat.U[i] = flo(0.0000000000000001);
		a = Mat.U[0];
		b = sqrt(_quad(Mat.U[1]) + _quad(Mat.U[3]));
		flo CosG = Mat.U[3] / b; //CosG
		flo SinG = sqrt(1 - _quad(CosG)); //SinG
		gamma = acos(CosG);

		flo temp = Mat.U[4] * CosG + Mat.U[5] * SinG;
		c = sqrt(_quad(Mat.U[2] * SinG + _quad(temp) - 2 * CosG*temp*Mat.U[4] + _quad(Mat.U[4]))) / SinG;
		beta = acos(Mat.U[4] / c);
		alpha = acos(cos(beta)*temp);
		toGrad();
	}
	Cell() {

	}
	void toGrad()
	{
		alpha *= RadtoGrad;
		beta *= RadtoGrad;
		gamma *= RadtoGrad;
	}
	void toRad()
	{
		alpha *= GradtoRad;
		beta *= GradtoRad;
		gamma *= GradtoRad;
	}
	Matrix CreateMatrix()
	{
		toRad();
		Matrix res;
		res.U[0] = a;
		res.U[1] = b*sin(gamma);
		res.U[2] = c*sqrt(1 - _quad(cos(alpha)) - _quad(cos(beta)) - _quad(cos(gamma)) + 2 * cos(alpha)*cos(beta)*cos(gamma)) / sin(gamma);
		res.U[3] = b*cos(gamma);
		res.U[4] = c*cos(beta);
		res.U[5] = c / sin(gamma)*(cos(alpha) - cos(beta)*cos(gamma));
		for (int i = 3; i < 6; i++)
			if (abs(res.U[i]) < (1E-10)) res.U[i] = 0;
		toGrad();
		return res;
	}
};
#endif
typedef Point Line;

const Point OX(1, 0, 0), OY(0, 1, 0), OZ(0, 0, 1);

const int MAX_LINE = 255;
#ifndef PARAM_H
class Param
{
public:
	string filename;
	Param() {
		filename = "a.txt";
	}
	void Take(int argn, char * argv[]) {
		if (argn == 1) return;
		filename = string(argv[1]);
	}
};
#endif

struct VecPoints
{
	VecPoints() {}
	VecPoints(char * a) {
		strcpy_s(name, a);
	}
	char name[3];
	vector<Point> points;

};
struct Supercell
{
	vector<VecPoints> vp;
	Supercell(Matrix & mat, vector<VecPoints> & vpIN)
	{
		vp = vpIN;
		for (int i = 0; i < (int)vpIN.size(); i++) {
			vp[i].points.clear();
			if (vpIN[i].name[0] == '\0') break;
			for (int j = 0; j < (int)vpIN[i].points.size(); j++) {
				for (int q = -1; q <= 1; q++) {
					for (int w = -1; w <= 1; w++) {
						for (int e = -1; e <= 1; e++) {
							//if (q == 0 && w == 0 && e == 0) continue;
							Point temp = vpIN[i].points[j] + Point((float)q, (float)w, (float)e);
							vp[i].points.push_back(mat.Transform(temp));
						}
					}
				}




			}
		}
		NCalls = 0;
	}
	int NCalls;
	void Uniq(int a, int b, flo bond)
	{
		b = b*27 + 13;
		vector<VecPoints> tempVP, outVP;
		tempVP = vp;
		for (int i = 0; i < (int)tempVP.size(); i++) {
			tempVP[i].points.clear();
		}
		outVP = tempVP;
		tempVP[a].points.push_back(vp[a].points[b]);
		ClearPoint(vp[a].points, b);
		int Count = 1;
		
		while (Count != 0) {
			for (int i = 0; i < (int)tempVP.size(); i++) {
				if (tempVP[i].name[0] == '\0') break;
				for (int j = 0; j < (int)tempVP[i].points.size(); j++) {

					for (int i1 = 0; i1 < (int)vp.size(); i1++) {
						if (vp[i1].name[0] == '\0') break;
						for (int j1 = 0; j1 < (int)vp[i1].points.size(); j1++) {
							if ((vp[i1].points[j1] - tempVP[i].points[j]).r() > (float)bond) continue;
							tempVP[i1].points.push_back(vp[i1].points[j1]);
							ClearPoint(vp[i1].points, j1);
							j1--;
							Count++;
						}
					}
					outVP[i].points.push_back(tempVP[i].points[j]);
					ClearPoint(tempVP[i].points, j);
					j--;
					Count--;
				}
			}
		}
		outVP.swap(vp);
	}
private:
	void ClearPoint(vector<Point> & a, int b)
	{
		NCalls++;
		vector<Point> temp;
		a.swap(temp);
		for (int i = 0; i < (int)temp.size(); i++) {
			if (i == b) continue;
			a.push_back(temp[i]);
		}
	}
};