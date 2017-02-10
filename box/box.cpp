// box.cpp: определ€ет точку входа дл€ консольного приложени€.
//

#include "stdafx.h"
#include "Classes.h"
using namespace std;

int main(int argn, char * argv[])
{
	flo _ce = 10.0;
	flo bond = 2.2;

	int _type=1, _number=1;

	char buf[MAX_LINE];
	vector<VecPoints> vp;
	Param param;
	param.Take(argn, argv);
	ifstream file(param.filename);
	Matrix mat;
	bool defenedCell = false;
	file >> buf; // —читывание €чейки
	if (strcmp(buf, "CELL") == 0) {
		Cell cell;
		file >> buf;

		file >> buf;
		cell.a = flo(atof(buf));

		file >> buf;
		cell.b = flo(atof(buf));

		file >> buf;
		cell.c = flo(atof(buf));


		file >> buf;
		cell.alpha = flo(atof(buf));

		file >> buf;
		cell.beta = flo(atof(buf));

		file >> buf;
		cell.gamma = flo(atof(buf));
		mat = cell.CreateMatrix();
		defenedCell = true;
	}
	else {
		file.close();
		cout << "Error! Can't read CELL line." << endl;
		exit(-1);
	}

	// ÷икл считывани€ атомов
	while (!file.eof())
	{
		file >> buf;
		if (buf[0] == '=')
			continue;
		int i = 0;
		for (; i < MAX_LINE; i++) {
			if ((buf[i]<'A' || buf[i]>'Z') && (buf[i]<'a' || buf[i]>'z')) {
				buf[i] = '\0';
				break;
			}
		}
		if (i > 3) exit(1);
		char atom[3];
		strcpy_s(atom, buf);

		file >> buf;
		if (buf[0] == ',')
			file >> buf;

		Point k;
		// чтение координат атома
		file >> buf;
		if (buf[1] == ',') buf[1] = '.';
		if (buf[2] == ',') buf[2] = '.';
		k.x = flo(atof(buf));

		file >> buf;
		if (buf[1] == ',') buf[1] = '.';
		if (buf[2] == ',') buf[2] = '.';
		k.y = flo(atof(buf));

		file >> buf;
		if (buf[1] == ',') buf[1] = '.';
		if (buf[2] == ',') buf[2] = '.';
		k.z = flo(atof(buf));
		for (i = 0; i < vp.size(); i++) {
			if (strcmp(vp[i].name, atom) == 0) {
				break;
			}
		}
		if (vp.size() == i) vp.push_back(VecPoints(atom));
		vp[i].points.push_back(k);
	}

	file.close();
	ofstream file2("POSCAR");
	file2 << "TITLE\n";
	file2 << "1.000\n";
	file2 << setprecision(6);
	file2.setf(ios::fixed, ios::floatfield);

	Point shift;
	int check = 0;
	bool _ip = false;
	Supercell supCell(mat, vp);
	supCell.Uniq(_type, _number, bond);
	vp = supCell.vp;
	//{
	//	_ip = false;
	//	shift = vp[_type].points[_number]-Point(0.5,0.5,0.5);
	//	for (int i = 0; i < vp.size(); i++) {
	//		if (vp[i].name[0] == '\0') break;
	//		for (int j = 0; j < vp[i].points.size(); j++) {
	//			vp[i].points[j] -= shift;
	//			if (vp[i].points[j].x > 1.0) { vp[i].points[j].x -= 1; _ip = true; }
	//			if (vp[i].points[j].y > 1.0) { vp[i].points[j].y -= 1; _ip = true; }
	//			if (vp[i].points[j].z > 1.0) { vp[i].points[j].z -= 1; _ip = true; }
	//			if (vp[i].points[j].x < 0) { vp[i].points[j].x += 1; _ip = true; }
	//			if (vp[i].points[j].y < 0) { vp[i].points[j].y += 1; _ip = true; }
	//			if (vp[i].points[j].z < 0) { vp[i].points[j].z += 1; _ip = true; }
	//		}
	//	}
	//}

	
	if (defenedCell) {
		file2 << _ce * 2 << " 0.00000 0.00000" << '\n';
		file2 << "0.00000 " << _ce * 2 << " 0.00000" << '\n';
		file2 << "0.00000 0.00000 " << _ce * 2 << '\n';
	}
	for (int i = 0; i < vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		file2 << vp[i].name << ' ';
	}
	file2 << '\n';
	for (int i = 0; i < vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		file2 << vp[i].points.size() << ' ';
	}
	file2 << '\n';
	file2 << "Cartesian\n";

	shift = Point();
	check = 0;
	for (int i = 0; i < vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		for (int j = 0; j < vp[i].points.size(); j++) {
			//vp[i].points[j] = mat.Transform(vp[i].points[j]);
			shift += vp[i].points[j];
			check++;
		}
	}
	shift = shift * (1.0 / check);
	shift -= Point(_ce, _ce, _ce);
	file2 << setprecision(6);
	file2.setf(ios::fixed, ios::floatfield);
	for (int i = 0; i < vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		for (int j = 0; j < vp[i].points.size(); j++) {
			vp[i].points[j] -= shift;
			file2 << vp[i].points[j].x << ' '
				<< vp[i].points[j].y << ' '
				<< vp[i].points[j].z << '\n';
		}
	}



	return 0;
}

