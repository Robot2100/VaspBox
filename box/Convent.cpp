// Convent.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "Classes.h"

int old_main(int argn, char * argv[])
{
	char buf[MAX_LINE];
	vector<VecPoints> vp;
	Param param;
	param.Take(argn, argv);
	ifstream file(param.filename);
	Matrix mat;
	bool defenedCell = false;
	file >> buf;
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
		file.open(param.filename);
	}


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
	if (defenedCell) {
		file2 << mat.U[0] << ' ' << mat.U[3] << ' ' << mat.U[4] << '\n';
		file2 << "0.000 "  << mat.U[1] << ' ' << mat.U[5] << '\n';
		file2 << "0.000 0.000 " << mat.U[2] << '\n';
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
	for (int i = 0; i < vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		for (int j = 0; j < vp[i].points.size(); j++) {
			Point temp = mat.Transform(vp[i].points[j]);
			file2 << temp.x << ' '
				<< temp.y << ' '
				<< temp.z << '\n';
		}
	}
	

    return 0;
}

