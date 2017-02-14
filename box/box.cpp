// box.cpp: определ€ет точку входа дл€ консольного приложени€.
//

#include "stdafx.h"
#include "Classes.h"
using namespace std;

int main(int argn, char * argv[])
{
	string atomname("");
	bool first_atom = true;
	string filename("a.txt");
	flo _ce = 20.0;
	flo bond = 2.0;
	{
		Param param;
		param.AddParam({ "-h","-help" }, 0);
		param.AddParam({ "-f","-filename" }, 1);
		param.AddParam({ "-c","-cell" }, 2);
		param.AddParam({ "-b","-bond" }, 3);
		param.AddParam({ "-a","-atom" }, 4);
		string strParam = param.TakeAgrs(argn, argv);
		if (strParam.compare("")) 
			cout << strParam << endl;
		vector<string> vStrParam;
		int res;
		while ((res = param.ReadNextParametr(vStrParam)) != -1) {
			switch (res) {
			case 0:
				cout << "Parameters:" << endl;
				cout << "     -h or -help - show this help" << endl;
				cout << "     -f or -filename <filename> - change input filename [a.txt]" << endl;
				cout << "     -c or -cell <vector length> - change cell vector length [20.0]" << endl;
				cout << "     -b or -bond <cutoff length> - change cutoff length [2.0]" << endl;
				cout << "     -a or -atom <AtomName> - determine molecule by atom [first founded]" << endl;
				system("pause");
				exit(0);
			case 1:
				if (vStrParam.size() != 1) {
					cout << "Parameter -f/-filename has wrong number of args!" << endl;
					cout << "Program terminated!" << endl;
					system("pause");
					exit(-1);
				}
				filename = vStrParam[0];
				break;
			case 2:
				if (vStrParam.size() != 1) {
					cout << "Parameter -c/-cell has wrong number of args!" << endl;
					cout << "Program terminated!" << endl;
					system("pause");
					exit(-1);
				}
				{
					flo temp = (flo)atof(vStrParam[0].data());
					if (temp > 0) _ce = temp;
					else {
						cout << "Parameter -c/-cell has unknown argument! " << endl;
						cout << "Program terminated!" << endl;
						system("pause");
						exit(-1);
					}
				}
				break;
			case 3:
				if (vStrParam.size() != 1) {
					cout << "Parameter -b/-bond has wrong number of args!" << endl;
					cout << "Program terminated!" << endl;
					system("pause");
					exit(-1);
				}
				{
					flo temp = (flo)atof(vStrParam[0].data());
					if (temp > 0) bond = temp;
					else {
						cout << "Parameter -b/-bond has unknown argument! " << endl;
						cout << "Program terminated!" << endl;
						system("pause");
						exit(-1);
					}
				}
				break;
			case 4:
				if (vStrParam.size() != 1) {
					cout << "Parameter -a/-atom has wrong number of args!" << endl;
					cout << "Program terminated!" << endl;
					system("pause");
					exit(-1);
				}
				atomname = vStrParam[0];
				first_atom = false;
				break;
			}
		}
	}

	int _type=0, _number=0;

	cout << "Launch parameters:" << endl;
	cout << "\t filename: " << filename << endl;
	cout << "\t vector length: " << _ce << endl;
	cout << "\t cutoff length: " << bond << endl;
	if(first_atom) cout << "\t atom name: [first founded]" << endl;
	else cout << "\t atom name: " << atomname << endl;
	char buf[MAX_LINE];
	vector<VecPoints> vp;
	ifstream file(filename);
	if (!file.is_open()) exit(-1);
	Matrix mat;
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
	}
	else {
		file.close();
		cout << "Error! Can't read CELL line." << endl;
		exit(-1);
	}
	bool catched = false;
	// ÷икл считывани€ атомов
	while (!file.eof())
	{
		bool monitor = false;
		file >> buf;
		if (buf[0] == '=')
			continue;
		int i = 0;
		if (!first_atom && !catched) {
			if (atomname.compare(buf) == 0) {
				monitor = true;
				catched = true;
			}
		}
		for (; i < MAX_LINE; i++) {
			if ((buf[i]<'A' || buf[i]>'Z') && (buf[i]<'a' || buf[i]>'z')) {
				buf[i] = '\0';
				break;
			}
		}
		if (i > 3) 
			exit(1);
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
		for (i = 0; i < (int)vp.size(); i++) {
			if (strcmp(vp[i].name, atom) == 0) {
				break;
			}
		}
		if (vp.size() == i) vp.push_back(VecPoints(atom));
		vp[i].points.push_back(k);
		if (monitor) {
			_type = i;
			_number = vp[i].points.size()-1;
		}
	}
	if (!(catched || first_atom)) {
		cout << "Atom  is not founded!" << endl;
		cout << "Program terminated!" << endl;
		system("pause");
		exit(-1);
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

	
	file2 << _ce << " 0.00000 0.00000" << '\n';
	file2 << "0.00000 " << _ce << " 0.00000" << '\n';
	file2 << "0.00000 0.00000 " << _ce << '\n';
	for (int i = 0; i < (int)vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		file2 << vp[i].name << ' ';
	}
	file2 << '\n';
	for (int i = 0; i < (int)vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		file2 << vp[i].points.size() << ' ';
	}
	file2 << '\n';
	file2 << "Cartesian\n";

	shift = Point();
	check = 0;
	for (int i = 0; i < (int)vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		for (int j = 0; j < (int)vp[i].points.size(); j++) {
			//vp[i].points[j] = mat.Transform(vp[i].points[j]);
			shift += vp[i].points[j];
			check++;
		}
	}
	shift = shift * (1.0f / check);
	shift -= Point(_ce/2, _ce/2, _ce/2);
	file2 << setprecision(6);
	file2.setf(ios::fixed, ios::floatfield);
	for (int i = 0; i < (int)vp.size(); i++) {
		if (vp[i].name[0] == '\0') break;
		for (int j = 0; j < (int)vp[i].points.size(); j++) {
			vp[i].points[j] -= shift;
			file2 << vp[i].points[j].x << ' '
				<< vp[i].points[j].y << ' '
				<< vp[i].points[j].z << '\n';
		}
	}



	return 0;
}

