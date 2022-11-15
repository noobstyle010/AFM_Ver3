#include <vector>
#include <stdio.h>
#include <random>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <unistd.h>

using namespace std;

/*
    単位について
    m : 10 ** -26 (Kg)
    eV : 1.609 * 10 ** -19 (J)
    d : 10 ** -10 (m) で無次元化
*/
double KB = 1.38062 * 1e-23; // [J/K]
double EV = 1.609 * 1e-19; // [J]
double D = 1e-10; // [m]

/*
    Lennard-Jones Parameterについて以下の近似式を利用
    SigmaA-B = (SigmaA-A + sigmaB-B) / 2
    EpsilonA-B = (EpsilonA-A * EpsilonB-B) ** 0.5
    各パラメーターは
    UFF, a full perodic table force field for molecular mechanics and molecular dynaics simulations
    より
*/

// Carbon
double C_SIGMA = 3.43085; // Å
double C_EPSILON = 52.83807; // K

// Oxygen
double O_SIGMA = 3.16438;
double O_EPSILON = 30.19318;

// Au
double AU_SIGMA = 2.93373; 
double AU_EPSILON = 19.62557;

/*
長い距離はカットする
*/
double CUT_OFF = 15.0;


struct Atom{
    double x;
    double y;
    double z;
    int number;
};

class Atoms{
    public:
    vector<Atom> atoms;
    
    Atoms(string model);
    void add_atom(double x, double y, double z, int number);

    

    void print(){
        for(const auto&atom : atoms){
            cout<<atom.x<<" "<<atom.y<<" "<<atom.z<<" "<<atom.number<<endl;
        }
    }
};

void Atoms::add_atom(double x, double y, double z, int number){
    atoms.push_back({x, y, z, number});
}

Atoms::Atoms(string model){
    if(model == "surface0"){
        add_atom(0, 0, 0, 79);
    }else if(model == "surface1"){
        add_atom(23, 23, 0, 79);
        add_atom(23, 23, 0, 79);
        add_atom(23, 23, 0, 79);
        add_atom(23, 23, 0, 79);
    }
};



/*
具体例
　　プローブ　試料　モード　曲率　高さ　x　y　長さ 
    .\main.exe probe0 surface0 scan1 70 3 10 20 5
*/
int main(int argc, char* argv[]){
    // main.exe Probe_Model Surface_Model Scan_Type その他
    string surface_model, probe_model, scan_type;
    double probe_radius, scan_height, scan_x, scan_y;
    int molecule_height;

    surface_model = argv[1];
    probe_model = argv[2];
    scan_type = argv[3];
    probe_radius = stod(argv[4]);
    scan_height = stod(argv[5]);
    scan_x = stod(argv[6]);
    scan_y = stod(argv[7]);
    molecule_height = stoi(argv[8]);

    // 試料の作成
    Atoms surface_atom(surface_model);

    //テスト
    surface_atom.add_atom(0, 0, 0, 10);
    surface_atom.print();
    return 0;
}

