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
        Atoms(string model, double probe_radius, double molecule_height);
        void add_atom(double x, double y, double z, int number);
        void face_centered_cubic(double x, double y, double z);
        void cut_out_hemishphere(double r);
        void print();
};
void Atoms::print(){
    for(const auto&atom : atoms){
        cout<<atom.x<<" "<<atom.y<<" "<<atom.z<<" "<<atom.number<<endl;
    }
}
void Atoms::add_atom(double x, double y, double z, int number){
    atoms.push_back({x, y, z, number});
}

void Atoms::face_centered_cubic(double x, double y, double z){
    // FOR GOLD
    float unit_cell_length = 4.07864;
    int lx = x / unit_cell_length;
    int ly = y / unit_cell_length;
    int lz = z / unit_cell_length;
    for(int ix = 0; ix<= lx; ix++){
        for(int iy = 0; iy<=ly; iy++){
            for(int iz = 0; iz<=lz; iz++){
                add_atom(unit_cell_length * ix, unit_cell_length * iy , unit_cell_length * iz, 79);
                add_atom(unit_cell_length * (double)(ix + 0.5), unit_cell_length * (double)(iy + 0.5) , unit_cell_length * iz, 79);
                add_atom(unit_cell_length * (double)(ix + 0.5), unit_cell_length * iy , unit_cell_length * (double)(iz + 0.5), 79);
                add_atom(unit_cell_length * ix, unit_cell_length * (double)(iy + 0.5) , unit_cell_length * (double)(iz + 0.5), 79);
            }
        }
    }
}

void Atoms::cut_out_hemishphere(double probe_radius){
    vector<Atom> new_atoms(0);
    double under_limit = min<double>({0, 15 - probe_radius});
    for(auto const atom : atoms){
        double x = atom.x - probe_radius;
        double y = atom.y - probe_radius;
        double z = atom.z - probe_radius;
        if(x*x + y*y + z*z < probe_radius * probe_radius && z <= under_limit){
            Atom tmp = {x, y, z, atom.number};                
            new_atoms.push_back(tmp);
        }
    }
    atoms = new_atoms;
}

Atoms::Atoms(string model, double probe_radius = 70, double molecule_height = 3){
    // single atom
    if(model == "surface0"){
        add_atom(0, 0, 0, 79);
    // four atoms
    }else if(model == "surface1"){
        add_atom(23, 23, 0, 79);
        add_atom(23, 27, 0, 79);
        add_atom(27, 23, 0, 79);
        add_atom(27, 27, 0, 79);
    // floating atom
    }else if(model == "probe0"){
        add_atom(0, 0, 0, 79);
    // hemishphere gold atoms
    }else if(model == "probe1"){
        // Create Big Cubic of Atoms x:0~2*r+20,y;:0~2*r+20,,z:0~2*r+20
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20);
        cut_out_hemishphere(probe_radius);
    // hemishphere gold atoms + line gold atoms
    }else if(model == "probe2"){
// Create Big Cubic of Atoms x:0~2*r+20,y;:0~2*r+20,,z:0~2*r+20
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20);
        cut_out_hemishphere(probe_radius);
        /*
        まだだよ
        */
    // hemishphere gold atoms + triangle gold atoms
    }else if(model == "probe3"){
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20);
        cut_out_hemishphere(probe_radius);
        /*
        まだだよ
        */
    // hemishphere gold atoms + CO tip
    }else if(model == "probe4"){
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20);
        cut_out_hemishphere(probe_radius);
        /*
        まだだよ
        */
    }else{
        cout<<"MISS MODEL_NAME"<<endl;
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
    
    probe_model = argv[1];
    surface_model = argv[2];
    scan_type = argv[3];
    probe_radius = stod(argv[4]);
    scan_height = stod(argv[5]);
    scan_x = stod(argv[6]);
    scan_y = stod(argv[7]);
    molecule_height = stoi(argv[8]);

    // 試料とプローブの作成
    Atoms surface_atom(surface_model);
    Atoms probe_atom(probe_model, probe_radius, molecule_height);

    // 計算を行う
    /*
    まだだよ
    */
   
    //テスト
    surface_atom.add_atom(0, 0, 0, 10);
    surface_atom.print();
    probe_atom.print();
    return 0;
}

