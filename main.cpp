#include <vector>
#include <stdio.h>
#include <random>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <map>

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

map<int, pair<double, double>> LUT = {
    {79, {AU_SIGMA, AU_EPSILON}},
    {16, {O_SIGMA, O_EPSILON}},
    {12, {C_SIGMA, C_EPSILON}}
};

/*
長い距離はカットする
*/
double CUT_OFF = 15.0;

int random_needle_num = 1000;


struct Atom{
    double x;
    double y;
    double z;
    int number;

    bool operator<(const Atom &another) const{
        return z < another.z;
    }
};

double atomic_d2(Atom a, Atom b){
    return pow((a.x-b.x) , 2) + pow((a.y - b.y) , 2) + pow((a.z - b.z) , 2);
};

double atomic_d1(Atom a, Atom b){
    return pow(pow((a.x-b.x) , 2) + pow((a.y - b.y) , 2) + pow((a.z - b.z) , 2), 0.5);
};

double LennardJonesForce(Atom s, Atom p, double dx = 0, double dy = 0, double dz = 0){
    p.x += dx;
    p.y += dy;
    p.z += dz;
    double SIGMA = (LUT[s.number].first + LUT[p.number].first) / 2;
    double EPSILON = pow(LUT[s.number].second * LUT[p.number].second, 0.5);
    double D1 = atomic_d1(s, p);
    double F = -48 * EPSILON * (pow(SIGMA / D1,12) - 0.5 * pow(SIGMA / D1, 6)) * abs(s.z - p.z) / pow(D1, 2);
    return F;
};

class Atoms{
    public:
        vector<Atom> atoms;
        Atoms(string model, double probe_radius, double molecule_height, double molecule_distance);
        void add_atom(double x, double y, double z, int number);
        void face_centered_cubic(double x, double y, double z);
        void cut_out_hemishphere(double r);
        void make_needle(double probe_radius, int molecule_height, double molecule_distance);
        void print();
};

void Atoms::print(){
    for(const auto&atom : atoms){
        cout<<atom.x<<" "<<atom.y<<" "<<atom.z<<" "<<atom.number<<endl;
    }
};

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

void Atoms::make_needle(double probe_radius, int molecule_height, double molecule_distance){
    srand(time(NULL));
    vector<Atom> random_points(0);
    for(int i = 0; i < random_needle_num; i++){
        double t = (double)rand()/RAND_MAX;
        double u = (double)rand()/RAND_MAX;
        t = -1 * asin(1 - 2 * t);
        u = u * 2 * M_PI - M_PI;
        double random_x = probe_radius * cos(t) * cos(u);
        double random_y = probe_radius * cos(t) * sin(u);
        double random_z = -1 * abs(probe_radius * sin(t)); 
        random_points.push_back({random_x, random_y, random_z, 79});
    }
    // 刈り取り
    vector<Atom> base_points(0);
    while(1){
        if(random_points.size()== 1){
            atoms.push_back(random_points[0]);
            break;
        }
        if(random_points.size() == 0){
            break;
        }
        Atom tmp_atom = random_points.back();
        vector<Atom> tmp_atoms(0);
        random_points.pop_back();
        for(const auto&atom : random_points){
            if(atomic_d2(atom, tmp_atom) > molecule_distance * molecule_distance){
                tmp_atoms.push_back(atom);
            }
        }
        random_points = tmp_atoms;
        base_points.push_back(tmp_atom);
    }
    for(const auto&base_point : base_points){
        for(int i = 1; i <= molecule_height; i++){
            add_atom(base_point.x * (1 + AU_SIGMA * (double)i/probe_radius), base_point.y * (1 + AU_SIGMA * (double)i/probe_radius), base_point.z * (1 + AU_SIGMA * (double)i/probe_radius), base_point.number);
        }
    }
}

Atoms::Atoms(string model, double probe_radius, double molecule_height, double molecule_distance){
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
        make_needle(probe_radius, molecule_height, molecule_distance);
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
　　プローブ　試料　モード　曲率　高さ　x　y　長さ 間隔
    .\main.exe probe0 surface0 scan1 70 3 10 20 5 10
    #針一本くらい
    .\main.exe probe1 surface0 scan1 1 3 0 0 1 10
*/
int main(int argc, char* argv[]){
    // main.exe Probe_Model Surface_Model Scan_Type その他
    string surface_model, probe_model, scan_type;
    double probe_radius, scan_height, scan_x, scan_y, molecule_distance;
    int molecule_height;
    
    probe_model = argv[1];
    surface_model = argv[2];
    scan_type = argv[3];
    probe_radius = stod(argv[4]);
    scan_height = stod(argv[5]);
    scan_x = stod(argv[6]);
    scan_y = stod(argv[7]);
    molecule_height = stoi(argv[8]);
    molecule_distance = stod(argv[9]);

    // 試料とプローブの作成
    Atoms surface_atom(surface_model, 0, 0, 0);
    Atoms probe_atom(probe_model, probe_radius, molecule_height, molecule_distance);
    // SCAN XY
    if(scan_type == "scan0"){
        vector<Atom> surface = surface_atom.atoms;
        vector<Atom> probe = probe_atom.atoms;
        int nx = 50;
        int ny = 50;
        double dx = 1;
        double dy = 1;
        // sort by z
        sort(surface.begin(), surface.end());
        sort(probe.begin(), probe.end());
        double probe_minimum_z = probe[0].z;
        double surface_max_z = surface[surface.size()-1].z;
        // z方向の補正用
        double helper_z = surface_max_z - probe_minimum_z + scan_height;
        for(int i=0;i<ny;i++){
            for(int j=0;j<nx;j++){
                double local_atomic_force = 0;
                for(Atom s: surface){
                    for(Atom p: probe){
                        p.x += dx * j;
                        p.y += dy * i;
                        p.z += helper_z;
                        if(atomic_d1(s, p) > 15){break;}
                        local_atomic_force += LennardJonesForce(s, p);
                    }
                }
                cout<<local_atomic_force<<" ";
            }
            cout<<endl;
        }
    // scan along z
    }else if(scan_type == "scan1"){

    }
    //テスト
    return 0;
}

