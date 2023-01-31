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
    LJポテンシャルの形がやや異なるので変形が必要
*/
double LJ_CONVERTER = pow(2,1/6);
// Carbon
double C_SIGMA = 3.851/LJ_CONVERTER; // Å
double C_EPSILON = 0.105; // kcal/mol

// Oxygen
double O_SIGMA = 3.500/LJ_CONVERTER;
double O_EPSILON = 0.060;

// Auの結晶とみなす
double AU_SIGMA = 2.951; 
double AU_EPSILON = 5.29;

// Br
double BR_SIGMA = 4.189/LJ_CONVERTER;
double BR_EPSILON = 0.251;


// S 
double S_SIGMA = 4.035/LJ_CONVERTER;
double S_EPSILON = 0.274;

// H
double H_SIGMA = 2.886/LJ_CONVERTER;
double H_EPSILON = 0.044;


// SIGMAを原子の直径として見ることができる??
map<int, pair<double, double>> LUT = {
    {79, {AU_SIGMA, AU_EPSILON}},
    {8, {O_SIGMA, O_EPSILON}},
    {12, {C_SIGMA, C_EPSILON}},
    {35, {BR_SIGMA, BR_EPSILON}},
    {16, {S_SIGMA, S_EPSILON}},
    {1, {H_SIGMA, H_EPSILON}},
};

/*
長い距離はカットする
*/
double CUT_OFF = 20.0;

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

vector<double> NormV(double a1, double a2, double a3, double b1, double b2, double b3){
    double c1 = a2 * b3 - a3 * b2;
    double c2 = a3 * b1 - a1 * b3;
    double c3 = a1 * b2 - a2 * b1;
    double N = pow(pow(c1, 2) + pow(c2, 2) + pow(c3, 2), 0.5);
    vector<double> norm{c1 / N, c2 / N, c3 / N};
    return norm;
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
    double F = -48 * EPSILON * (pow(SIGMA / D1,12) - 0.5*pow(SIGMA / D1, 6)) * abs(s.z - p.z) / pow(D1, 2);
    return F;
};

class Atoms{
    public:
        vector<Atom> atoms;
        Atoms(string model, double probe_radius, double molecule_height, double molecule_distance);
        void add_atom(double x, double y, double z, int number);
        void face_centered_cubic(double x, double y, double z, int number);
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

void Atoms::face_centered_cubic(double x, double y, double z, int number){
    // FOR GOLD
    float unit_cell_length = 4.07864;
    int lx = x / unit_cell_length;
    int ly = y / unit_cell_length;
    int lz = z / unit_cell_length;
    for(int ix = 0; ix<= lx; ix++){
        for(int iy = 0; iy<=ly; iy++){
            for(int iz = 0; iz<=lz; iz++){
                add_atom(unit_cell_length * ix, unit_cell_length * iy , unit_cell_length * iz, number);
                add_atom(unit_cell_length * (double)(ix + 0.5), unit_cell_length * (double)(iy + 0.5) , unit_cell_length * iz, number);
                add_atom(unit_cell_length * (double)(ix + 0.5), unit_cell_length * iy , unit_cell_length * (double)(iz + 0.5), number);
                add_atom(unit_cell_length * ix, unit_cell_length * (double)(iy + 0.5) , unit_cell_length * (double)(iz + 0.5), number);
            }
        }
    }
}

void Atoms::cut_out_hemishphere(double probe_radius){
    vector<Atom> new_atoms(0);
    double under_limit = min<double>({0, CUT_OFF - probe_radius});
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
    if(model == "SingleGold"){
        add_atom(0, 0, 0, 79);
    // four atoms
    }else if(model == "SingleSi"){
        add_atom(0, 0, 0, 14);
    }else if(model == "FourGold"){
        add_atom(23, 23, 0, 79);
        add_atom(23, 27, 0, 79);
        add_atom(27, 23, 0, 79);
        add_atom(27, 27, 0, 79);
    // floating atom
    }else if(model == "GoldPlane"){
        face_centered_cubic(500,500,30,79);
        for(auto atom : atoms){
            atom.x -= 250;
            atom.y -= 250;
            atom.z -= 30;
        }
    }else if(model == "SingleGold"){
        add_atom(0, 0, 0, 79);
    // hemishphere gold atoms
    }else if(model == "GoldTip"){
        // Create Big Cubic of Atoms x:0~2*r+20,y;:0~2*r+20,,z:0~2*r+20
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20,79);
        cut_out_hemishphere(probe_radius);
    // hemishphere gold atoms + line gold atoms
    }else if(model == "FunctionalizedTip"){
// Create Big Cubic of Atoms x:0~2*r+20,y;:0~2*r+20,,z:0~2*r+20
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20,79);
        cut_out_hemishphere(probe_radius);
        make_needle(probe_radius, molecule_height, molecule_distance);
        /*
        まだだよ
        */
    // hemishphere gold atoms + triangle gold atoms
    }else if(model == "probe3"){
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20,79);
        cut_out_hemishphere(probe_radius);
        add_atom(0,0,-probe_radius-1,79);
        /*
        まだだよ
        */
    // hemishphere gold atoms + CO tip
    }else if(model == "COTip"){
        face_centered_cubic(2 * probe_radius + 20, 2 * probe_radius + 20, 2 * probe_radius + 20,79);
        cut_out_hemishphere(probe_radius);
        double C_O = (C_SIGMA + O_SIGMA) / 2;
        double C_AU = (C_SIGMA + AU_SIGMA) / 2;
        add_atom(0, 0, -probe_radius - C_AU, 12);
        add_atom(0, 0, -probe_radius - C_AU - C_O, 16);
        // Cが根本,Oが先っぽ
    }else if(model == "Tripod"){
        // Tripodを表現
        Atom BaseVirtualPoint{(-3.479000+2.527000+0.049000)/3,(-0.933000+3.306000-0.839000)/3,(-1.885500-1.647500+3.271500)/3, 1};
        vector<Atom> Tripod{
            
        }
        add_atom(1.814, -3.664, -1.2685, 35);
        add_atom(-3.479000,-0.933000,-1.885500, 16); // これを起点としたときS-S-S平面の法線ベクトルが[-0.547322 0.75522 0.360667]
        add_atom(2.527000, 3.306000,-1.647500, 16);  //
        add_atom(0.049000,-0.839000, 3.271500, 16);  //
        // double a1 = 0.049000 - (-3.479000);
        // double a2 = -0.839000 - (-0.933000);
        // double a3 = 3.271500 - (-1.885500);
        // double b1 = 2.527000 - (-3.479000);
        // double b2 = 3.306000 - (-0.933000);
        // double b3 = -1.647500 - (-1.885500);
        // vector<double> WOP = NormV(a1, a2, a3, b1, b2, b3);
        // cout<<WOP[0] <<" "<<WOP[1]<<" "<<WOP[2]<<endl;
        add_atom( 1.780000,  0.601000, -1.278500, 12);
        add_atom( 1.357000, -1.601000,  0.565500, 12);
        add_atom(-0.229000, -1.609000, -1.723500, 12);
        add_atom( 1.881000,  0.815000,  0.255500, 12);
        add_atom(-0.537000,  0.071000,  0.044500, 12);
        add_atom( 1.239000, -1.794000, -0.958500, 12);
        add_atom( 0.293000,  0.802000, -2.025500, 12);
        add_atom( 2.754000,  1.516000, -1.855500, 12);
        add_atom(-2.122000,  0.109000, -2.285500, 12);
        add_atom( 2.199000, -0.849000, -1.504500, 12);
        add_atom( 1.105000,  0.129000,  2.342500, 12);
        add_atom( 0.936000, -0.153000,  0.813500, 12);
        add_atom(-0.666000, -0.165000, -1.491500, 12);
        add_atom(-3.702000,-0.616000,-0.794500,1);
        add_atom( 1.496000, 3.664000,-2.398500,1);
        add_atom(-1.075000,-0.362000, 3.089500,1);
        add_atom( 2.308000,-1.765000, 1.052500,1);
        add_atom( 0.751000,-2.239000, 0.913500,1);
        add_atom(-0.836000,-2.250000,-1.381500,1);
        add_atom(-0.297000,-1.773000,-2.710500,1);
        add_atom( 1.624000, 1.748000, 0.421500,1);
        add_atom( 2.835000, 0.673000, 0.743500,1);
        add_atom(-0.818000, 0.997000, 0.199500,1);
        add_atom(-1.153000,-0.548000, 0.400500,1);
        add_atom( 0.011000, 1.734000,-1.893500,1);
        add_atom( 0.221000, 0.655000,-3.017500,1);
        add_atom( 2.647000, 1.325000,-2.844500,1);
        add_atom( 3.702000, 1.280000,-1.413500,1);
        add_atom(-2.128000, 0.011000,-3.271500,1);
        add_atom(-2.345000, 1.053000,-2.129500,1);
        add_atom( 2.142000,-1.013000,-2.490500,1);
        add_atom( 3.153000,-1.005000,-1.019500,1);
        add_atom( 0.911000, 1.089000, 2.457500,1);
        add_atom( 2.073000,-0.025000, 2.778500,1);
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
    if(scan_type == "Scan"){
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
        double helper_x = probe[0].x;
        double helper_y = probe[0].y;
        for(int i=0;i<ny;i++){
            for(int j=0;j<nx;j++){
                double local_atomic_force = 0;
                for(Atom s: surface){
                    for(Atom &p: probe){
                        p.x += dx * j - helper_x;
                        p.y += dy * i - helper_y;
                        p.z += helper_z;
                        if(atomic_d1(s, p) > CUT_OFF){
                            p.x -= dx * j - helper_x;
                            p.y -= dy * i - helper_y;
                            p.z -= helper_z;
                            break;
                        }
                        local_atomic_force += LennardJonesForce(s, p);
                        p.x -= dx * j - helper_x;
                        p.y -= dy * i - helper_y;
                        p.z -= helper_z;
                    }
                }
                cout<<local_atomic_force<<" ";
            }
            cout<<endl;
        }
    // scan along z
    // .\main.exe SingleGold GoldPlane ForceCurve 70 3 0 0 3 5.7
    }else if(scan_type == "ForceCurve"){
        vector<Atom> surface = surface_atom.atoms;
        vector<Atom> probe = probe_atom.atoms;
        int nz = 100;
        int z_off_set = 2;
        double dz = 0.05;
        // sort by z
        sort(surface.begin(), surface.end());
        sort(probe.begin(), probe.end());
        double probe_minimum_z = probe[0].z;
        double surface_max_z = surface[surface.size()-1].z;
        // z方向の補正用
        double helper_z = surface_max_z - probe_minimum_z + z_off_set;
        double helper_x = probe[0].x;
        double helper_y = probe[0].y;
        for(int i=0;i<nz;i++){
            double local_atomic_force = 0;
             for(Atom s: surface){
                    for(Atom &p: probe){
                        p.x += scan_x - helper_x;
                        p.y += scan_y - helper_y;
                        p.z += helper_z + i*dz;
                        if(atomic_d1(s, p) > CUT_OFF){
                            p.x -= scan_x - helper_x;
                            p.y -= scan_y - helper_y;
                            p.z -= helper_z + i*dz;
                            break;
                        }
                        local_atomic_force += LennardJonesForce(s, p);
                        p.x -= scan_x - helper_x;
                        p.y -= scan_y - helper_y;
                        p.z -= helper_z + i*dz;
                    }
             }
             cout<<local_atomic_force<<endl;
        }
    }else if(scan_type == "Model"){
        vector<Atom> surface = surface_atom.atoms;
        vector<Atom> probe = probe_atom.atoms;
        int z_off_set = 2;
        double dz = 0.05;
        // sort by z
        sort(surface.begin(), surface.end());
        sort(probe.begin(), probe.end());
        double probe_minimum_z = probe[0].z;
        double surface_max_z = surface[surface.size()-1].z;
        // z方向の補正用
        double helper_z = surface_max_z - probe_minimum_z + scan_height;
        double helper_x = probe[0].x;
        double helper_y = probe[0].y;
        cout<<surface.size()<<endl;
        cout<<probe.size()<<endl;
        cout<<probe_radius<<endl;
        for(Atom p: probe){
            cout<<p.x-helper_x<<" "<<p.y-helper_y<<" "<<p.z+helper_z<<" "<<p.number<<endl;
        }
        for(Atom s: surface){
            cout<<s.x - 25<<" "<<s.y-25<<" "<<s.z<<" "<<s.number<<endl;
        }
        
    }
    //テスト
    return 0;
}

