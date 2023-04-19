#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <string>
using namespace std;

bool isBalance (double Fsize, double Tsize){
    return abs(Fsize-Tsize)*10<(Fsize+Tsize);
}

struct cell{
    bool lock = 0, inA;
    int id, gain = 0;
    double asize = 0, bsize = 0;
};

struct net{
    bool lock[2] = {0,0};
    int id, A = 0, B = 0;
    vector<cell*> cells, acells, bcells;
    int cell_size = 0;
};




bool Net_cmp(net* T, net* F){return T->cell_size < F->cell_size;}

bool Cell_cmp(cell* T, cell* F){return T->id < F->id;}

bool gain_cmp(cell* T, cell* F){return T->gain > F->gain;}


void read_arg(int argc, char *argv[], ifstream &cells_file, ifstream &nets_file, ofstream &return_cut, bool &message){
    cells_file.open(argv[1]);
    if(argc<3){
        cout << "Need 3~4 Args. Error!" << '\n';
        exit(0);
    }

    if (!cells_file.is_open())
    {
        cout << "Input Cells file Path Error" << '\n';
        exit(0);
    }

    nets_file.open(argv[2]);
    if (!nets_file.is_open())
    {
        cout << "Input Nets file Path Error" << '\n';
        exit(0);
    }

    return_cut.open(argv[3]);
    if (!return_cut.is_open())
    {
        cout << "Output Path Error" << '\n';
        exit(0);
    }
}


