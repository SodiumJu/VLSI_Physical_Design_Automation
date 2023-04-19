//YR v122
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <string>

#include "main.h"

using namespace std;

auto t_begin = chrono::steady_clock::now(), t_end = chrono::steady_clock::now(), t_now = chrono::steady_clock::now();
//unordered_map<int, bool> InSetA;
vector<bool> InSetA; //recording cell is in A or not
pair<double, double> size_AB = {0, 0}; //area size for A B set
pair<int, int> cnt_AB = {0, 0}; //cell count size for A B set

bool better_cut=false; //if true, it will sort partial range of the gain vector
//just like the effect in the update bucketlist data structure
//but it will take more time to execute the program
//so the default = false

void reading_cells_nets(ifstream &cell_files, ifstream &net_files, vector<net*> &N,
 vector<cell*> &C, vector<vector<net*>> &Cell_Nets){
    string tmp;
    net *n = nullptr;
    cell *c = nullptr;
    int id = 0;
    int tmp_size=0;
    int tmp_cell_size = 0;

    while (cell_files >> tmp)
    {
        if (tmp[0] == 'c') 
        {   
            c = new cell();
            c->id = stoi(tmp.substr(1,tmp.length()-1));
        }
        else
        {
            tmp_size = stoi(tmp);
            c->asize = c->asize+tmp_size;
            cell_files >> tmp;
            tmp_size = stoi(tmp);
            c->bsize = c->bsize+tmp_size;
            C.push_back(c);
            c = nullptr;
        }
    }
    Cell_Nets.resize(C.size());
    InSetA.resize(C.size());
    sort(C.begin(), C.end(), Cell_cmp);
    
    while (net_files >> tmp)
    {
        if (tmp[0] == 'N') n = new net();
        else if (tmp[0] == 'n') n->id = stoi(tmp.substr(1,tmp.length()-1));
        else if (tmp[0] == 'c') 
        {
            id = stoi(tmp.substr(1,tmp.length()-1));
            n->cells.push_back(C[id-1]);
            tmp_cell_size = tmp_cell_size + C[id-1]->asize;
            Cell_Nets[id-1].push_back(n);
        }
        else if (tmp[0] == '}')
        {
            n->cell_size = tmp_cell_size;
            tmp_cell_size = 0;
            N.push_back(n);
            n = nullptr;
        }
    }
    sort(N.begin(), N.end(), Net_cmp);
}

void initiaPartition(vector<net*> &N, vector<cell*> &C, vector<vector<net*>> &Cell_Nets){
	//put all cells in B set first
    for(auto c: C){
        c->inA = 0; //in B set
        InSetA[c->id-1]=0;
        size_AB.second += c->bsize;
        cnt_AB.second += 1;

        for(auto n: Cell_Nets[c->id-1]){
            n->B++;
            n->bcells.push_back(c);
        }
    }

    for(auto n: N){
        if(isBalance(size_AB.first,size_AB.second)){
            break;
        }
        for(auto c: n->cells){
            if(!c->inA){
                c->inA = 1;
                InSetA[c->id-1]=1;

                size_AB.second -= c->bsize;
                cnt_AB.second -= 1;

                size_AB.first += c->asize;
                cnt_AB.first += 1;
                
                for(auto cn: Cell_Nets[c->id-1]){
                    cn->B--;
                    cn->A++;
                    cn->bcells.erase(find(cn->bcells.begin(), cn->bcells.end(), c));
                    cn->acells.push_back(c);
                }
            }
        }
    }
    

}

int get_cut_size(vector<net*> &N){
	int cutsize=0;
	for(auto n: N){
        if(n->A > 0 && n->B > 0)cutsize++;
    }
    return cutsize;
}

void Initial_gain(vector<vector<net*>> &Cell_Nets, vector<cell*> &C){
    for(auto c: C){
        for(auto n: Cell_Nets[c->id-1]){
            if(c->inA){
                if(n->A==1)c->gain++;
                if(n->B==0)c->gain--;
            }
            else{
                if(n->B==1)c->gain++;
                if(n->A==0)c->gain--;
            }
        }
    }
}
void move_gain_vec(vector<cell*> &gain_vector, bool add, cell* c){
    if(better_cut){
        vector<cell*>::iterator it;
        vector<cell*>::iterator jt;
        it=find(gain_vector.begin(), gain_vector.end(),c);
        jt=it;
        if(add){      
            for(jt; jt != gain_vector.end(); jt++){
                if( (*jt)->gain > (*it)->gain ){
                    break;
                }
            }
            sort(it,jt,gain_cmp);
        }else{
            for (jt; jt != gain_vector.begin(); jt--){
                if( (*jt)->gain < (*it)->gain ){
                    break;
                }
            }
            sort(jt,it,gain_cmp);
        } 
    }
}

void Update_gain(cell* target, vector<net*> target_net, vector<cell*> &gain_vector){
    target->lock = 1;
    if(target->inA){
        for(auto n: target_net){
            if(n->lock[0] == 1 && n->lock[1] == 1)continue; //if there are 2 cells locked in AB
                
            if(n->B == 0){
                for(auto c: n->acells){
                    if(c->lock)continue;
                    c->gain++;
                    move_gain_vec(gain_vector, 1, c);
                }
            }
            else if(n->B == 1){
                if(!n->bcells[0]->lock){
                    n->bcells[0]->gain--;
                    move_gain_vec(gain_vector, 0, n->bcells[0]);
                }
            }
            n->bcells.push_back(target);
            n->acells.erase(find(n->acells.begin(), n->acells.end(), target));
            n->A--;
            n->B++;
            n->lock[1] = 1;
        }

        InSetA[target->id-1]=0; //inB

        size_AB.first -= target->asize;
        cnt_AB.first -= 1;

        size_AB.second += target->bsize;
        cnt_AB.second += 1;
        for(auto n: target_net){
            if(n->lock[0] == 1 && n->lock[1] == 1)continue;
            if(n->A == 0){
                for(auto c: n->bcells){
                    if(c->lock)continue;
                    c->gain--;
                    move_gain_vec(gain_vector, 0, c);
                }
            }
            else if(n->A == 1){
                if(n->acells[0]->lock)continue;
                n->acells[0]->gain++;
                move_gain_vec(gain_vector, 1, n->acells[0]);
            }
        }
    }
    else{
        for(auto n: target_net){
            if(n->lock[0] == 1 && n->lock[1] == 1)continue;
            if(n->A == 0){
                for(auto c: n->bcells){
                    if(c->lock)continue;
                    c->gain++;
                    move_gain_vec(gain_vector, 1, c);
                }
            }
            else if(n->A == 1){
                if(!n->acells[0]->lock){
                    n->acells[0]->gain--;
                    move_gain_vec(gain_vector, 0, n->acells[0]);
                }
            }
            n->acells.push_back(target);
            n->bcells.erase(find(n->bcells.begin(), n->bcells.end(), target));
            n->B--;
            n->A++;
            n->lock[0] = 1;
        }

        InSetA[target->id-1]=1; //erase in B
        size_AB.second -= target->bsize;
        cnt_AB.second -= 1;
        size_AB.first += target->asize;
        cnt_AB.first += 1;

        for(auto n: target_net){
            if(n->lock[0] == 1 && n->lock[1] == 1)continue;
            if(n->B == 0){
                for(auto c: n->acells){
                    if(c->lock)continue;
                    c->gain--;
                    move_gain_vec(gain_vector, 0, c);
                }
            }
            else if(n->B == 1){
                if(n->bcells[0]->lock)continue;
                n->bcells[0]->gain++;
                move_gain_vec(gain_vector, 1, n->bcells[0]);
            }
        }
    }
}

void move_cell(vector<bool> &Ans_InSetA, pair<int, int> &Ans_cnt_AB, cell* target_cell, int &max_gain, vector<vector<net*>> &Cell_Nets,
    vector<cell*> &gain_vector, int &bestG, pair<double, double> &Ans_size_AB){
    max_gain += target_cell->gain;
    Update_gain(target_cell, Cell_Nets[target_cell->id-1], gain_vector);
    //change gain vector
    if(max_gain > bestG){
        bestG = max_gain;
        Ans_InSetA = InSetA;
        Ans_size_AB = size_AB;
        Ans_cnt_AB = cnt_AB;
    }
}

int FM(vector<net*> &N, vector<cell*> &C, vector<vector<net*>> &Cell_Nets, vector<cell*> &gain_vector, int cutsize,
    vector<bool> &Ans_InSetA, pair<double, double> &Ans_size_AB, pair<int, int> &Ans_cnt_AB){
	int max_gain = 0;
	int bestG = 0;

	do{
        if(cutsize > 394528){
            t_now = chrono::steady_clock::now();
            double computation = chrono::duration_cast<chrono::microseconds>(t_now - t_begin).count();
            //cout << computation << endl;
            if(computation > 290000000.0){
                //cout << "time out~" << computation << endl;
                break;
            }
        }
        
        int i;
        for(i = 0; i < gain_vector.size(); i++){
            if(gain_vector[i]->inA){
            	if(isBalance(size_AB.first-gain_vector[i]->asize,size_AB.second+gain_vector[i]->bsize)){
                    cell* target_cell=gain_vector[i];
                    gain_vector.erase(gain_vector.begin()+i);
                    move_cell(Ans_InSetA, Ans_cnt_AB, target_cell, max_gain, Cell_Nets, gain_vector, bestG, Ans_size_AB);
                    break;
                }
            }
            else{
            	if(isBalance(size_AB.second-gain_vector[i]->bsize, size_AB.first+gain_vector[i]->asize)){
                    cell* target_cell=gain_vector[i];
                    gain_vector.erase(gain_vector.begin()+i);
                    move_cell(Ans_InSetA, Ans_cnt_AB, target_cell, max_gain, Cell_Nets, gain_vector, bestG, Ans_size_AB);
                    break;
                }
            }
        }
        if(i == gain_vector.size())break;
    }while(max_gain > 0);
    //}while(true);
    cutsize=cutsize - bestG;
    return cutsize;
}

void savefile(ofstream &final_cut, vector<bool> &Ans_InSetA,
 vector<cell*> &C, int cutsize, pair<int, int> Ans_cnt_AB)
{
    final_cut << "cut_size " << cutsize << '\n';
    final_cut << "A " << Ans_cnt_AB.first << '\n';

    for (int i=0; i<Ans_InSetA.size(); i++){
        if(Ans_InSetA[i]==1) final_cut << 'c' << i+1 << '\n';
    }
    final_cut << "B " << Ans_cnt_AB.second << '\n';

    for (int i=0; i<Ans_InSetA.size(); i++){
        if(Ans_InSetA[i]==0) final_cut << 'c' << i+1 << '\n';
    }
    final_cut.close();
}

void print_log(bool &msg, int &initcutSize, int &cutSize,
    double &io_time, double &computation_time){
    if (msg)
    {
        cout << "------------------ Result ------------------\n";
        cout << "  Initial CutSize: " << initcutSize << '\n';
        cout << "  Final CutSize: " << cutSize << '\n';
        cout << "  I/O time: " << io_time << " microseconds\n";
        cout << "  Computation time: " << computation_time << " microseconds\n";
        cout << "  Execution time: " << io_time+computation_time << " microseconds\n";
        cout << "--------------------------------------------\n";
    }
}

int main(int argc, char *argv[]){
	bool msg=false;
	t_begin = chrono::steady_clock::now(); //time start

	ifstream cells_file, nets_file;
    ofstream final_file;
    vector<vector<net*>> Cell_Nets; //recording Nets vector for each cell
    read_arg(argc, argv, cells_file, nets_file, final_file, msg);
    vector<net*> N; vector<cell*> C, gain_vector;
    reading_cells_nets(cells_file, nets_file, N, C, Cell_Nets);

    vector<bool> Ans_InSetA;
	pair<double, double> Ans_size_AB = {0, 0};
	pair<int, int> Ans_cnt_AB = {0, 0};

    int cutsize=0, initcutSize=0;

    t_end = chrono::steady_clock::now();
    double io_time = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();

    t_begin = chrono::steady_clock::now();
    initiaPartition(N, C, Cell_Nets);
    initcutSize=get_cut_size(N);
    cutsize=initcutSize;
    Initial_gain(Cell_Nets, C);
    gain_vector = C;
    sort(gain_vector.begin(), gain_vector.end(), gain_cmp);

    cutsize=FM(N, C, Cell_Nets, gain_vector, cutsize, Ans_InSetA, Ans_size_AB, Ans_cnt_AB);
    t_end = chrono::steady_clock::now();
    double computation_time = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();

    t_begin = chrono::steady_clock::now();
    savefile(final_file, Ans_InSetA, C, cutsize, Ans_cnt_AB);
    t_end = chrono::steady_clock::now();
    io_time += chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();
    print_log(msg, initcutSize, cutsize, io_time, computation_time);
    return 0;
}
