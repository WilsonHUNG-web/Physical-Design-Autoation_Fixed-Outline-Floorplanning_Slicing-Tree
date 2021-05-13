#include<iostream>
#include<tuple>
#include<list>
#include<sstream>
#include<cstdlib>
#include<cstring>
#include<chrono>
#include<unistd.h>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<algorithm>
#include<ctime>
#include<cstdio>
#include<string>
#include<stack>
#include<vector>
#include<omp.h>
//Parameters
#define initTcount 1000
#define weight_outline_W 50000  //dimension control W
#define weight_outline_H 50000 //dimension control H
#define Tlimit0 100
#define Tlimit1 0.0000000001
#define Tlimit2 0.000000000000001
#define HPWLweight 0.001// lumbda*HPWL
#define areaweight 1 //AREA
#define seed 691600
#define MinTemperature 0.00000000000001
#define epsilon 10e-30
#define terminatecount 100
#define probability 0.9
//#define aspectweight 0
#define kmove 6
#define V -3
#define H -4
#define none -10
#define coolrate 0.95

using namespace std;

class net {
public:
	net() {}
	int deg;
	vector<int> ListBlock, ListTermi;
};


class terminal {
public:
	terminal() {
		id = none;
	}
	int id;
	float x, y;
};

class block {
public:
	block() {}
	int area;
	bool rotate;
	int w, h;
	float x, y;
};

class node {
public:
	node() {}
	int id;
	int idx;
	int Parent = none;
	int Lchild = none, Rchild = none;
	int w = 0, h = 0;
	int newH, newW;
	int x, y;
	vector<tuple<int, int, int, int> > listWH;
};

//
int NumBlocks, NumTermi, NumNets, NumPin;
int total_area, FixOutline;
double HPWL;
int total_W, total_H;
vector<block*> dicBlock;
vector<net> dicNet;
vector<terminal*> terminals;
node* iniTree, *bestTree, *dummyTree;
//

chrono::high_resolution_clock::time_point time_record() {
	return chrono::high_resolution_clock::now();
}

double time_output(chrono::high_resolution_clock::time_point start_time, chrono::high_resolution_clock::time_point end_time) {
	return chrono::duration<double>(end_time - start_time).count();
}


void parse_block(FILE* input) {

	fscanf(input, "NumHardRectilinearBlocks : %d\n", &NumBlocks);
	fscanf(input, "NumTerminals : %d\n\n", &NumTermi);


	for (int i = 0; i < NumBlocks; i++) {
		int a;
		block* temp = new block();
		fscanf(input, " sb%d hardrectilinear 4 (%f, %f) (%d, %d) (%d, %d) (%d, %d)\n", &a, &temp->x, &temp->y, &a, &a, &temp->w, &temp->h, &a, &a);

		temp->area = temp->w * temp->h;
		total_area += temp->area;
		temp->rotate = false;
		dicBlock.push_back(temp);
	}
	

	return;
}

void parse_terminal(FILE* input) {
	terminal* temp1 = new terminal();
	terminals.push_back(temp1);
	for (int i = 0; i < NumTermi; i++) {
		terminal* temp2 = new terminal();
		fscanf(input, " p%d %f %f\n", &temp2->id, &temp2->x, &temp2->y);
		terminals.push_back(temp2);
	}
	return;
}

void parse_net(FILE* input) {
	fscanf(input, " NumNets : %d\n", &NumNets);
	fscanf(input, " NumPins : %d\n", &NumPin);

	for(int i = 0; i < NumNets; i++){
		net temp;
		string str;
		char c[40];
		fscanf(input, " NetDegree : %d\n", &temp.deg);
		for (int i = 0; i < temp.deg; i++) {
			fscanf(input, " %s\n", c);
			str = c;
			if (str[0] == 'p') {
				str.erase(0, 1);
				temp.ListTermi.push_back(stoi(str));
			}
			else {
				str.erase(0, 2);
				temp.ListBlock.push_back(stoi(str));
			}
		}
		dicNet.push_back(temp);
	}
	return;
}

vector<int> ini_NPE() {

	vector<int> NPE;
	NPE.push_back(0);
	for (int i = 1; i < NumBlocks; i++) {
		NPE.push_back(i);
		int m = rand() % 2;
		if (m)NPE.push_back(V);
		else NPE.push_back(V);
	}
	return NPE;
}


void ini_Tree(vector<int> NPE) {
	stack<int> stk;
	iniTree = new node[2 * NumBlocks - 1]();
	int index1, index2;

	for (int i = 0; i < (int)NPE.size(); i++) {
		iniTree[i].id = NPE[i];
		if (NPE[i] >= 0) {
			stk.push(i);
			if (dicBlock[NPE[i]]->w > dicBlock[NPE[i]]->h) {
				iniTree[i].w = dicBlock[NPE[i]]->w;
				iniTree[i].h = dicBlock[NPE[i]]->h;
			}
			else {
				iniTree[i].w = dicBlock[NPE[i]]->h;
				iniTree[i].h = dicBlock[NPE[i]]->w;
			}

			iniTree[i].listWH.push_back(make_tuple(iniTree[i].w, iniTree[i].h, -1, -1));
			if (iniTree[i].w != iniTree[i].h) {
				iniTree[i].listWH.push_back(make_tuple(iniTree[i].h, iniTree[i].w, -1, -1));
			}
		}
		else if (NPE[i] == V) { 
			index2 = iniTree[i].Rchild = stk.top();
			stk.pop();
			index1 = iniTree[i].Lchild = stk.top();
			stk.pop();
			stk.push(i);
			iniTree[iniTree[i].Rchild].Parent = iniTree[iniTree[i].Lchild].Parent = i;
			int j = 0, k = 0;
			while (j < (int)iniTree[index1].listWH.size() && k < (int)iniTree[index2].listWH.size()) {
				int w = 0, h = 0;
				w = get<0>(iniTree[index1].listWH[j]) + get<0>(iniTree[index2].listWH[k]);
				h = max(get<1>(iniTree[index1].listWH[j]), get<1>(iniTree[index2].listWH[k]));
				iniTree[i].listWH.push_back(make_tuple(w, h, j, k));
				if (get<1>(iniTree[index1].listWH[j]) > get<1>(iniTree[index2].listWH[k])) {
					j++;
				}
				else if (get<1>(iniTree[index1].listWH[j]) < get<1>(iniTree[index2].listWH[k])) {
					k++;
				}
				else {
					j++;
					k++;
				}
			}

		}
		else if (NPE[i] == H) { 
			index2 = iniTree[i].Rchild = stk.top();
			stk.pop();
			index1 = iniTree[i].Lchild = stk.top();
			stk.pop();
			stk.push(i);
			iniTree[iniTree[i].Rchild].Parent = iniTree[iniTree[i].Lchild].Parent = i;
			int j = (int)iniTree[index1].listWH.size()-1, k = (int)iniTree[index2].listWH.size()-1;
			while (j >=0  && k >=0 ) {
				int w = 0, h = 0;
				w = max(get<0>(iniTree[index1].listWH[j]), get<0>(iniTree[index2].listWH[k]));
				h = get<1>(iniTree[index1].listWH[j]) + get<1>(iniTree[index2].listWH[k]);
				iniTree[i].listWH.push_back(make_tuple(w, h, j, k));
				if (get<0>(iniTree[index1].listWH[j]) < get<0>(iniTree[index2].listWH[k])) {
					k--;
				}
				else if (get<0>(iniTree[index1].listWH[j]) > get<0>(iniTree[index2].listWH[k])) {
					j--;
				}
				else {
					j--;
					k--;
				}
			}
		}

	}
	return;
}

void Compute_coord(node* tree, int idx, int x, int y) {
	tree[idx].x = x;
	tree[idx].y = y;
	int lc_x, lc_y, rc_x, rc_y;
	if (tree[idx].id == V) {
		lc_x = tree[idx].x;
		lc_y = tree[idx].y;
		rc_x = tree[idx].x + tree[tree[idx].Lchild].newW;
		rc_y = tree[idx].y;
		Compute_coord(tree, tree[idx].Lchild, lc_x, lc_y);
		Compute_coord(tree, tree[idx].Rchild, rc_x, rc_y);
	}
	else if (tree[idx].id == H) {
		lc_x = tree[idx].x;
		lc_y = tree[idx].y;
		rc_x = tree[idx].x;
		rc_y = tree[idx].y + tree[tree[idx].Lchild].newH;
		Compute_coord(tree, tree[idx].Lchild, lc_x, lc_y);
		Compute_coord(tree, tree[idx].Rchild, rc_x, rc_y);
	}
	else if (tree[idx].id >= 0) { //other than V, H 
		if (dicBlock[tree[idx].id]->w != tree[idx].newW) {
			int temp;
			temp = dicBlock[tree[idx].id]->w;
			dicBlock[tree[idx].id]->w = dicBlock[tree[idx].id]->h;
			dicBlock[tree[idx].id]->h = temp;
			dicBlock[tree[idx].id]->rotate = !dicBlock[tree[idx].id]->rotate;
		}
		dicBlock[tree[idx].id]->x = x;
		dicBlock[tree[idx].id]->y = y;
	}
	return;
}

void Compute_WH(node* tree, int picked, int idx) {
	int lc_idx, rc_idx, lc_picked, rc_picked;
	tree[idx].newW = get<0>(tree[idx].listWH[picked]);
	tree[idx].newH = get<1>(tree[idx].listWH[picked]);
	lc_picked = get<2>(tree[idx].listWH[picked]);
	rc_picked = get<3>(tree[idx].listWH[picked]);
	lc_idx = tree[idx].Lchild;
	rc_idx = tree[idx].Rchild;
	if (lc_idx != none) Compute_WH(tree, lc_picked, lc_idx);
	if (rc_idx != none) Compute_WH(tree, rc_picked, rc_idx);
	return;
}

void compute_area(stack<int> bigstack, node* tree, int root_index) {
	int temp, index1, index2;
	while (!bigstack.empty()) {
		temp = bigstack.top();
		bigstack.pop();
		if (tree[temp].id == V) {
			index2 = tree[temp].Rchild;
			index1 = tree[temp].Lchild;

			tree[temp].listWH.clear();
			int j = 0, k = 0;
			while (j < (int)tree[index1].listWH.size() && k < (int)tree[index2].listWH.size()) {
				int w = 0, h = 0;
				w = get<0>(tree[index1].listWH[j]) + get<0>(tree[index2].listWH[k]);
				h = max(get<1>(tree[index1].listWH[j]), get<1>(tree[index2].listWH[k]));
				tree[temp].listWH.push_back(make_tuple(w, h, j, k));
				if (get<1>(tree[index1].listWH[j]) > get<1>(tree[index2].listWH[k])) {
					j++;
				}
				else if (get<1>(tree[index1].listWH[j]) < get<1>(tree[index2].listWH[k])) {
					k++;
				}
				else {
					j++;
					k++;
				}
			}
		}
		//else if (tree[temp].id == H) {

		//	index2 = tree[temp].Rchild;
		//	index1 = tree[temp].Lchild;
		//	//cout << "hello1" << endl;
		//	tree[temp].listWH.clear();
		//	//cout << "hello2" << endl;
		//	int j = (int)tree[index1].listWH.size()-1, k = (int)tree[index2].listWH.size()-1;
		//	//cout << "hello3" << endl;
		//	while (j >=0 && k >=0) {
		//		int w = 0, h = 0;
		//		//cout << "hello w = 0, h = 0;" << endl;
		//		w = max(get<0>(tree[index1].listWH[j]), get<0>(tree[index2].listWH[k]));
		//		//cout << "hello4" << endl;
		//		h = get<1>(tree[index1].listWH[j]) + get<1>(tree[index2].listWH[k]);
		//		//cout << "hello5" << endl;
		//		tree[temp].listWH.push_back(make_tuple(w, h, j, k));
		//		if (get<0>(tree[index1].listWH[j]) < get<0>(tree[index2].listWH[k])) {
		//		
		//			k--;
		//			//cout << "hello6 k:" << k << endl;
		//		}
		//		else if (get<0>(tree[index1].listWH[j]) > get<0>(tree[index2].listWH[k])) {
		//			
		//			j--;
		//			//cout << "hello7 j:" << j << endl;
		//		}
		//		else {
		//			j--;
		//			k--; 
		//			//cout << "hello8 j:" << j<<" k:" <<k<< endl;
		//		}
		//	}
		//}


		else if (tree[temp].id == H) {

			index2 = tree[temp].Rchild;
			index1 = tree[temp].Lchild;
			//cout << "hello1" << endl;
			tree[temp].listWH.clear();
			//cout << "hello2" << endl;
			int j = 0, k = 0;
			//cout << "hello3" << endl;
			while (j < (int)tree[index1].listWH.size() && k < (int)tree[index2].listWH.size()) {
				int w = 0, h = 0;
				//cout << "hello w = 0, h = 0;" << endl;
				w = max(get<0>(tree[index1].listWH[j]), get<0>(tree[index2].listWH[k]));
				//cout << "hello4" << endl;
				h = get<1>(tree[index1].listWH[j]) + get<1>(tree[index2].listWH[k]);
				//cout << "hello5" << endl;
				tree[temp].listWH.push_back(make_tuple(w, h, j, k));
				if (get<0>(tree[index1].listWH[j]) < get<0>(tree[index2].listWH[k])) {

					j++;
					//cout << "hello6 k:" << k << endl;
				}
				else if (get<0>(tree[index1].listWH[j]) > get<0>(tree[index2].listWH[k])) {

					k++;
					//cout << "hello7 j:" << j << endl;
				}
				else {
					j++;
					k++;
					//cout << "hello8 j:" << j<<" k:" <<k<< endl;
				}
			}
		}



		else {
			cout << tree[temp].id << " " << tree[tree[temp].Parent].id << " " << tree[temp].Lchild << " " << tree[temp].Rchild << endl;
		}
	}

	int min = 99999999;
	int w, h, picked = 0;
	for (int i = 0; i < (int)tree[root_index].listWH.size(); i++) {

		w = get<0>(tree[root_index].listWH[i]);
		h = get<1>(tree[root_index].listWH[i]);
		if (w*h < min) {
			min = w * h;
			total_W = w;
			total_H = h;
			picked = i;
		}
	}
	Compute_WH(tree, picked, root_index);
	Compute_coord(tree, root_index, 0, 0);
	return;
}

void Compute_HPWL() {

	HPWL = 0; //clean HPWL

	for (int i = 0; i < NumNets; i++) {
		int minX = 9999999, minY = 9999999, maxX = 0, maxY = 0;
		int X, Y;


		for (auto j = dicNet[i].ListBlock.begin(); j != dicNet[i].ListBlock.end(); j++) {
			X = dicBlock[*j]->x + dicBlock[*j]->w / 2;
			Y = dicBlock[*j]->y + dicBlock[*j]->h / 2;
			maxX = X > maxX ? X : maxX;
			maxY = Y > maxY ? Y : maxY;
			minX = X < minX ? X : minX;
			minY = Y < minY ? Y : minY;
		}

		for (auto j = dicNet[i].ListTermi.begin(); j != dicNet[i].ListTermi.end(); j++) {
			X = terminals[*j]->x;
			Y = terminals[*j]->y;
			maxX = X > maxX ? X : maxX;
			maxY = Y > maxY ? Y : maxY;
			minX = X < minX ? X : minX;
			minY = Y < minY ? Y : minY;
		}
		HPWL += (maxX - minX) + (maxY - minY);
	}
	return;
}

int compute_cost() {//COST = A + W + Aspecratio
	int cost = 0;

	cost += (total_H * total_W)*areaweight; //A
	

	if (total_H > FixOutline) cost += (total_H - FixOutline)* weight_outline_H;
	if (total_W > FixOutline) cost += (total_W - FixOutline) * weight_outline_W;


	//float aspect = (float)total_H / (float)total_W;
	//if (aspect > 1.02 || aspect < 0.98) cost += (total_H - total_W)*(total_H - total_W)*aspectweight;


	Compute_HPWL();  //update HPWL for cost usage
	
	cost += (int)(HPWLweight * HPWL);
	return cost;
}


vector<int> M1(vector<int> NPE, node* tree) {
	
	int min = 1, cnt = 0;
	int picked = rand() % (NumBlocks - min) + min, picked2 = picked;
	while (picked == picked2) {
		picked2 = rand() % (NumBlocks - min) + min;
	}

	int index1 = -1, index2 = -1;

	for (int i = 0; i < (int)NPE.size(); i++) {
		if (NPE[i] >= 0) {
			cnt++;
			if (cnt == picked) {
				index1 = i;
			}
			if (cnt == picked2) {
				index2 = i;
			}
		}
	}
	if (index1 == -1 || index2 == -1) cout << "idx == -1" << endl;
	if (NPE[index1] < 0 || NPE[index2] < 0) cout << "NPE[idx]<0" << endl;

	int temp = NPE[index1];
	NPE[index1] = NPE[index2];
	NPE[index2] = temp;

	temp = tree[index1].Parent;
	tree[index1].Parent = tree[index2].Parent;
	tree[index2].Parent = temp;

	node n;
	n = tree[index1];
	tree[index1] = tree[index2];
	tree[index2] = n;


	stack<int> s1, s2, bigstack;
	n = tree[index1];
	while (n.Parent != none) {
		s1.push(n.Parent);
		n = tree[n.Parent];
	}

	n = tree[index2];
	while (n.Parent != none) {
		s2.push(n.Parent);
		n = tree[n.Parent];
	}
	while (!s1.empty() && !s2.empty()) {
		if (s1.top() == s2.top()) {
			bigstack.push(s1.top());
			s1.pop();
			s2.pop();
		}
		else {
			break;
		}
	}
	while (!s1.empty()) {
		bigstack.push(s1.top());
		s1.pop();
	}
	while (!s2.empty()) {
		bigstack.push(s2.top());
		s2.pop();
	}
	compute_area(bigstack, tree, NPE.size() - 1);

	return NPE;
}

vector<int> M2(vector<int>NPE, node* tree) {
	int cnt = 0, idx = 0;
	int picked = rand() % (NumBlocks - 1) + 1;
	stack<int> stk, bigstack;
	for (int i = 0; i < (int)NPE.size(); i++) {
		if (NPE[i] < 0) {
			cnt++;
			if (cnt == picked) {
				idx = i;
				break;
			}
		}
	}
	while (NPE[idx - 1] < 0) {
		idx = idx - 1;
	}
	while (NPE[idx] == V || NPE[idx] == H) { 
		tree[idx].id = tree[idx].id == V ? H : V;
		NPE[idx] = NPE[idx] == V ? H : V;
		stk.push(idx);
		if (idx == (int)NPE.size() - 1) break;
		if (NPE[idx + 1] == V || NPE[idx + 1] == H) idx = idx + 1;
		else break;
	}
	while (tree[idx].Parent != none && (tree[tree[idx].Parent].id == V || tree[tree[idx].Parent].id == H)) {
		stk.push(tree[idx].Parent);
		idx = tree[idx].Parent;
	}

	while (!stk.empty()) {
		bigstack.push(stk.top());
		stk.pop();
	}

	compute_area(bigstack, tree, NPE.size() - 1);
	return NPE;
}


void cpyTree(node* t1, node* t2, int k) {
	for (int i = 0; i < k; i++) {
		t1[i] = t2[i];
	}
	return;
}

vector<int> M3(vector<int>NPE, node* tree) {
	int picked, NumOptr, NumOpnd, index1, index2, terminate = 0, temp;
	bool done = false;
	while (!done) {
		if (terminate++ > terminatecount) {
			return NPE;
		}
		NumOptr = NumOpnd = 0;
		picked = rand() % (NumBlocks)+1;

		for (int i = 0; i < (int)NPE.size(); i++) {
			if (NPE[i] == V || NPE[i] == H) {
				NumOptr++;
			}
			else if (NPE[i] >= 0) {
				NumOpnd++;
				if (NumOpnd == picked) {
					if (i + 1 < (2 * NumBlocks - 1) && NPE[i + 1] < 0 && (NPE[i + 1] != NPE[i - 1]) && (NumOpnd-1 > NumOptr+1 )) {
						index2 = i + 1;
						index1 = i;
				
						done = true;

						temp = NPE[index1];
						NPE[index1] = NPE[index2];
						NPE[index2] = temp;
						//cout<<"1";
						//change tree
						temp = index2; // H/V
						int flag = 0;
						while (tree[tree[temp].Parent].Lchild == temp) {
							temp = tree[temp].Parent;
							flag = 1;
						}
						temp = tree[temp].Parent;

						tree[tree[temp].Lchild].Parent = index1;
						tree[index2].Rchild = tree[index2].Lchild;
						tree[index2].Lchild = tree[temp].Lchild;
						tree[tree[index2].Rchild].Parent = index1;
						tree[temp].Lchild = index1;
						tree[index1].Parent = tree[index2].Parent;
						if (flag == 0)tree[tree[index2].Parent].Rchild = index2;
						else tree[tree[index2].Parent].Lchild = index2;
						tree[index2].Parent = temp;
						//						
												//cout<<"3";
						node n;
						n = tree[index1];
						tree[index1] = tree[index2];
						tree[index2] = n;
					}
					else if (i - 1 > 0 && NPE[i - 1] < 0 && (NPE[i + 1] != NPE[i - 1])) {
							//cout<<2;
							index1 = i;
							index2 = i - 1;
							done = true;
							//change order in NPE
							temp = NPE[index1];
							NPE[index1] = NPE[index2];
							NPE[index2] = temp;
							//change H/V's left child's Parent to its Parent
							temp = tree[index2].Lchild;
							tree[tree[index2].Lchild].Parent = tree[index2].Parent;
							tree[tree[index2].Parent].Lchild = tree[index2].Lchild;
							//change H/V's left child to right child
							tree[tree[index2].Rchild].Parent = index1;
							tree[index2].Lchild = tree[index2].Rchild;
							temp = tree[index1].Parent;

							tree[index1].Parent = tree[index2].Parent;
							tree[index2].Parent = temp;
							tree[index2].Rchild = index2;
							tree[index1].Parent = index1;
							//change order in tree
							node n;
							n = tree[index1];
							tree[index1] = tree[index2];
							tree[index2] = n;
						}
					break;
				}
			}
		}
	}
	stack<int> s1, s2, bigstack;
	temp = index1;
	s1.push(temp);
	while (tree[temp].Parent != none && (tree[tree[temp].Parent].id == V || tree[tree[temp].Parent].id == H)) {

		s1.push(tree[temp].Parent);
		temp = tree[temp].Parent;



	}

	temp = index2;
	while (tree[temp].Parent != none && (tree[tree[temp].Parent].id == V || tree[tree[temp].Parent].id == H)) {
		s2.push(tree[temp].Parent);
		temp = tree[temp].Parent;
	}

	while (!s1.empty() && !s2.empty()) {
		if (s1.top() == s2.top()) {
			bigstack.push(s1.top());
			s1.pop();
			s2.pop();
		}
		else {
			break;
		}
	}
	while (!s1.empty()) {
		bigstack.push(s1.top());
		s1.pop();
	}
	while (!s2.empty()) {
		bigstack.push(s2.top());
		s2.pop();
	}
	compute_area(bigstack, tree, NPE.size() - 1);

	return NPE;
}


double Compute_iniT(float p, vector<int> NPE, node* tree) {
	node* temp_tree;
	temp_tree = new node[2 * NumBlocks - 1]();
	vector<int> temp_NPE = NPE;
	cpyTree(temp_tree, tree, NPE.size());

	int cost = compute_cost();
	long long move, sum = 0, temp_cost;
	int cnt = initTcount;
	while (cnt-- > 0) {
		move = rand() % 3;
		if (move == 0) temp_NPE = M1(temp_NPE, temp_tree);
		else if (move == 1) temp_NPE = M2(temp_NPE, temp_tree);
		else if (move == 2) temp_NPE = M3(temp_NPE, temp_tree);
	
		temp_cost = compute_cost();
		if (temp_cost > cost) sum += (temp_cost - cost);
		cost = temp_cost;
	}
	return (double)(1.0*sum / initTcount) / (-1.0 * log(p));

}

vector<int> SA_floorplanning(vector<int> NPE, node* tree, float p, float e) {
	chrono::high_resolution_clock::time_point begin_time = time_record();
	chrono::high_resolution_clock::time_point end_time;
	
	vector<int> E = NPE;
	vector<int> NE;
	vector<int> bestE = E;

	bestTree = new node[2 * NumBlocks - 1]();
	dummyTree = new node[2 * NumBlocks - 1]();

	int done_count=0;
	cpyTree(bestTree, tree, NPE.size());
	cpyTree(dummyTree, tree, NPE.size());

	stack<int> nullstack;
	compute_area(nullstack, bestTree, E.size() - 1);
	double T = Compute_iniT(p, NPE, tree);
	cout << "initial T: " << T << endl;
	

	int  MT = 0, uphill = 0, reject = 0;
	bool successflag = false;
	int N = kmove * NumBlocks;

	cpyTree(tree, bestTree, NPE.size());
	compute_area(nullstack, bestTree, E.size() - 1);
	int cost = compute_cost();
	int dif_cost = 0;
	int newCost;
	int best_cost = cost;

	while (true) {
		MT = uphill = reject = 0;
		while (true) {
			
			//int move = 0;
/*		if (T > Tlimit0) move = rand() % 3;
		else if (T > Tlimit1) move = 0;
		else if (T > Tlimit2) move = rand() % 3;
		else move = 0;*/
			int move = T < Tlimit0 ? 0: rand() % 3;



			if (move == 0) NE = M1(E, tree);
			else if (move == 1) NE = M2(E, tree);
			else if (move == 2) NE = M3(E, tree);
			MT++;

			newCost = compute_cost();
			dif_cost = newCost - cost;
			cout << "T: " << T << " w: " << total_W << " h: " << total_H << " HPWL: " << HPWL << " cost: " << newCost << "\r";
		
			if (dif_cost <= 0 || (double)rand() / (RAND_MAX + 1.0) < exp(-1.0*dif_cost / T)) {
			
				if (dif_cost > 0) uphill++;
				E = NE;
				cpyTree(dummyTree, tree, NPE.size());
				cost = newCost;
				if (total_H < FixOutline && total_W < FixOutline && !successflag) {  successflag = true; }
				if (cost < best_cost) {
					if (!successflag || (total_H < FixOutline && total_W < FixOutline)) {
						bestE = E;
						best_cost = cost;
						if (successflag) {
							done_count++;
							cout << "\rsuccess2!!                                             " << HPWL << endl; }
						cpyTree(bestTree, tree, NPE.size());
					}
				}
			}
			else {
				reject++;
				cpyTree(tree, dummyTree, NPE.size());
			}
			if (uphill > N || MT > 2 * N) break;
		}
		T = coolrate * T;
		cout << "T: " << T << " W: " << total_W << " H: " << total_H << " HPWL: " << HPWL << endl;
		//if ((1.0*reject / MT) > 0.99) {
		end_time = time_record();
		if (  T< epsilon ||done_count>=1 || time_output(begin_time, end_time) > 1000) break;
			
		
		if (T < 1000) N = kmove *NumBlocks*3;
		else if (T < 10) N = kmove *NumBlocks*4;
		else if (T < 0.1) N = kmove *NumBlocks*5;
	}

	return bestE;
}


void output_function(char* arg) {
	ofstream output;
	output.open(arg);
	output << "Wirelength " << HPWL << endl;
	output << "Blocks" << endl;;
	for (int i = 0; i < NumBlocks; i++) {
		output << "sb" << i << " " << dicBlock[i]->x << " " << dicBlock[i]->y << " " << dicBlock[i]->rotate << endl;
	}
	output.close();
	return;
}


int main(int argc, char *argv[]) {
	srand(seed);
	int max_theads = omp_get_max_threads();
	omp_set_num_threads(max_theads);

	chrono::high_resolution_clock::time_point I_start = time_record();
	FILE* input_block = fopen(argv[1], "r");
	FILE* input_net = fopen(argv[2], "r");
	FILE* input_pl = fopen(argv[3], "r");
	srand(time(NULL));
	float ds_ratio = atof(argv[5]);
	parse_block(input_block);
	parse_net(input_net);
	parse_terminal(input_pl);
	chrono::high_resolution_clock::time_point I_end = time_record();
	double I_time = time_output(I_start, I_end);

	chrono::high_resolution_clock::time_point Ini_start = time_record();
	vector<int> NPE = ini_NPE();
	ini_Tree(NPE);
	stack<int> nullstack;
	compute_area(nullstack, iniTree, NPE.size() - 1);
	FixOutline = (int)sqrt(total_area*(1 + ds_ratio));
	chrono::high_resolution_clock::time_point Ini_end = time_record();
	double Ini_time = time_output(Ini_start, Ini_end);

	chrono::high_resolution_clock::time_point Compu_start = time_record();
	NPE = SA_floorplanning(NPE, iniTree, probability, MinTemperature);
	chrono::high_resolution_clock::time_point Compu_end = time_record();
	double Compu_time = time_output(Compu_start, Compu_end);

	compute_area(nullstack, bestTree, NPE.size() - 1);
	compute_cost();
	cout << "\n[BestE] W:" << total_W << " H:" << total_H << " wirelenth: " << HPWL << endl;
	cout << "Boundary= " << FixOutline << endl;
	cout <<"total_area= "<< total_area << endl;
	cout << "my deadspace ratio = " << (float)total_area /( (float)FixOutline * (float)FixOutline) << endl;
	if (total_H <= FixOutline && total_W <= FixOutline) cout << "[[success fit!!]]" << endl;

	Compute_HPWL();
	chrono::high_resolution_clock::time_point O_start = time_record();
	output_function(argv[4]);
	chrono::high_resolution_clock::time_point O_end = time_record();
	double O_time = time_output(O_start, O_end);

	//ofstream output2;
	//string str;
	//str.assign(argv[6]);
	//str.erase(0, 12);
	//output2.open(argv[6]);
	//output2 << str << " deadspace ratio: "<<argv[5] << endl;
	//output2 << "HPWL: " << HPWL << endl;
	//output2 << "FixOutline: " << FixOutline<<" W: "<<total_W<<" H:"<< total_H << endl;
	//if (total_H <= FixOutline && total_W <= FixOutline) output2 << "Great fit!!" << endl;
	//if (total_H <= FixOutline && total_W > FixOutline) output2 << "W failed" << endl;
	//if (total_H > FixOutline && total_W <= FixOutline) output2 << "H failed" << endl;
	//if (total_H > FixOutline && total_W > FixOutline) output2 << "W and H both failed" << endl;
	//output2 << "\n[Parameters]" << endl;
	//output2 << "weight_outline_H: " << weight_outline_H << "\n";
	//output2 << "weight_outline_W: "<< weight_outline_W << "\n";
	//output2 << "deadweight: "<< HPWLweight <<"(no * deadratio)\n";
	//output2 << "Areaweight: "<<areaweight<<"\n";
	//output2 << "Seed "<<seed<<"\n";

	//output2.close();
	cout<<"I/O time: "<<O_time + I_time <<endl;
	cout << "Ini_time: " << Ini_time << endl;
	cout<<"Compute time: "<< Compu_time <<endl;
	cout << "Total time: " << O_time + I_time + Compu_time+ Ini_time << endl;

	return 0;
}

