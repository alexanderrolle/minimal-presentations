#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <iterator>
#include <algorithm>
#include <array>
#include <chrono>
#include <cassert>
#include <unordered_map>
#include <limits>
#include <iomanip>

using namespace std;

using chrono::high_resolution_clock;
using chrono::duration_cast;
typedef chrono::milliseconds TimeUnit;

//------------ENTRIES OF THE BOUNDARY------------------------------//

struct entry{
  int simplex = -1;
  int coefficient = -1;

  entry() {};
  entry(int s, int c) : simplex(s), coefficient(c) {};
};

struct compareEntry{
  bool operator ()(const entry &a, const entry &b) const
  {
    return a.simplex < b.simplex;
  }
};

//------------ArrayHasher------------------------------//

struct ArrayHasher {
    std::size_t operator()(const std::array<int, 2>& a) const {
        std::size_t h = 0;

        for (auto e : a) {
            h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

//------------COLUMNS------------------------------//

class column {
public:
  int index;
  int dim;
  array<float, 2> value;
  priority_queue<entry, vector<entry>, compareEntry> boundary;
  int local;   //  -1 local negative, 0 global, 1 local positive
  int pair;
};

struct compare{
  bool operator ()(const column &a, const column &b) const
  {
    if(a.value[0] == b.value[0]) {
      if(a.value[1] == b.value[1])
        {return a.dim < b.dim;}
      return a.value[1] < b.value[1];
    }
    return a.value[0] < b.value[0];
  }
};

bool hasPivot(vector<column> &v, int i)
{
  if(v[i].boundary.empty()==false)
    {
      if(v[v[i].boundary.top().simplex].value==v[i].value)
      {
        return true;
      }
    }
  return false;
}

void compress(priority_queue <entry, vector<entry>, compareEntry> &b)
{
  priority_queue <entry, vector<entry>, compareEntry> tempb;
  while(b.empty()==false)
  {
    entry e, f;
    e=b.top();
    b.pop();
    if(b.empty()==false)
      {
        f=b.top();
        if(e.simplex==f.simplex)
          {
            e.coefficient=(e.coefficient+f.coefficient)%2;
            b.pop();
          }
      }
    if(e.coefficient!=0)
      {
        tempb.push(e);
      }
  }
  b=tempb;
}

void addColumns(vector<column> &v, int i, int j, int l) // i -> i + l x j
{
  priority_queue <entry, vector<entry>, compareEntry> bj;
  bj=v[j].boundary;
  while(bj.empty()==false)
  {
    entry e;
    e.simplex=bj.top().simplex;
    e.coefficient=l*bj.top().coefficient;
    v[i].boundary.push(e);
    bj.pop();
  }
  compress(v[i].boundary);
}

int checkBefore(vector<column> &v, vector<int> vv, int i)
{
  int c=vv[i];
  for(int j=0; j<i; j++ )
    {
      int cc=vv[j];
      if(v[c].boundary.empty()==false && v[cc].boundary.empty()==false && v[c].boundary.top().simplex==v[cc].boundary.top().simplex)
      {
        return cc;
      }
    }
  return(-1);
}

//------------PRINT FUNCTIONS------------------------------//

template <class myType>

void printVector(vector<myType> v)
{
  for(auto it = v.begin() ; it!=v.end() ; it++)
  {
    cout << *it << ' ';
  }
  cout << '\n';
}

void printEntry(entry e)
{
    cout << "("<< e.simplex << "," << e.coefficient  << ")"<< '\t';
}

void printQueue(priority_queue <entry, vector<entry>, compareEntry> q)
{
    priority_queue <entry, vector<entry>, compareEntry> g = q;
    while (!g.empty())
    {
        // printEntry(g.top());
        cout << g.top().simplex << ' ';
        g.pop();
    }
    cout << '\n';
}

// void printColumns(vector<column> v)
// {
//   for(auto it = v.begin() ; it!=v.end() ; it++)
//   {
//     cout << it->dim << ' ';
//     printQueue(it->boundary);
//     cout << ':' << ' '<< it->value[0] << ' ' << it->value[1] << '\n';
//   }
// }

void print_firep(vector<column>& v) {

  int dim_0=0,dim_1=0,dim_2=0;

  int index_0=0,index_1=0,index_2=0;

  std::unordered_map<int,int> re_index_0;
  std::unordered_map<int,int> re_index_1;
  std::unordered_map<int,int> re_index_2;

  

  for(auto it = v.begin() ; it!=v.end() ; it++)
  {
    if(it->dim==0) {
      re_index_0[it->index]=index_0++;
      dim_0++;
    }
    if(it->dim==1) {
      re_index_1[it->index]=index_1++;
      dim_1++;
    }
    if(it->dim==2) {
      re_index_2[it->index]=index_2++;
      dim_2++;
    }
  }
  cout << "firep\nx-coordinate\ny-coordinate\n" << dim_2 << " "<< dim_1 << " " << dim_0 << std::endl;

   for(auto it = v.begin() ; it!=v.end() ; it++)
   {
     if(it->dim==2) {
       cout << fixed << setprecision(12) << it->value[0] << " " <<fixed << setprecision(12) << it->value[1] << " ; ";
       priority_queue <entry, vector<entry>, compareEntry> g = it->boundary;
       vector<int> new_entries;
       while (!g.empty()) {
	 new_entries.push_back(re_index_1[g.top().simplex]);
	 g.pop();
       }
       sort(new_entries.begin(),new_entries.end());
       for(auto vit = new_entries.begin();vit!=new_entries.end();vit++) {
	 cout << *vit << " ";
       }
       cout << endl;
    }
   }
   for(auto it = v.begin() ; it!=v.end() ; it++)
   {
     if(it->dim==1) {
       cout << it->value[0] << " " << it->value[1] << " ; ";
       priority_queue <entry, vector<entry>, compareEntry> g = it->boundary;
       vector<int> new_entries;
       while (!g.empty()) {
	 new_entries.push_back(re_index_0[g.top().simplex]);
	 g.pop();
       }
       sort(new_entries.begin(),new_entries.end());
       for(auto vit = new_entries.begin();vit!=new_entries.end();vit++) {
	 cout << *vit << " ";
       }
       cout << endl;
     }
   }

}

void printSurvivedColumns(vector<column>& v)
{
  for(auto it = v.begin() ; it!=v.end() ; it++)
  {
    cout << "Index: " <<  it->index << ", ";
    cout << "Value: " << "(" <<  it->value[0] << ", " << it->value[1] << ")" << ", ";
    cout << "Dimension: " <<  it->dim << ", ";
    cout << "Boundary: ";
    printQueue(it->boundary);
  }
}

//------------------------------------------//

int main (int argc, char* argv[])
{
  // cout << "Where not specified time is misured in milliseconds\n";
  high_resolution_clock::time_point tpStart, tpTmp;
  tpStart = tpTmp = high_resolution_clock::now();
  //MemoryUsage mem;
  // time_t tTmp, tStart;
  // tTmp = tStart = time(0);

  //--------------READING THE INPUT----------------------------//

  string line;
  ifstream infile(argv[1]);
  vector<column> InputSimplices;
  vector<int> Input;
  int d=0;
  int NV;
  int NF;
  if (infile.is_open())
  {

    // getline (infile,line);
    //
    getline (infile,line);
    getline (infile,line);
    stringstream ss(line);
    ss >> NV;
    ss >> NF;
    for(int i=0; i<NV; i++)
    {
      column c;
      c.index=i;
      c.dim=0;
      c.local=0;
      getline (infile,line);
      istringstream iss(line);
      iss >> c.value[0];
      iss >> c.value[1];
      InputSimplices.push_back(c);
    }

    for(int i=0; i<NF; i++)
    {
      getline (infile,line);
      istringstream iss(line);
      int j;
      iss >> j;
      for(int k=0; k<j; k++)
      {
        int v;
        iss >> v;
        Input.push_back(v);
      }

    }

  }

  else cout << "Unable to open file";

  infile.close();

// printSurvivedColumns(InputSimplices);


//--------------CREATING THE BOUNDARY MATRIX (AN UNORDERED VERSION OF THAT)----------------------------//

  auto tpStartMeshConstr=high_resolution_clock::now();
  //cout << "\n";
  //mem.getValue_in_MB(true);

  unordered_map< array<int,2> , int, ArrayHasher >  umap;
  int index=NV;

  for(int i=0; i<NF; i++)
  {
    column c;
    c.dim=2;
    d=max(d,c.dim);
    c.local=0;
    c.value[0]=numeric_limits<int>::min();
    c.value[1]=numeric_limits<int>::min();
    // c.value[0]=-1000;
    // c.value[1]=-1000;
    for(int m=0; m<3; m++)
    {
      array<int,2> e;
      int edgeIndex;
      array<float, 2> maxEdgeValue;
      e[0]=Input[3*i+(m%3)];
      e[1]=Input[3*i+((m+1)%3)];
      sort(e.begin(), e.end());
      if (umap.find(e) == umap.end())
      {
        // cout << "new" << "\n";
        umap[e]=index;
        c.boundary.emplace(index, 1);
        column dd;
        dd.dim=1;
        dd.index=index;
        dd.local=0;
        // cout << index << "\n";
        dd.boundary.emplace(e[0], 1);
        dd.boundary.emplace(e[1], 1);
        dd.value[0]=max(InputSimplices[e[0]].value[0], InputSimplices[e[1]].value[0]);
        dd.value[1]=max(InputSimplices[e[0]].value[1], InputSimplices[e[1]].value[1]);
        InputSimplices.push_back(dd);
        // cout << c.value[0] << " " << dd.value[0] << "\n";
        c.value[0]=max(c.value[0],dd.value[0]);
        c.value[1]=max(c.value[1],dd.value[1]);
        // cout << c.value[0] << "\n";
        index++;
      }
      else
      {
        // cout << e[0] << " " << e[1] << "\n";
        // cout << umap.find(e)->second << "\n";
        c.boundary.emplace(umap.find(e)->second,1);
        c.value[0]=max(c.value[0],InputSimplices[umap.find(e)->second].value[0]);
        c.value[1]=max(c.value[1],InputSimplices[umap.find(e)->second].value[1]);
        // cout << c.value[0] << "qui" <<"\n";
      }
    }

    c.index=index;
    index++;

    // iss >> c.value[0];
    // cout << c.value[0] << "\n";
    // iss >> c.value[1];
    // cout << c.value[1] << "\n";
    InputSimplices.push_back(c);
  }

  // printSurvivedColumns(InputSimplices);

  //------REORDER THE COLUMNS ACCORDING WITH A TOTAL ORDER (respecting the filtration values)---------------//

  auto tpStartReorder=high_resolution_clock::now();
  //mem.getValue_in_MB(true);

  vector<column> &OrderedSimplices=InputSimplices;
  sort(OrderedSimplices.begin(), OrderedSimplices.end(), compare());

  vector<int> mappa(OrderedSimplices.size());
  for(size_t i=0; i<OrderedSimplices.size(); i++)
  {
    mappa[OrderedSimplices[i].index]=i;
  }

  for(size_t i=0; i<OrderedSimplices.size(); i++)
  {
    OrderedSimplices[i].index=i;
    priority_queue <entry, vector<entry>, compareEntry> tempQueue;
    while (!OrderedSimplices[i].boundary.empty())
    {
        entry e;
        e.simplex=mappa[OrderedSimplices[i].boundary.top().simplex];
        e.coefficient=1;
        OrderedSimplices[i].boundary.pop();
        tempQueue.push(e);
    }
    OrderedSimplices[i].boundary=tempQueue;
  }

  vector< vector<int> > SimplicesOfDim;
  for(int k=0; k<d+1; k++)
  {
    vector<int> myvector;
    for(auto it = OrderedSimplices.begin() ; it!=OrderedSimplices.end() ; it++)
    {
      if(it->dim==k)
      {
        myvector.push_back( it->index );
      }
    }
    SimplicesOfDim.push_back(myvector);
  }

  // cout << "\n";
  // cout << "--# Columns in Input: "<< OrderedSimplices.size() <<"-----------------\n";
  //printSurvivedColumns(OrderedSimplices);
  print_firep(OrderedSimplices);

  return 0;
}
