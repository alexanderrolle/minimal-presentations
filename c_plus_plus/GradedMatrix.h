#include<phat/boundary_matrix.h>
#include<algorithm>

#include<unordered_map>
#include <boost/functional/hash.hpp>

#include <string>
#include <cstdio>
#include <cerrno>

// Note: If USE_DOUBLE is switched off, a boost number type is used which
// in principle results in more robust code. However, we experienced
// problems with that code (when converting the float type into a
// rational type), so we represent the grade coordinates with doubles,
// but store the exact input string as well to be able to convert it
// to a raional later (using rationals throughout yields in a performance
// penalty)
#define USE_DOUBLE 1

#if !USE_DOUBLE
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#endif

#include<gmpxx.h>

namespace phat {

  typedef std::pair<index,index> index_pair;

  typedef std::priority_queue<index,std::vector<index>,std::greater<index>> PQ;

  // This map is artificial, but needed in order to produce the rational representation as rivet.
  std::unordered_map<double,std::string> double_map;

  mpq_class rat_representation(double d) {
    auto it = double_map.find(d);
    if(it==double_map.end()) {
      std::cout << "DOUBLE " << d << " not found" << std::endl;
      std::exit(1);
    }
    const char* str = it->second.c_str();
    //std::cout << "Converting " << str << std::endl;
    int sign=+1;
    int start=0;
    if(str[0]=='-') {
      start++;
      sign=-1;
    }
    int n = strlen(str);
    mpz_class num(0), denom(1);
    bool dp_detected=false;
    for(int i = start;i<n;i++) {
      if(str[i]=='.') {
	dp_detected=true;
      } else {
	int digit = (int)str[i]-48;
	num=10*num+digit;
	if(dp_detected) {
	  denom*=10;
	}
      }
    }
    mpq_class result(num*sign,denom);
    result.canonicalize();
    return result;
  }




#if USE_DOUBLE

  typedef double Coordinate;

  Coordinate from_string(char* str) {
    double result = atof(str);
    double_map[result]=std::string(str);
    //std::cout << "Mapping " << result << " to " << str << std::endl;
    return result;
    
  }
#else
  
#if 0

  // This seems to work, but it is infintiely slow

  // For precise comparison with Rivet output, this needs to be used
  typedef boost::multiprecision::cpp_rational Coordinate;
  
  Coordinate from_string(char* str) {
    //std::cout << "Converting " << str << " to mpq_class" << std::endl;
    boost::multiprecision::cpp_dec_float_100 dec_f(str);
    //std::cout << "cpp_dec: " << dec_f << std::endl;
    Coordinate result = dec_f.convert_to<boost::multiprecision::cpp_rational>();
    //std::cout << "Converted " << str << " to " << result << std::endl;
    return result;
  }

#else
  // Using slow boost types
  typedef boost::multiprecision::cpp_dec_float_100 Coordinate;
  
  Coordinate from_string(char* str) {
    //std::cout << "Converting " << str << " to mpq_class" << std::endl;
    Coordinate result(str);
    return result;
  }

#endif  


#endif // of USE_DOUBLE

  struct Grade {
    Coordinate first_val;
    index first_index;
    Coordinate second_val;
    index second_index;
    Grade() {}
    Grade(Coordinate& x, Coordinate& y) : first_val(x), second_val(y) {}
    Grade(index ind_x, index ind_y, Coordinate& val_x, Coordinate& val_y) : first_val(val_x), first_index(ind_x), second_val(val_y), second_index(ind_y) {}
    Grade(const Grade& other) : first_val(other.first_val), first_index(other.first_index),
				second_val(other.second_val), second_index(other.second_index) {}
    bool operator== (const Grade& other) {
      return this->first_val==other.first_val && this->second_val==other.second_val;
    }
  };

  struct pre_column {
    index idx;
    Grade grade;
    std::vector<index> boundary;
    pre_column() {}
    pre_column(index idx, Grade& grade, std::vector<index>& boundary)
      : idx(idx), grade(grade), boundary(boundary) {
      std::sort(boundary.begin(),boundary.end());
    }
  };

  template<typename Column>
  struct Sort_pre_column {
    bool operator() (const Column& c1, const Column& c2) {
      if(c1.grade.second_val < c2.grade.second_val) {
	return true;
      }
      if(c1.grade.second_val > c2.grade.second_val) {
	return false;
      }
      return c1.grade.first_val < c2.grade.first_val;
    }
  };

  struct Sort_grades {
    bool operator() (const index_pair& c1, const index_pair& c2) {
      if(c1.first > c2.first) {
	return true;
      }
      if(c1.first < c2.first) {
	return false;
      }
      return c1.second > c2.second;
    }
  };

#if SPARSE_GRID_TRAVERSAL
  class Grid_scheduler {

  public:

    std::priority_queue<index_pair,std::vector<index_pair>,Sort_grades> grades;

    std::map<index_pair,index_pair> index_range;

    index_pair curr_grade;

    Grid_scheduler() {}

    // It is assumed that columns with the same grade appear in M consecutively
    template<typename GradedMatrix>
      Grid_scheduler(GradedMatrix& M) {
      
      index_pair last_pair=std::make_pair(-1,-1);
      index curr_start=-1;
      for(int i=0;i<M.get_num_cols();i++) {
	index curr_x=M.grades[i].first_index;
	index curr_y=M.grades[i].second_index;
	if(curr_x!=last_pair.first || curr_y !=last_pair.second) {
	  // New grade
	  if(curr_start!=-1) {
	    index_range[last_pair]=std::make_pair(curr_start,i);
	  }
	  curr_start=i;
	  last_pair = std::make_pair(curr_x,curr_y);
	  grades.push(last_pair);
	}
      }
      if(curr_start!=-1) {
	index_range[last_pair]=std::make_pair(curr_start,M.get_num_cols());
      }
      curr_grade=std::make_pair(-1,-1);
	  
    }

    bool at_end() {
      return grades.empty();
    }

    index_pair next_grade() {
      index_pair result = grades.top();
      grades.pop();
      while(!grades.empty() && grades.top()==result) {
	grades.pop();
      }
      curr_grade=result;
      return result;
    }

    index_pair index_range_at(index x, index y) {
      auto find_grade = index_range.find(std::make_pair(x,y));
      if(find_grade==index_range.end()) {
	return std::make_pair(0,0);
      }
      return find_grade->second;
    }

    void notify(index x,index y) {
      //std::cout << "Got notified about " << x << " " << y << std::endl;
      if(curr_grade.first!=x || curr_grade.second!=y) {
	grades.push(std::make_pair(x,y));
      }
    }
    
    // extended_index_range_at not needed, because SMART_REDUCTION is required

    

  };
#else

  class Grid_scheduler {
    
    index_pair _curr_pair;
    bool _at_end;

    //std::unordered_map<index_pair,index,boost::hash<index_pair> > start_index_of_pair;
    //std::map<index_pair,index> start_index_of_pair;
    std::vector<std::vector<index>> start_index_of_pair;

    index num_grades_x;
    index num_grades_y;

    index num_cols;

    

  public:

    Grid_scheduler() {}
    
    template<typename GradedMatrix>
      Grid_scheduler(GradedMatrix& M) {
      
      num_cols = M.get_num_cols();
      
      num_grades_x=M.num_grades_x;
      num_grades_y=M.num_grades_y;

      index run_x=0;
      index run_y=0;
      
      index curr_x,curr_y;
      
      start_index_of_pair.resize(num_grades_x);
      for(int i=0;i<M.num_grades_x;i++) {
	start_index_of_pair[i].resize(num_grades_y);
      }
      
      for(int i=0;i<=num_cols;i++) {
	
	if(i<num_cols) {
	  curr_x=M.grades[i].first_index;
	  curr_y=M.grades[i].second_index;
	} else {
	  curr_x=num_grades_x-1;
	  curr_y=num_grades_y-1;
	}
	
	//assert(run_x<=curr_x);
	//assert(run_y<=curr_y);
	
	while(run_y<curr_y || (run_y==curr_y && run_x <= curr_x)) {
	  //std::cout << "M1: Set index of " << run_x << " " << run_y << " to " << i << std::endl;
	  //std::cout << "Curr: " << curr_x << " " << curr_y << std::endl;
	  start_index_of_pair[run_x][run_y]=i;
	  run_x++;
	  if(run_x==num_grades_x) {
	    run_y++;
	    run_x=0;
	  }
	}
      }
      _curr_pair=std::make_pair(0,0);
      _at_end=(num_grades_x==0 && num_grades_y==0);
    }
    
    bool at_end() {
      return _at_end;
    }

    index_pair next_grade() {
      index_pair result=_curr_pair;

      if(_curr_pair.first==num_grades_x-1 && _curr_pair.second==num_grades_y-1) {
	_at_end=true;
      } else {
	if(_curr_pair.second==num_grades_y-1) {
	  _curr_pair.second=0;
	  _curr_pair.first++;
	} else {
	  _curr_pair.second++;
	}
      }	

      return result;
    }

    index_pair index_range_at(index x, index y) {
      index end_xy;
      if(x<num_grades_x-1) {
	end_xy = start_index_of_pair[x+1][y];
      } else {
	end_xy = start_index_of_pair[0][y+1];
      }
      if(x==num_grades_x-1 && y==num_grades_y-1) {
	end_xy=num_cols;
      }
      index start_xy = start_index_of_pair[x][y];
      return std::make_pair(start_xy,end_xy);
    }
    
    // Needed for the original version (without SMART_REDUCTION)
    void extended_index_range_at(index x,index y,
				 index& start_xy,
				 index& end_xy,
				 index& start_xy_of_grade) {

      index_pair index_range=this->index_range_at(x,y);
      start_xy_of_grade=index_range.first;
      end_xy=index_range.second;
      start_xy = start_index_of_pair[0][y];
    }

    
    void notify(index x,index y) {
      // Do nothing
    }

  };
#endif // of SPARSE_GRID_TRAVERSAL

  template<typename Representation=vector_heap>
    class GradedMatrix  : public boundary_matrix<Representation> {
    
    public:

    std::vector<Coordinate> x_vals,y_vals;

    index num_grades_x;
    index num_grades_y;

    index num_rows;
    
    Grid_scheduler grid_scheduler;
    
    std::vector<Grade> grades;

    std::vector<Grade> row_grades;
    
    std::vector<index> pivots;

#if SMART_REDUCTION
    std::vector<PQ> pq_row;
#endif

#if CLEARING
    std::unordered_map<index,index> clearing_info;
#endif

    boundary_matrix<Representation> slave;

    /*
    void clear() {
      x_vals.clear();
      y_vals.clear();
      start_index_of_pair.clear();
      grades.clear();
      row_grades.clear();
      pivots.clear();
      for(int i=0;i<this->get_num_cols();i++) {
	this->clear(i);
      }
    }
    */

#if !LAZY_MINIMIZATION

    // Needed only for vector_vector column type
    // This is not following the phat interface, but uses the internal structure of columns (pure hack to avoid excessive copying of columns)
    bool contains(index i, index m) {
      return std::binary_search(this->rep.matrix[i].indices.begin(),this->rep.matrix[i].indices.end(),m);
    }

#endif

    void print(bool print_row_grades=false,bool print_indices=true) {
      std::cout << "Number of columns: " << this->get_num_cols() << std::endl;
      std::cout << "Number of rows: " << this->num_rows << std::endl;

      for(int i=0;i<this->get_num_cols();i++) {
	if(print_indices) {
	  std::cout << grades[i].first_index <<  " " << grades[i].second_index << " -- ";
	} else {
	  std::cout << grades[i].first_val <<  " " << grades[i].second_val << " -- ";
	}
	column col;
	this->get_col(i,col);
	for(std::size_t j=0;j<col.size();j++) {
	  std::cout << col[j] << " ";
	}
	std::cout << std::endl;

      }

      if(print_row_grades) {
	std::cout << "Row grades:" << std::endl;
	for(index j=0;j<this->num_rows;j++) {
	  std::cout << j << " : ";
	  if(print_indices) {
	    std::cout << row_grades[j].first_index <<  " " << row_grades[j].second_index << std::endl;
	  } else {
	    std::cout << row_grades[j].first_val <<  " " << row_grades[j].second_val << std::endl;
	  }
	}
      }
      
    }

    bool is_local(index i) {
      index p = this->get_max_index(i);
      //std::cout << "Info: " << i << " "<<  p  << std::endl;
      return p!=-1 && (this->grades[i].first_index == this->row_grades[p].first_index) &&  (this->grades[i].second_index == this->row_grades[p].second_index) ;
    }

    bool pivot_is_dominating(index i) {
      if(this->is_empty(i)) {
	return false;
      }
      if(this->is_local(i)) {
	return true;
      }
      index p = this->get_max_index(i);
      Grade& pgr = this->row_grades[p];
      //std::cout << "#### NEW COLUMN ####" << std::endl;
      //std::cout << "pivot grade is " << pgr.first_index << ", " << pgr.second_index << std::endl;
      std::vector<index> col;
      this->get_col(i,col);
      for(int j=0;j<col.size();j++) {
	Grade& curr=this->row_grades[col[j]];
	//std::cout << "next grade is " << curr.first_index << ", " << curr.second_index << std::endl;
	if(curr.first_index>pgr.first_index || curr.second_index>pgr.second_index) {
	  //std::cout << "NOT HERE" << std::endl;
	  return false;
	}
      }
      
      return true;
    }

    void reduce_column(index i, index_pair& curr_gr, bool use_slave=false, bool notify_pq=false) {
      

      /*
      std::cout << "Pivots vector: " << std::endl;
      for(int i=0;i<this->num_rows;i++) {
	std::cout << i << " -> " << this->pivots[i] << std::endl;
      }
      */

	

      index p = this->get_max_index(i);
      
      //std::cout << "p=" << p << std::endl;

      while(p!=-1 && pivots[p]!=-1 && pivots[p]<i) {
	index k = pivots[p];

	//std::cout << "Adding column " << k << std::endl;
	

	this->add_to(k,i);
	gl_no_column_additions++;

	if(use_slave) {
	  slave.add_to(k,i);
	  gl_no_column_additions++;
	}
	p=this->get_max_index(i);
      }
#if SMART_REDUCTION
      if(notify_pq && p!=-1 && pivots[p]>i) {

	index j = pivots[p];
	index gr_y_index = this->grades[j].second_index;
	//std::cout << "SCHEDULING COLUMN FOR LATER " << i << " " << j << " " << this->grades[j].first_index<< " " << gr_y_index << std::endl;
	this->pq_row[gr_y_index].push(j);

#if SPARSE_GRID_TRAVERSAL
	index gr_x_index = curr_gr.first;
	this->grid_scheduler.notify(gr_x_index,gr_y_index);
#endif
      }
#endif
      if(p!=-1 && (pivots[p]==-1 || pivots[p]>i)) {
	pivots[p]=i;
      }
    }

    

    void assign_pivots() {
      for(index i=0;i<this->num_rows;i++) {
	pivots.push_back(-1);
      }
    }



    void assign_slave_matrix() {
      index n = this->get_num_cols();
      slave.set_num_cols(n);
      for(int i=0;i<n;i++) {
	std::vector<index> col;
	col.push_back(i);
	slave.set_col(i,col);
      }
    }

  };


  template<typename GradedMatrix>
    void check_boundaries(GradedMatrix& M, std::string msg) {
    for(index i=0;i<M.get_num_cols();i++) {
      std::vector<index> col;
      M.get_col(i,col);
      for(index j=1;j<col.size();j++) {
	if(col[j-1]>=col[j]) {
	  std::cout << "Bad boundary " << col[j-1] << " " << col[j] << " at column " << i << "(" << msg << " matrix)" << std::endl;
	  std::exit(1);
	}
      }
    }
  }

  template<typename GradedMatrix>
    void check_grade_sanity(GradedMatrix& M) {
    
    for(index i=0;i<M.get_num_cols();i++) {
      Grade& col_gr = M.grades[i];
      std::vector<index> col;
      M.get_col(i,col);
      for(index j=0;j<col.size();j++) {
	Grade& row_gr = M.row_grades[col[j]];
	if(row_gr.first_index>col_gr.first_index || row_gr.second_index>col_gr.second_index) {
	  std::cout << "Bad grading at: " << col[j] << "( " << row_gr.first_index << ", " << row_gr.second_index << ") " << i << "( " << col_gr.first_index << ", " << col_gr.second_index << ")" << std::endl;
	  std::exit(1);
	}
	assert(row_gr.first_index<=col_gr.first_index);
	assert(row_gr.second_index<=col_gr.second_index);
      }
    }
    

  }

  
  template<typename GradedMatrix>
    void assign_index_range(GradedMatrix& M) {
   int n = M.get_num_cols();

    index run_x=0;
    index run_y=0;
    
    index curr_x,curr_y;

    M.start_index_of_pair.resize(M.num_grades_x);
    for(int i=0;i<M.num_grades_x;i++) {
      M.start_index_of_pair[i].resize(M.num_grades_y);
    }

    for(int i=0;i<=n;i++) {
      
      if(i<n) {
	curr_x=M.grades[i].first_index;
	curr_y=M.grades[i].second_index;
      } else {
	curr_x=M.num_grades_x-1;
	curr_y=M.num_grades_y-1;
      }
      
	//assert(run_x<=curr_x);
	//assert(run_y<=curr_y);
      
      while(run_y<curr_y || (run_y==curr_y && run_x <= curr_x)) {
	//std::cout << "M1: Set index of " << run_x << " " << run_y << " to " << i << std::endl;
	//std::cout << "Curr: " << curr_x << " " << curr_y << std::endl;
	M.start_index_of_pair[run_x][run_y]=i;
	run_x++;
	if(run_x==M.num_grades_x) {
	  run_y++;
	  run_x=0;
	}
      }
    }
 
    
  }

  // Assumes that columns are already added to matrix in lex order
  template<typename GradedMatrix>
    void assign_grade_indices(GradedMatrix& M1, GradedMatrix& M2) {
    
    int n1=M1.get_num_cols();
    int n2=M2.get_num_cols();

    //std::cout << "Assgin grade indices with " << n1 << ", " << n2 << " columns" << std::endl;
    if(n1==0 && n2==0) {
      return;
    }
    
#if USE_DOUBLE
    std::unordered_map<Coordinate,index> val_to_index_x, val_to_index_y;
#else
    std::map<Coordinate,index> val_to_index_x, val_to_index_y;
#endif
    
    std::vector<Coordinate> x_vals,y_vals;
    
    for(int i=0;i<n1;i++) {
      x_vals.push_back(M1.grades[i].first_val);
      y_vals.push_back(M1.grades[i].second_val);
    }
    for(int i=0;i<n2;i++) {
      x_vals.push_back(M2.grades[i].first_val);
      y_vals.push_back(M2.grades[i].second_val);
    }
    
    std::sort(x_vals.begin(),x_vals.end());
    auto last_x = std::unique(x_vals.begin(),x_vals.end());
    x_vals.erase(last_x,x_vals.end());
    std::sort(y_vals.begin(),y_vals.end());
    auto last_y = std::unique(y_vals.begin(),y_vals.end());
    y_vals.erase(last_y,y_vals.end());
    
    std::cout << "Found " << x_vals.size() << " different x-values and " << y_vals.size() << " different y-values" << std::endl;
    M1.x_vals=x_vals;
    M1.y_vals=y_vals;
    M1.num_grades_x = x_vals.size();
    M1.num_grades_y = y_vals.size();
    M2.x_vals=x_vals;
    M2.y_vals=y_vals;
    M2.num_grades_x = x_vals.size();
    M2.num_grades_y = y_vals.size();
    
    for(index i=0;i<x_vals.size();i++) {
      val_to_index_x[x_vals[i]]=i;
    }
    for(index i=0;i<y_vals.size();i++) {
      val_to_index_y[y_vals[i]]=i;
    }
    
    for(int i=0;i<n1;i++) {
      M1.grades[i].first_index=val_to_index_x[M1.grades[i].first_val];
      M1.grades[i].second_index=val_to_index_y[M1.grades[i].second_val];
      }
    for(int i=0;i<n2;i++) {
      M2.grades[i].first_index=val_to_index_x[M2.grades[i].first_val];
	M2.grades[i].second_index=val_to_index_y[M2.grades[i].second_val];
    }
  
    
    std::cout << "Num grades x " << M1.num_grades_x << std::endl;
    std::cout << "Num grades y " << M1.num_grades_y << std::endl;
    
    //std::cout << "n1=" << n1 << std::endl;
    //std::cout << "n2=" << n2 << std::endl;
    std::cout << "N is " << n1+n2 << std::endl;
  }

  template<typename GradedMatrix>
    void assign_grade_indices(GradedMatrix& M1) {

    GradedMatrix dummy;
    assign_grade_indices(M1,dummy);
  }

  struct File_reader {
    
    int pos;
    
    std::string data;
    
    File_reader(const char *filename) : pos(0) {
      // Copied from http://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
      std::FILE *fp = std::fopen(filename, "rb");
      if (fp) {
	std::fseek(fp, 0, SEEK_END);
	data.resize(std::ftell(fp));
	std::rewind(fp);
	// Return value is stored t supress compiler warning
	std::size_t no_read=std::fread(&data[0], 1, data.size(), fp);
	std::fclose(fp);
      } else {
	throw(errno);
      }
    }
    
    std::string next_line() {
      int old_pos=pos;
      while(data[pos]!='\n') {
	pos++;
      }
      if(pos!=data.length()) {
	pos++;
      }
      return data.substr(old_pos,pos-old_pos-1);
    }
  };
 
  // The sign is just a hack for the perturbation scenario and normally has no effect
  void load_prematrix_contents(File_reader& reader,
			       std::vector<pre_column>& pre_matrix,
			       int n,
			       int sign) {


#if 1

#if PERTURB
    std::uniform_real_distribution<double> unif(0,0.01);
    std::default_random_engine re;
#endif
    pre_matrix.resize(n);
    std::vector<std::string> lines;
    lines.resize(n);
    //std::string line;
    for(int i=0;i<n;i++) {

      lines[i]=reader.next_line();
    }
#if STRINGTOK_R_AVAILABLE && PARALLEL_FOR_LOOPS
#pragma omp parallel for schedule(guided,1)
#endif
    for(int i=0;i<n;i++) {
      pre_column& curr = pre_matrix[i];
      curr.idx=i;
      char* line = (char*)lines[i].c_str();
#if STRINGTOK_R_AVAILABLE
      char* saveptr;
      char* token = strtok_r(line," ",&saveptr);
#else
      char* token = strtok(line," ");
#endif
      
      curr.grade.first_val = from_string(token);
      //std::cout << "value: " <<  curr.grade.first_val << std::endl;
#if PERTURB
      curr.grade.first_val+=sign*unif(re);
#endif
#if STRINGTOK_R_AVAILABLE
      token = strtok_r(NULL," ",&saveptr);
#else
      token = strtok(NULL," ");
#endif
      curr.grade.second_val = from_string(token);
#if PERTURB
      curr.grade.second_val+=sign*unif(re);
#endif
#if SWAP_GRADE
      std::swap(curr.grade.first_val,curr.grade.second_val);
#endif
#if STRINGTOK_R_AVAILABLE
      token = strtok_r(NULL," ",&saveptr);
#else
      token = strtok(NULL," ");
#endif
      if(strcmp(token,";")) {
	std::cerr << "Semicolon missing" << std::endl;
      }
      std::vector<index> indices;
#if STRINGTOK_R_AVAILABLE
      while(token = strtok_r(NULL," ",&saveptr) ) {
#else
      while(token = strtok(NULL," ") ) {
#endif
	curr.boundary.push_back(atoi(token));
      }
      std::sort(curr.boundary.begin(),curr.boundary.end());
      //std::cout << "Read " << indices.size() << " indices" << std::endl;
    }
#else
    std::string line,next;           
    for(int i=0;i<n;i++) {
      std::getline(instr,line);
      std::stringstream sstr(line);
      Coordinate x,y;
      sstr >> x >> y;
#if SWAP_GRADE
      Grade grade(y,x);
#else
      Grade grade(x,y);
#endif
      sstr >> next;
      if(next!=";") {
	std::cerr << "Semicolon missing" << std::endl;
      }
      std::vector<index> indices;
      int next_id;
      while(sstr.good()) {
	sstr >> next_id;
	indices.push_back(next_id);
      }
      pre_matrix.push_back(pre_column(i,grade,indices));
    }
#endif
  }


  
  void load_data_into_prematrix(File_reader& reader,
				std::vector<pre_column>& pre_matrix1,
				std::vector<pre_column>& pre_matrix2,
				int& r) {



    std::string next=reader.next_line();
    if(next!="firep") {
      std::cerr << "Keyword 'firep' expected" << std::endl;
      std::exit(1);
    }
    std::string line=reader.next_line();
    line=reader.next_line();
    //std::cout << line << std::endl;

    int t,s;

    {
      line = reader.next_line();
      std::stringstream sstr(line);
      sstr >> t >> s >> r;
    }

    std::cout << "t,s,r=" << t << " " << s << " " << r << std::endl;
    
    load_prematrix_contents(reader,pre_matrix1,t,+1);
    load_prematrix_contents(reader,pre_matrix2,s,-1);

  }

  template<typename GrMat>
    void write_matrix_to_phat(char* filename, 
			      GrMat& matrix) {
    std::ofstream out(filename);

    for(int i=0;i<matrix.get_num_cols();i++) {
      std::vector<index> col;
      matrix.get_col(i,col);
      if(col.size()==0) {
	out << "0 " << std::endl;
      } else {
	index s = col.size();
	out << s-1 << " ";
	for(int j=0;j<s;j++) {
	  out << col[j] << " ";
	}
	out << std::endl;
      }
    }
  }
    

  template<typename GrMat>
    void write_matrix_to_firep(char* filename, 
				  GrMat& matrix1, 
				  GrMat& matrix2) {

    std::ofstream out(filename);
    out << "firep" << std::endl;
    out << "x-coordinate" << std::endl;
    out << "y-coordinate" << std::endl;
    out << matrix1.get_num_cols() << " " << matrix2.get_num_cols() << " " << matrix2.num_rows << std::endl;
    for(int i=0;i<matrix1.get_num_cols();i++) {
      out << std::fixed << std::setprecision(12) << matrix1.grades[i].first_val << " " << std::fixed << std::setprecision(12) << matrix1.grades[i].second_val << " ; ";
      std::vector<index> col;
      matrix1.get_col(i,col);
      for(index j=0;j<col.size();j++) {
	out << col[j] << " ";
      }
      out << std::endl;
    }
    for(int i=0;i<matrix2.get_num_cols();i++) {
      out << matrix2.grades[i].first_val << " " << matrix2.grades[i].second_val << " ; ";
      std::vector<index> col;
      matrix2.get_col(i,col);
      for(index j=0;j<col.size();j++) {
	out << col[j] << " ";
      }
      out << std::endl;
    }
  }
  
  template<typename Matrix>
    void create_matrix_from_firep(char* filename, 
				  Matrix& matrix1, 
				  Matrix& matrix2) {
    std::vector<pre_column> pre_matrix1, pre_matrix2;
    int r;
    std::cout << "Loading data into string..." << std::flush;
    File_reader reader(filename);
    std::cout << "done" << std::endl;
    test_timer1.resume();
    load_data_into_prematrix(reader,pre_matrix1,pre_matrix2,r);
    test_timer1.stop();
    test_timer2.resume();
    Sort_pre_column<pre_column> sort_pre_column;
    std::sort(pre_matrix1.begin(),pre_matrix1.end(),sort_pre_column);
    std::sort(pre_matrix2.begin(),pre_matrix2.end(),sort_pre_column);
    std::vector<index> re_index;
    re_index.resize(pre_matrix2.size());
    test_timer2.stop();
    test_timer3.resume();
    for(index i=0;i<pre_matrix2.size();i++) {
      re_index[pre_matrix2[i].idx]=i;
    }

    for(index i=0;i<pre_matrix1.size();i++) {
      for(index j=0;j<pre_matrix1[i].boundary.size();j++) {
	pre_matrix1[i].boundary[j]=re_index[pre_matrix1[i].boundary[j]];
      }
      std::sort(pre_matrix1[i].boundary.begin(),
		pre_matrix1[i].boundary.end());
    } 
    test_timer3.stop();
    test_timer4.resume();
    {
      int n = pre_matrix1.size();
      NUM_ROWS=pre_matrix2.size();
      matrix1.set_num_cols(n);
      NUM_ROWS=1;
      matrix1.grades.resize(n);

      matrix1.sync();
#if PARALLEL_FOR_LOOPS
#pragma omp parallel for schedule(guided,1)
#endif
      for(int i=0;i<n;i++) {
	pre_column& pcol = pre_matrix1[i];
        matrix1.grades[i]=pcol.grade;
	matrix1.set_col(i,pcol.boundary);
      }
      matrix1.sync();
      matrix1.assign_slave_matrix();
      matrix1.num_rows=pre_matrix2.size();
      matrix1.assign_pivots();
    }
    {
      int n = pre_matrix2.size();
      NUM_ROWS=r;
      matrix2.set_num_cols(n);
      NUM_ROWS=-1;
      matrix2.grades.resize(n);
      matrix2.sync();
#if PARALLEL_FOR_LOOPS
#pragma omp parallel for schedule(guided,1)
#endif
      for(int i=0;i<n;i++) {
	pre_column& pcol = pre_matrix2[i];
        matrix2.grades[i]=pcol.grade;
	matrix2.set_col(i,pcol.boundary);
      }
      matrix2.sync();
      matrix2.assign_slave_matrix();
      matrix2.num_rows=r;
      matrix2.assign_pivots();
    }
    test_timer4.stop();
    test_timer5.resume();
    assign_grade_indices(matrix1,matrix2);
    test_timer5.stop();
    /*
#if !CHUNK_PREPROCESSING
    matrix1.grid_scheduler=Grid_scheduler(matrix1);
    matrix2.grid_scheduler=Grid_scheduler(matrix2);
#endif
    */
    for(index i=0;i<matrix1.num_rows;i++) {
      matrix1.row_grades.push_back(matrix2.grades[i]);
    }
#if SMART_REDUCTION
    matrix1.pq_row.resize(matrix1.num_grades_y);
    matrix2.pq_row.resize(matrix2.num_grades_y);
#endif

#if !NDEBUG
    check_grade_sanity(matrix1);
#endif

  }

  template<typename GradedMatrix>
    void chunk_preprocessing(GradedMatrix& M1, GradedMatrix& M2,
			     GradedMatrix& result1, GradedMatrix& result2) {

    std::vector<int> local_pivots;
    for(index i=0;i<M1.num_rows;i++) {
      local_pivots.push_back(-1);
    }
    std::cout << "Local reduction" << std::endl;
    std::vector<index> global_indices;
    std::vector<char> local_vec;
    local_vec.resize(M1.get_num_cols());
    //test_timer1.start();
    for(index i=0;i<M1.get_num_cols();i++) {
      while(M1.is_local(i)) {
	index p = M1.get_max_index(i);
	index j = local_pivots[p];
	if(j!=-1) {
	  M1.add_to(j,i);
	  gl_no_column_additions++;
	} else {
	  local_pivots[p]=i;
	  local_vec[i]=1;
	  break;
	}
      }
      if(!M1.is_local(i)) {
	global_indices.push_back(i);
	//std::cout << i << " is global" << std::endl;
	local_vec[i]=0;
      }
    }
    //test_timer1.stop();
    std::cout << "Sparsification" << std::endl;
    //test_timer2.start();
    /*
    for(index i=0;i<M1.get_num_cols();i++) {
      if(local_vec[i]) {
	index p = M1.get_max_index(i);
	assert(p!=-1);
	index j = local_pivots[p];
	assert(j==i);
      }
    }
    */
    std::vector<index> col;
    M1.sync();
#if PARALLEL_FOR_LOOPS
#pragma omp parallel for schedule(guided,1), private(col)
#endif
    for(index r=0;r<global_indices.size();r++) {
      col.clear();
      index i = global_indices[r];
      //std::cout << "Fetching " << i << std::endl;
      assert(!local_vec[i]);
      
      //std::cout << "i=" << i << std::endl;
      while(!M1.is_empty(i)) {
	index p = M1.get_max_index(i);
	index j = local_pivots[p];
	if(j!=-1) {
	  assert(local_vec[j]);
	  assert(!local_vec[i]);
	  M1.add_to(j,i);
	  gl_no_column_additions++;
	} else {
	  col.push_back(p);
	  M1.remove_max(i);
	}
      }
      std::reverse(col.begin(),col.end());
      M1.set_col(i,col);
    }
    M1.sync();
    /*
    for(index i=0;i<M1.get_num_cols();i++) {
      std::cout << "XXX ";
      std::vector<index> col;
      M1.get_col(i,col);
      for(index j : col) {
	std::cout << j << " ";
      }
      std::cout << std::endl;
    }
    */

    //test_timer2.stop();
    std::cout << "Build up smaller matrices" << std::endl;
    //test_timer3.start();
    std::vector<int> new_row_index;
    new_row_index.resize(M1.num_rows);
    index row_count=0;
    for(index i=0;i<M1.num_rows;i++) {
      if(local_pivots[i]==-1) {
	new_row_index[i]=row_count++;
      } else {
	new_row_index[i]=-1;
      }
    }
    std::vector<int> new_col_index;
    new_col_index.resize(M1.get_num_cols());
    index col_count=0;
    for(index i=0;i<M1.get_num_cols();i++) {
      if(M1.is_empty(i) || M1.is_local(i)) {
	new_col_index[i]=-1;
      } else {
	new_col_index[i]=col_count++;
      }
    }
    //test_timer3.stop();
    //test_timer4.start();
    NUM_ROWS=row_count;
    result1.set_num_cols(col_count);
    NUM_ROWS=-1;
    for(index i=0;i<M1.get_num_cols();i++) {
      if(new_col_index[i]!=-1) {
	result1.grades.push_back(M1.grades[i]);
	std::vector<index> col,new_col;
	M1.get_col(i,col);
	for(index j=0;j<col.size();j++) {
	  col[j]=new_row_index[col[j]];
	}
	result1.set_col(new_col_index[i],col);
      }
    }
    
    result1.num_rows=row_count;
    NUM_ROWS=M2.num_rows;
    result2.set_num_cols(row_count);
    NUM_ROWS=-1;
    for(index i=0;i<M1.num_rows;i++) {
      if(local_pivots[i]==-1) {
	result1.row_grades.push_back(M1.row_grades[i]);
	result2.grades.push_back(M1.row_grades[i]);
	std::vector<index> col;
	M2.get_col(i,col);
	result2.set_col(new_row_index[i],col);
      } 
    }
    result2.num_rows=M2.num_rows;
    
    result1.assign_slave_matrix();
    result1.assign_pivots();   
    result2.assign_slave_matrix();
    result2.assign_pivots();   


    result1.num_grades_x=M1.num_grades_x;
    result1.num_grades_y=M1.num_grades_y;
    std::copy(M1.x_vals.begin(),M1.x_vals.end(),std::back_inserter(result1.x_vals));
    std::copy(M1.y_vals.begin(),M1.y_vals.end(),std::back_inserter(result1.y_vals));
    

    result2.num_grades_x=M2.num_grades_x;
    result2.num_grades_y=M2.num_grades_y;
    std::copy(M2.x_vals.begin(),M2.x_vals.end(),std::back_inserter(result2.x_vals));
    std::copy(M2.y_vals.begin(),M2.y_vals.end(),std::back_inserter(result2.y_vals));

    std::cout << "Setting up grid scheduler..." << std::flush;

    /*
    result1.grid_scheduler=Grid_scheduler(result1);
    result2.grid_scheduler=Grid_scheduler(result2);
    */    

    std::cout << "done" << std::endl;

#if SMART_REDUCTION
    result1.pq_row.resize(result1.num_grades_y);
    result2.pq_row.resize(result2.num_grades_y);
#endif    

    std::cout << "After chunk reduction, matrix has " << result1.get_num_cols() << " columns and " << result1.num_rows << " rows" << std::endl;
    std::cout << "N' is " << result1.get_num_cols()+ result1.num_rows << std::endl;
    std::cout << "Ratio N'/N: " << (double)(result1.get_num_cols()+ result1.num_rows)/(M1.get_num_cols()+M1.num_rows) << std::endl;

    //test_timer4.stop();

    //result1.print(true,true);
    //result2.print(false,true);

  }



#if SMART_REDUCTION
  template<typename GradedMatrix>
    void min_gens(GradedMatrix& M, GradedMatrix& result) {
    
    index count=0;

    std::vector<Grade> new_grades;
    std::vector<std::vector<index>> new_cols;

    Grid_scheduler& grid = M.grid_scheduler;

    while(! grid.at_end()) {
      
      index_pair new_grade = grid.next_grade();
      
      index x = new_grade.first;
      index y = new_grade.second;
	
      //std::cout << "Min gens " << x << " " << y << std::endl;
	
      PQ& pq = M.pq_row[y];

      index_pair range_xy = grid.index_range_at(x,y);
      
      index start_xy = range_xy.first;
      index end_xy = range_xy.second;
      
      
      //std::cout << "Min gens for " << x << " " << y << " traverses through index range " << start_xy << " " << end_xy << std::endl;
      assert(start_xy<=end_xy);
      //std::cout << "Before adding, pq of row has size " << pq.size() << std::endl;
      for(index i = start_xy;i<end_xy;i++) {
	pq.push(i);
      }
      //std::cout << "After adding, pq of row has size " << pq.size() << std::endl;
      /*
      if(!pq.empty()) {
	std::cout << "Min gens for " << x << " " << y << " has to do work " << pq.size() << " grade cols: " << end_xy-start_xy << std::endl;

      }
      */
      while(!pq.empty()) {
	index i = pq.top();
	//std::cout << "Index " << i << std::endl;
	// Remove duplicates
	while(!pq.empty() && i==pq.top()) {
	  pq.pop();
	}
	assert(M.grades[i].first_index<=x);
	assert(M.grades[i].second_index==y);
	M.reduce_column(i,new_grade,false,true);
	if(!M.is_empty(i) && i>=start_xy&& i<end_xy) {
	  std::vector<index> col;
	  M.get_col(i,col);
	  //std::cout << "NEW MIN GENERATOR Count" << count << " index " << i << " grade " << x << " " << y << std::endl; 
	  new_grades.push_back(M.grades[i]);
	  new_cols.push_back(col);
	}
      }
    }
    NUM_ROWS=M.num_rows;
    result.set_num_cols(new_cols.size());
    NUM_ROWS=-1;
    for(index i=0;i<new_cols.size();i++) {
      result.grades.push_back(new_grades[i]);
      result.set_col(i,new_cols[i]);
    }
    result.num_rows=M.num_rows;
    result.row_grades=M.row_grades;
    
    
#if CLEARING
    // check potential of clearing
    int no_local_pairs=0;
    std::set<int> columns_saved;
    int no_cols_with_dominating_pivots=0;
    for(int i=0;i<result.get_num_cols();i++) {
      if(result.is_empty(i)) {
	continue;
      }
      if(result.is_local(i)) {
	index p = result.get_max_index(i);
	columns_saved.insert(p);
	no_local_pairs++;
      }
      if(result.pivot_is_dominating(i)) {
	no_cols_with_dominating_pivots++;
	index p = result.get_max_index(i);
	columns_saved.insert(p);
	result.clearing_info[p]=i;
      }
      
    }
    std::cout << "Min gens matrix has " << result.get_num_cols() << " columns, out of which " << no_local_pairs << " are local and " << no_cols_with_dominating_pivots << " have a dominating pivot (" << double(no_cols_with_dominating_pivots*100)/result.get_num_cols() << "%)" << std::endl;
    std::cout << "Columns replaced: " << columns_saved.size() << std::endl;
#endif // of CLEARING

  }
#else

  template<typename GradedMatrix>
    void min_gens(GradedMatrix& M, GradedMatrix& result) {
    
    std::vector<Grade> new_grades;
    std::vector<std::vector<index>> new_cols;

    Grid_scheduler& grid = M.grid_scheduler;

    while(! grid.at_end()) {
      
      index_pair new_grade = grid.next_grade();
      
      index x = new_grade.first;
      index y = new_grade.second;
      
      index start_xy,end_xy,start_xy_of_grade;
      
      grid.extended_index_range_at(x,y,start_xy,end_xy,start_xy_of_grade);
	
      //std::cout << "Min gens for " << x << " " << y << " traverses through index range " << start_xy << " (" << start_xy_of_grade << ") " << end_xy << std::endl;
      assert(start_xy<=start_xy_of_grade);
      assert(start_xy_of_grade<=end_xy);
      for(index i = start_xy;i<end_xy;i++) {
	M.reduce_column(i,new_grade);
	if(!M.is_empty(i) && i>=start_xy_of_grade&& i<end_xy) {
	  std::vector<index> col;
	  M.get_col(i,col);
	  //std::cout << "NEW MIN GENERATOR Count" << count << " index " << i << " grade " << x << " " << y << std::endl; 
	  new_grades.push_back(M.grades[i]);
	  new_cols.push_back(col);
	}
      }
    }
    NUM_ROWS=M.num_rows;
    result.set_num_cols(new_cols.size());
    NUM_ROWS=-1;
    for(index i=0;i<new_cols.size();i++) {
      result.grades.push_back(new_grades[i]);
      result.set_col(i,new_cols[i]);
    }
    result.num_rows=M.num_rows; 
  }

#endif

#if SMART_REDUCTION


  // If no clearing is used, the argument mingens is ignored
  template<typename GradedMatrix>
    void ker_basis(GradedMatrix& M, GradedMatrix& result, GradedMatrix& mingens) {
    
    std::vector<Grade> new_grades;
    std::vector<std::vector<index>> new_cols;
    
    std::vector<bool> indices_in_kernel;
    indices_in_kernel.resize(M.get_num_cols());
    for(index i=0;i<M.get_num_cols();i++) {
      indices_in_kernel[i]=false;
    }

    //std::set<index> indices_in_kernel;

    Grid_scheduler& grid = M.grid_scheduler;
    
    while(! grid.at_end()) {
      
      index_pair new_grade = grid.next_grade();
      
      index x = new_grade.first;
      index y = new_grade.second;

      //std::cout << "At grade " << x << " " << y << std::endl;
      
      PQ& pq = M.pq_row[y];
      
      index_pair range_xy = grid.index_range_at(x,y);
      
      index start_xy = range_xy.first;
      index end_xy = range_xy.second;

      assert(start_xy<=end_xy);

      for(index i = start_xy;i<end_xy;i++) {
	pq.push(i);
      }
      //std::cout << "After adding, pq of row has size " << pq.size() << std::endl;
      //test_timer1.resume();
      while(!pq.empty()) {
	index i = pq.top();
	//std::cout << "next top" << i << std::endl;
	// Remove duplicates
	while(!pq.empty() && i==pq.top()) {
	  pq.pop();
	}
	assert(M.grades[i].first_index<=x);
	assert(M.grades[i].second_index==y);
#if CLEARING
	if(mingens.clearing_info.count(i)!=0) {
	  std::vector<index> col;
	  mingens.get_col(mingens.clearing_info[i],col);
	  new_cols.push_back(col);
	  new_grades.push_back(Grade(x,y,M.x_vals[x],M.y_vals[y]));
	  indices_in_kernel[i]=true;;
	  //indices_in_kernel.insert(i);
	  M.clear(i);
	}
#endif
	//test_timer3.resume();
	M.reduce_column(i,new_grade,true,true);
	//test_timer3.stop();
	//test_timer2.resume();
	if(!indices_in_kernel[i] && M.is_empty(i)) {
	  //if(M.is_empty(i) && indices_in_kernel.count(i)==0) {
	  //std::cout << "NEW KERNEL ELEMENT " << i << " Count: " << count << " Grade " << x << " " << y << std::endl;
	  {
	    std::vector<index> col;
	    M.slave.get_col(i,col);
	    new_cols.push_back(col);
	    new_grades.push_back(Grade(x,y,M.x_vals[x],M.y_vals[y]));
	    indices_in_kernel[i]=true;
	    //indices_in_kernel.insert(i);
	  }
	}
	//test_timer2.stop();
      }
      //test_timer1.stop();
    }
    NUM_ROWS=M.get_num_cols();
    result.set_num_cols(new_cols.size());
    NUM_ROWS=-1;
    for(index i=0;i<new_cols.size();i++) {
      result.grades.push_back(new_grades[i]);
      result.set_col(i,new_cols[i]);
    }
    result.num_rows=M.get_num_cols();
    result.slave.set_num_cols(result.get_num_cols());
    result.assign_pivots();
    for(index i=0;i<result.get_num_cols();i++) {
      result.pivots[result.get_max_index(i)]=i;
      std::vector<index> slave_col;
      slave_col.push_back(i);
      result.slave.set_col(i,slave_col);
    }
    result.num_grades_x=M.num_grades_x;
    result.num_grades_y=M.num_grades_y;
    std::copy(M.x_vals.begin(),M.x_vals.end(),std::back_inserter(result.x_vals));
    std::copy(M.y_vals.begin(),M.y_vals.end(),std::back_inserter(result.y_vals));
  }


#else

  template<typename GradedMatrix>
    void ker_basis(GradedMatrix& M, GradedMatrix& result,GradedMatrix& mingens) {

    std::vector<Grade> new_grades;
    std::vector<std::vector<index>> new_cols;

    std::vector<bool> indices_in_kernel;
    indices_in_kernel.resize(M.get_num_cols());
    for(index i=0;i<M.get_num_cols();i++) {
      indices_in_kernel[i]=false;
    }

    Grid_scheduler& grid = M.grid_scheduler;
    
    while(! grid.at_end()) {
      
      index_pair new_grade = grid.next_grade();
      
      index x = new_grade.first;
      index y = new_grade.second;
      
      index start_xy,end_xy,start_xy_of_grade;
      
      grid.extended_index_range_at(x,y,start_xy,end_xy,start_xy_of_grade);
      
      assert(start_xy<=start_xy_of_grade);
      assert(start_xy_of_grade<=end_xy);
      
      //test_timer3.resume();
      for(index i = start_xy;i<end_xy;i++) {
	M.reduce_column(i,new_grade,true);
	if(!indices_in_kernel[i] && M.is_empty(i)) {
	  //std::cout << "NEW KERNEL ELEMENT " << i << " Count: " << count << " Grade " << x << " " << y << std::endl;
	  std::vector<index> col;
	  M.slave.get_col(i,col);
	  new_cols.push_back(col);
	  new_grades.push_back(Grade(x,y,M.x_vals[x],M.y_vals[y]));
	  indices_in_kernel[i]=true;
	}
      }
      //test_timer3.stop();
      
    }
    NUM_ROWS=M.get_num_cols();
    result.set_num_cols(new_cols.size());
    NUM_ROWS=-1;
    for(index i=0;i<new_cols.size();i++) {
      result.grades.push_back(new_grades[i]);
      result.set_col(i,new_cols[i]);
    }
    result.num_rows=M.get_num_cols();
    result.slave.set_num_cols(result.get_num_cols());
    result.assign_pivots();
    for(index i=0;i<result.get_num_cols();i++) {
      result.pivots[result.get_max_index(i)]=i;
      std::vector<index> slave_col;
      slave_col.push_back(i);
      result.slave.set_col(i,slave_col);
    }
    result.num_grades_x=M.num_grades_x;
    result.num_grades_y=M.num_grades_y;
    std::copy(M.x_vals.begin(),M.x_vals.end(),std::back_inserter(result.x_vals));
    std::copy(M.y_vals.begin(),M.y_vals.end(),std::back_inserter(result.y_vals));
  }


#endif

  template<typename GradedMatrix>
    void reparameterize(GradedMatrix& cols, GradedMatrix& ker, GradedMatrix& result) {
    index ker_cols = ker.get_num_cols();
    NUM_ROWS=ker.num_rows;
    ker.set_num_cols(ker_cols + cols.get_num_cols());
    NUM_ROWS=-1;
    NUM_ROWS=ker_cols;
    ker.slave.set_num_cols(ker_cols+cols.get_num_cols());
    NUM_ROWS=-1;
    std::vector<index> empty_vec;
    for(index i = ker_cols;i<ker.get_num_cols();i++) {
      ker.slave.set_col(i,empty_vec);
    }
    NUM_ROWS=ker_cols;
    result.set_num_cols(cols.get_num_cols());
    NUM_ROWS=-1;
    
    /*
    std::cout << "Pivots:" << std::endl;
    for(index i=0;i<ker.num_rows;i++) {
      std::cout << i << " -> " << ker.pivots[i] << std::endl;
    }
    */


    std::copy(cols.grades.begin(),cols.grades.end(),std::back_inserter(result.grades));

    //test_timer4.start();
    index_pair dummy;
    cols.sync();
    ker.sync();
    result.sync();
#if PARALLEL_FOR_LOOPS
#pragma omp parallel for schedule(guided,1)
#endif
    for(index i=0;i<cols.get_num_cols();i++) {
      //std::cout << "index" << i << std::endl;
      std::vector<index> col;
      cols.get_col(i,col);
      ker.set_col(ker_cols+i,col);
      //std::cout << "reduce" << std::endl;
      ker.reduce_column(ker_cols+i,dummy,true);
      //std::cout << "done" << std::endl;
      assert(ker.is_empty(ker_cols+i));
      std::vector<index> new_col;
      ker.slave.get_col(ker_cols+i,new_col);
      result.set_col(i,new_col);
    }
    cols.sync();
    ker.sync();
    result.sync();
    //test_timer4.stop();
    // Assign row grades
    result.num_rows=ker_cols;
    std::copy(ker.grades.begin(),ker.grades.end(),std::back_inserter(result.row_grades));
    
    // We should probably trim here, or later for the minimized matrix
    result.num_grades_x=ker.num_grades_x;
    result.num_grades_y=ker.num_grades_y;
    std::copy(ker.x_vals.begin(),ker.x_vals.end(),std::back_inserter(result.x_vals));
    std::copy(ker.y_vals.begin(),ker.y_vals.end(),std::back_inserter(result.y_vals));
#if !NDEBUG
    check_grade_sanity(result);
#endif
  }

  template<typename Rep>
    void convert_to_vec_vec(GradedMatrix<Rep>& M, GradedMatrix<phat::vector_vector>& result) {

    result.num_rows=M.num_rows;
    std::copy(M.pivots.begin(),M.pivots.end(),std::back_inserter(result.pivots));
    std::copy(M.grades.begin(),M.grades.end(),std::back_inserter(result.grades));
    std::copy(M.row_grades.begin(),M.row_grades.end(),std::back_inserter(result.row_grades));
    result.set_num_cols(M.get_num_cols());
    for(index i=0;i<result.get_num_cols();i++) {
      std::vector<index> col;
      M.get_col(i,col);
      result.set_col(i,col);
    }

  }

  // GradedMatrix must be vector_vector here!
  template<typename GradedMatrix>
    bool contains(GradedMatrix& M, index i, index m) {
    std::vector<index> col;
    M.get_col(i,col);
    return std::binary_search(col.begin(),col.end(),m);
  }



  template<typename GradedMatrix_>
    void minimize(GradedMatrix_& M,GradedMatrix_& result) {


    //test_timer1.start();
#if LAZY_MINIMIZATION
    GradedMatrix_& VVM=M;
#else    
    GradedMatrix<phat::vector_vector> VVM;
    
    convert_to_vec_vec(M,VVM);

#endif

    //test_timer1.stop();

    typedef std::vector<index> Column;


    VVM.assign_pivots();

    std::set<index> rows_to_delete;

    std::vector<index> cols_to_keep;
    std::vector<Grade> col_grades;
    for(index i=0;i<VVM.get_num_cols();i++) {
      while(! VVM.is_empty(i)) {
	
	index col_grade_x = VVM.grades[i].first_index;
	index col_grade_y = VVM.grades[i].second_index;
	index p = VVM.get_max_index(i);
	index row_grade_x = VVM.row_grades[p].first_index;
	index row_grade_y = VVM.row_grades[p].second_index;
	
	//std::cout << "Grades: " << col_grade_x << " " << col_grade_y <<  " -- " <<  row_grade_x << " " << row_grade_y << std::endl;

	//std::cout << "Pivot is " <<  p << std::endl;
	if(col_grade_x!=row_grade_x || col_grade_y!=row_grade_y) {
	  cols_to_keep.push_back(i);
	  col_grades.push_back(VVM.grades[i]);
	  break;
	}
	if(VVM.pivots[p]==-1) {
	  //std::cout << "Found removable pair " << p << "( " <<  VVM.row_grades[p].first_index << " " << VVM.row_grades[p].second_index << " )" << " " << i << "( " << VVM.grades[i].first_index << " " << VVM.grades[i].second_index << " )" << std::endl;
	  rows_to_delete.insert(p);
#if !LAZY_MINIMIZATION
	  for(index j=i+1;j<VVM.get_num_cols();j++) {
	    if(VVM.contains(j,p)) {
	      VVM.add_to(i,j);
	      gl_no_column_additions++;
	    }
	  }
#endif
	  VVM.pivots[p]=i;
	  break;
	} else {
	  //std::cout << "add a column" << VVM.pivots[p] << std::endl;
	  VVM.add_to(VVM.pivots[p],i);
	  gl_no_column_additions++;
	}
      }
      assert(! VVM.is_empty(i));
    }
    std::vector<Column> new_cols;
    new_cols.resize(cols_to_keep.size());


#if LAZY_MINIMIZATION
    //std::cout << "HERE I AM " << std::endl;
    //test_timer2.start();
    VVM.sync();
#if PARALLEL_FOR_LOOPS
#pragma omp parallel for schedule(guided,1)
#endif
    for(index k=0;k<cols_to_keep.size();k++) {
      Column& col=new_cols[k];
      index i = cols_to_keep[k];
      while(!VVM.is_empty(i)) {
	index p = VVM.get_max_index(i);
	if(VVM.pivots[p]==-1) {
	  col.push_back(p);
	  VVM.remove_max(i);
	} else {
	  VVM.add_to(VVM.pivots[p],i);
	  gl_no_column_additions++;
	}
      }
      std::reverse(col.begin(),col.end());
    }
    VVM.sync();
    //test_timer2.stop();
#else
    for(index i=0;i<cols_to_keep.size();i++) {
      Column& col=new_cols[i];
      VVM.get_col(cols_to_keep[i],new_cols[i]);
    }
#endif
    
    index nr = VVM.num_rows;
    //std::cout << "Number of rows of VVM: " << nr << std::endl;
    //std::cout << "Have removed " << rows_to_delete.size() << std::endl;
    //std::cout << "Kept " << cols_to_keep.size() << " columns" << std::endl;
    index count=0;
    std::unordered_map<index,index> index_map;
    std::vector<Grade> res_row_grades;
    // Build re-indexing map
    for(index i=0;i<nr;i++) {
      if(rows_to_delete.count(i)==0) {
	//std::cout << "reindex " << i << " -> " << count << std::endl; 
	index_map[i]=count++;
	res_row_grades.push_back(VVM.row_grades[i]);
      }
    }
    // Update columns using re-index map
    
    for(index i=0;i<cols_to_keep.size();i++) {
      //std::cout << "Iterating throgh column " << cols_to_keep[i] << std::endl;
      Column& col = new_cols[i];
      for(int j=0;j<col.size();j++) {
	//std::cout << "Lookup " << col[j] << ", count=" << index_map.count(col[j]) << std::endl;
	assert(index_map.count(col[j]));
	col[j]=index_map[col[j]];
      }
    }
    // Now build up the result
    NUM_ROWS=index_map.size();
    result.set_num_cols(cols_to_keep.size());
    NUM_ROWS=-1;
    result.num_rows=index_map.size();
    result.grades = col_grades;
    result.row_grades=res_row_grades;
    for(int i=0;i<cols_to_keep.size();i++) {
      result.set_col(i,new_cols[i]);
    }
    // We should trim the grades here
    result.num_grades_x = M.num_grades_x;
    result.num_grades_y = M.num_grades_y;
    std::copy(M.x_vals.begin(),M.x_vals.end(),std::back_inserter(result.x_vals));
    std::copy(M.y_vals.begin(),M.y_vals.end(),std::back_inserter(result.y_vals));
  }
  
  template<typename GradedMatrix>
    void swap_grades(GradedMatrix& M) {
    for(index i=0;i<M.get_num_cols();i++) {
      Grade& gr = M.grades[i];
      std::swap(gr.first_index,gr.second_index);
      std::swap(gr.first_val,gr.second_val);
      
    }
  }

  struct grade_colex_sort {

    template<typename T>
    bool operator()(T& t1, T& t2) {
      if(t1.first.second_index<t2.first.second_index) {
	return true;
      }
      if(t1.first.second_index>t2.first.second_index) {
	return false;
      }
      if(t1.first.first_index < t2.first.first_index) {
	return true;
      }
      if(t1.first.first_index > t2.first.first_index) {
	return false;
      }
      // flip the stable sort
      return t1.second < t2.second;
    }

  };

  template<typename GradedMatrix>
    void convert_to_colex_order(GradedMatrix& M) {
    
    std::vector<std::pair<Grade,index> > row_info;

    for(index i=0;i<M.num_rows;i++) {
      row_info.push_back(std::make_pair(M.row_grades[i],i));
    }
    std::sort(row_info.begin(),row_info.end(),grade_colex_sort());
    
    std::unordered_map<index,index> re_index_rows;

    for(index i=0;i<row_info.size();i++) {
      re_index_rows[row_info[i].second]=i;
    }

    for(index i=0;i<M.num_rows;i++) {
      M.row_grades[i]=row_info[i].first;
    }
    for(index i=0;i<M.get_num_cols();i++) {
      std::vector<index> col;
      M.get_col(i,col);
      for(index j=0;j<col.size();j++) {
	col[j]=re_index_rows[col[j]];
      }
      std::sort(col.begin(),col.end());
      M.set_col(i,col);
    }

    // Now columns
    std::vector<std::pair<Grade,index> > column_info;
    for(index i=0;i<M.get_num_cols();i++) {
      column_info.push_back(std::make_pair(M.grades[i],i));
    }
    std::sort(column_info.begin(),column_info.end(),grade_colex_sort());

    std::vector<std::vector<index>> cols;
    cols.resize(M.get_num_cols());

    for(index i=0;i<M.get_num_cols();i++) {
      M.grades[i]=column_info[i].first;
      M.get_col(column_info[i].second,cols[i]);
    }
    for(index i=0;i<M.get_num_cols();i++) {
      M.set_col(i,cols[i]);
    }
  }

  template<typename GradedMatrix,typename ofstr>
    void print_in_rivet_format(GradedMatrix& M,ofstr& out) {

    convert_to_colex_order(M);
    
    out << "x-grades" << std::endl;
    for(index i=0;i<M.num_grades_x;i++) {
#if USE_DOUBLE
      out << rat_representation(M.x_vals[i]) << std::endl;
#else
      Coordinate& c = M.x_vals[i];
      boost::multiprecision::cpp_rational to_print = c.convert_to<boost::multiprecision::cpp_rational>();
      out << to_print << std::endl;
#endif
      
    }
    out << std::endl;
    out << "y-grades" << std::endl;
    for(int i=0;i<M.num_grades_y;i++) {
#if USE_DOUBLE
      out << rat_representation(M.y_vals[i]) << std::endl;
#else
      Coordinate& c = M.y_vals[i];
      boost::multiprecision::cpp_rational to_print = c.convert_to<boost::multiprecision::cpp_rational>();
      out << to_print << std::endl;
#endif
    }
    out << std::endl;
    
    out << "MINIMAL PRESENTATION:" << std::endl;
    out << "Number of rows:" << M.num_rows << std::endl;
    out << "Row bigrades:" << std::endl;
    out << "| ";
    for(index i=0;i<M.num_rows;i++) {
      out<<"("<<M.row_grades[i].first_index<<","<<M.row_grades[i].second_index<<") ";
    }
    out << "|" << std::endl;
    out << "Number of columns:" << M.get_num_cols() << std::endl;
    out << "Column bigrades:" << std::endl;
    out << "| ";
    for(index i=0;i<M.get_num_cols();i++) {
      out<<"("<<M.grades[i].first_index<<","<<M.grades[i].second_index<<") ";
    }
    out << "|" << std::endl;
    for(index i=0;i<M.get_num_cols();i++) {
      std::vector<index> col;
      M.get_col(i,col);
      for(index j=0;j<col.size();j++) {
	out << col[j] << " ";
      }
      out << std::endl;
    }
  }


  
}//of namespace phat
