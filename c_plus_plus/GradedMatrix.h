#include<phat/boundary_matrix.h>

namespace phat {

  typedef std::pair<index,index> index_pair;

  struct Grade {
    double first_val;
    index first_index;
    double second_val;
    index second_index;
    Grade(double x, double y) : first_val(x), second_val(y) {}
    Grade(const Grade& other) : first_val(other.first_val), first_index(other.first_index),
				second_val(other.second_val), second_index(other.second_index) {}
    bool operator== (const Grade& other) {
      return this->first_val==other.first_val && this->second_val==other.second_val;
    }
  };

  struct pre_column {
    Grade grade;
    std::vector<index> boundary;
    pre_column(Grade& grade, std::vector<index>& boundary)
      : grade(grade), boundary(boundary) {
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

  template<typename Representation=vector_heap>
    class GradedMatrix  : public boundary_matrix<Representation> {
    
    public:

    index num_grades_x;
    index num_grades_y;

    index num_rows;
    
    std::map<index_pair,index> start_index_of_pair;

    std::vector<Grade> grades;

    std::vector<Grade> row_grades;
    
    std::vector<index> pivots;

    boundary_matrix<Representation> slave;

    void print(bool print_indices=true) {
      std::cout << "Number of columns: " << this->get_num_cols() << std::endl;
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
    }

    void reduce_column(index i, bool use_slave=false) {
      
      std::cout << "Reduce " << i << std::endl;

      std::vector<index> col;
      this->get_col(i,col);
      for(int i=0;i<col.size();i++) {
	std::cout << col[i] << " ";
      }
      std::cout << std::endl;

      index p = this->get_max_index(i);
      
      std::cout << "p=" << p << std::endl;

      while(p!=-1 && pivots[p]!=-1 && pivots[p]<i) {
	index k = pivots[p];

	this->add_to(k,i);
	
	if(use_slave) {
	  slave.add_to(k,i);
	}
	p=this->get_max_index(i);
      }
      if(p!=-1 && (pivots[p]==-1 || pivots[p]>i)) {
	pivots[p]=i;
      }
    }

    

    void assign_pivots() {
      for(index i=0;i<this->num_rows;i++) {
	pivots.push_back(-1);
      }
    }

    // Assumes that columns are already added to matrix in lex order
    void assign_grade_indices() {
      int n=this->get_num_cols();
      std::cout << "Assgin grade indices with " << n << " columns" << std::endl;
      if(n==0) {
	return;
      }
      
      std::map<double,index> val_to_index_x, val_to_index_y;
      std::vector<double> x_vals,y_vals;

      for(int i=0;i<n;i++) {
	x_vals.push_back(grades[i].first_val);
	y_vals.push_back(grades[i].second_val);
      }
      std::sort(x_vals.begin(),x_vals.end());
      auto last_x = std::unique(x_vals.begin(),x_vals.end());
      x_vals.erase(last_x,x_vals.end());
      std::sort(y_vals.begin(),y_vals.end());
      auto last_y = std::unique(y_vals.begin(),y_vals.end());
      y_vals.erase(last_y,y_vals.end());

      std::cout << "Found " << x_vals.size() << " different x-values and " << y_vals.size() << " different y-values" << std::endl;
      this->num_grades_x = x_vals.size();
      this->num_grades_y = y_vals.size();
      
      for(index i=0;i<x_vals.size();i++) {
	val_to_index_x[x_vals[i]]=i;
      }
      for(index i=0;i<y_vals.size();i++) {
	val_to_index_y[y_vals[i]]=i;
      }

      for(int i=0;i<n;i++) {
	grades[i].first_index=val_to_index_x[grades[i].first_val];
	grades[i].second_index=val_to_index_y[grades[i].second_val];
      }

      index run_x=0;
      index run_y=0;

      index curr_x,curr_y;

      for(int i=0;i<=n;i++) {
	
	if(i<n) {
	  curr_x=grades[i].first_index;
	  curr_y=grades[i].second_index;
	} else {
	  curr_x=this->num_grades_x-1;
	  curr_y=this->num_grades_y-1;
	}

	//assert(run_x<=curr_x);
	//assert(run_y<=curr_y);

	while(run_x <= curr_x && run_y <= curr_y) {
	  std::cout << "Set index of " << run_x << " " << run_y << " to " << i << std::endl;
	  start_index_of_pair[std::make_pair(run_x,run_y)]=i;
	  run_x++;
	  if(run_x==this->num_grades_x) {
	    run_y++;
	    run_x=0;
	  }
	}
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
  
  template<typename Instream, typename Matrix>
    void create_matrix_from_firep(Instream& instr, 
				  Matrix& matrix1, 
				  Matrix& matrix2) {
    
    std::vector<pre_column> pre_matrix1, pre_matrix2;

    std::string next;
    std::getline(instr,next);
    if(next!="firep") {
      std::cerr << "Keyword 'firep' expected" << std::endl;
      return;
    }
    std::string line;
    // Read over label info
    std::getline(instr,line);
    std::cout << line << std::endl;
    std::getline(instr,line);
    std::cout << line << std::endl;

    int t,s,r;

    {
      std::getline(instr,line);
      std::stringstream sstr(line);
      sstr >> t >> s >> r;
    }

    std::cout << "t,s,r=" << t << " " << s << " " << r << std::endl;
    
    for(int i=0;i<t;i++) {
      std::getline(instr,line);
      std::stringstream sstr(line);
      double x,y;
      sstr >> x >> y;
      Grade grade(x,y);
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
      pre_matrix1.push_back(pre_column(grade,indices));
    }
    for(int i=0;i<s;i++) {
      std::getline(instr,line);
      std::stringstream sstr(line);
      double x,y;
      sstr >> x >> y;
      Grade grade(x,y);
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
      pre_matrix2.push_back(pre_column(grade,indices));
    }
    Sort_pre_column<pre_column> sort_pre_column;
    std::sort(pre_matrix1.begin(),pre_matrix1.end(),sort_pre_column);
    std::sort(pre_matrix2.begin(),pre_matrix2.end(),sort_pre_column);

    {
      int n = pre_matrix1.size();
      matrix1.set_num_cols(n);
      
      for(int i=0;i<n;i++) {
	pre_column& pcol = pre_matrix1[i];
        matrix1.grades.push_back(pcol.grade);
	matrix1.set_col(i,pcol.boundary);
      }
      matrix1.assign_grade_indices();
      matrix1.assign_slave_matrix();
      matrix1.num_rows=s;
      matrix1.assign_pivots();
    }
    {
      int n = pre_matrix2.size();
      matrix2.set_num_cols(n);
      
      for(int i=0;i<n;i++) {
	pre_column& pcol = pre_matrix2[i];
        matrix2.grades.push_back(pcol.grade);
	matrix2.set_col(i,pcol.boundary);
      }
      matrix2.assign_grade_indices();
      matrix2.assign_slave_matrix();
      matrix2.num_rows=r;
      matrix2.assign_pivots();
    }

  }

  template<typename GradedMatrix>
    void min_gens(GradedMatrix& M, GradedMatrix& result) {
    
    index count=0;

    for(index x = 0; x < M.num_grades_x;x++) {
      for(index y = 0; y < M.num_grades_y;y++) {
	
	index start_xy = M.start_index_of_pair[std::make_pair(0,y)];
	index end_xy;
	if(x<M.num_grades_x-1) {
	  end_xy = M.start_index_of_pair[std::make_pair(x+1,y)];
	} else {
	  end_xy = M.start_index_of_pair[std::make_pair(0,y+1)];
	}
	if(x==M.num_grades_x-1 && y==M.num_grades_y-1) {
	  end_xy=M.get_num_cols();
	}
	std::cout << "Min gens for " << x << " " << y << " traverses through index range " << start_xy << " " << end_xy << std::endl;
	index start_xy_of_grade = M.start_index_of_pair[std::make_pair(x,y)];
	assert(start_xy<=start_xy_of_grade);
	assert(start_xy_of_grade<=end_xy);
	for(index i = start_xy;i<end_xy;i++) {
	  M.reduce_column(i);
	  if(!M.is_empty(i) && i>=start_xy_of_grade) {
	    std::vector<index> col;
	    M.get_col(i,col);
	    result.set_num_cols(count+1);
	    result.set_col(count++,col);
	    result.grades.push_back(M.grades[i]);
	  }
	}
      }
    }
  }

  template<typename GradedMatrix>
    void ker_basis(GradedMatrix& M, GradedMatrix& result) {
    
    index count=0;

    std::set<index> indices_in_matrix;

    for(index x = 0; x < M.num_grades_x;x++) {
      for(index y = 0; y < M.num_grades_y;y++) {
	
	index start_xy = M.start_index_of_pair[std::make_pair(0,y)];
	index end_xy;
	if(x<M.num_grades_x-1) {
	  end_xy = M.start_index_of_pair[std::make_pair(x+1,y)];
	} else {
	  end_xy = M.start_index_of_pair[std::make_pair(0,y+1)];
	}
	if(x==M.num_grades_x-1 && y==M.num_grades_y-1) {
	  end_xy=M.get_num_cols();
	}
	std::cout << "Min gens for " << x << " " << y << " traverses through index range " << start_xy << " " << end_xy << std::endl;
	index start_xy_of_grade = M.start_index_of_pair[std::make_pair(x,y)];
	assert(start_xy<=start_xy_of_grade);
	assert(start_xy_of_grade<=end_xy);
	for(index i = start_xy;i<end_xy;i++) {
	  M.reduce_column(i,true);
	  if(M.is_empty(i) && indices_in_matrix.count(i)==0) {
	    std::vector<index> col;
	    M.slave.get_col(i,col);
	    result.set_num_cols(count+1);
	    result.set_col(count++,col);
	    result.grades.push_back(Grade(x,y));
	  }
	}
      }
    }
    result.num_rows=count;
    result.slave.set_num_cols(result.get_num_cols());
    for(index i=0;i<result.get_num_cols();i++) {
      result.pivots.push_back(result.get_max_index(i));
      std::vector<index> slave_col;
      slave_col.push_back(i);
      result.slave.set_col(i,slave_col);
    }
  }
  template<typename GradedMatrix>
    void reparameterize(GradedMatrix& cols, GradedMatrix& ker, GradedMatrix& result) {
    index ker_cols = ker.get_num_cols();
    ker.set_num_cols(ker_cols + cols.get_num_cols());
    ker.slave.set_num_cols(ker_cols+cols.get_num_cols());
    std::vector<index> empty_vec;
    for(index i = ker_cols;i<ker.get_num_cols();i++) {
      ker.slave.set_col(i,empty_vec);
    }
    result.set_num_cols(cols.get_num_cols());
    

    for(index i=0;i<cols.get_num_cols();i++) {
      std::cout << "index" << i << std::endl;
      std::vector<index> col;
      cols.get_col(i,col);
      ker.set_col(ker_cols+i,col);
      std::cout << "reduce" << std::endl;
      ker.reduce_column(ker_cols+i,true);
      std::cout << "done" << std::endl;
      assert(ker.is_empty(ker_cols+i));
      std::vector<index> new_col;
      ker.slave.get_col(ker_cols+i,new_col);
      result.set_col(i,new_col);
      result.grades.push_back(cols.grades[i]);
    }
    // Assign row grades
    result.num_rows=ker_cols;
    std::copy(ker.grades.begin(),ker.grades.end(),std::back_inserter(result.row_grades));
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
    return std::bindary_search(col.begin(),col.end(),m);
  }

  template<typename GradedMatrix_>
    void minimize(GradedMatrix_& M) {
    
    GradedMatrix<phat::vector_vector> VVM;
    convert_to_vec_vec(M,VVM);

    for(index i=0;i<M.get_num_cols();i++) {
      if(M.is_empty(i)) {
	continue;
      }
      index col_grade_x = M.grades[i].first_index;
      index col_grade_y = M.grades[i].second_index;
      index p = M.get_max_index(i);
      index row_grade_x = M.row_grades[p].first_index;
      index row_grade_y = M.row_grades[p].second_index;
      if(col_grade_x==row_grade_x && col_grade_y==row_grade_y) {
	
      }
    }
    
  }
  
  
}//of namespace phat
