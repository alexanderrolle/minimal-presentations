#include<phat/boundary_matrix.h>
#include<algorithm>

namespace phat {

  typedef std::pair<index,index> index_pair;

  struct Grade {
    double first_val;
    index first_index;
    double second_val;
    index second_index;
    Grade(double x, double y) : first_val(x), second_val(y) {}
    Grade(index ind_x, index ind_y, double val_x, double val_y) : first_val(val_x), first_index(ind_x), second_val(val_y), second_index(ind_y) {}
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

  template<typename Representation=vector_heap>
    class GradedMatrix  : public boundary_matrix<Representation> {
    
    public:

    std::vector<double> x_vals,y_vals;

    index num_grades_x;
    index num_grades_y;

    index num_rows;
    
    std::map<index_pair,index> start_index_of_pair;

    std::vector<Grade> grades;

    std::vector<Grade> row_grades;
    
    std::vector<index> pivots;

    boundary_matrix<Representation> slave;

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

    void reduce_column(index i, bool use_slave=false) {
      
      //std::cout << "Reduce " << i << std::endl;

      /*
      std::cout << "Pivots vector: " << std::endl;
      for(int i=0;i<this->num_rows;i++) {
	std::cout << i << " -> " << this->pivots[i] << std::endl;
      }
      */

      /*
      std::vector<index> col;
      this->get_col(i,col);
      for(int i=0;i<col.size();i++) {
	std::cout << col[i] << " ";
      }
      std::cout << std::endl;
      */

      index p = this->get_max_index(i);
      
      //std::cout << "p=" << p << std::endl;

      while(p!=-1 && pivots[p]!=-1 && pivots[p]<i) {
	index k = pivots[p];

	//std::cout << "Adding column " << k << std::endl;

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


  // Assumes that columns are already added to matrix in lex order
  template<typename GradedMatrix>
    void assign_grade_indices(GradedMatrix& M1, GradedMatrix& M2) {
    int n1=M1.get_num_cols();
    int n2=M2.get_num_cols();
    //std::cout << "Assgin grade indices with " << n1 << ", " << n2 << " columns" << std::endl;
    if(n1==0 && n2==0) {
      return;
    }
    
    std::map<double,index> val_to_index_x, val_to_index_y;
    
    std::vector<double> x_vals,y_vals;
    
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
    
    index run_x=0;
    index run_y=0;
    
    index curr_x,curr_y;
    
    for(int i=0;i<=n1;i++) {
      
      if(i<n1) {
	curr_x=M1.grades[i].first_index;
	curr_y=M1.grades[i].second_index;
      } else {
	curr_x=M1.num_grades_x-1;
	curr_y=M1.num_grades_y-1;
      }
      
	//assert(run_x<=curr_x);
	//assert(run_y<=curr_y);
      
      while(run_y<curr_y || (run_y==curr_y && run_x <= curr_x)) {
	//std::cout << "M1: Set index of " << run_x << " " << run_y << " to " << i << std::endl;
	//std::cout << "Curr: " << curr_x << " " << curr_y << std::endl;
	M1.start_index_of_pair[std::make_pair(run_x,run_y)]=i;
	run_x++;
	if(run_x==M1.num_grades_x) {
	  run_y++;
	  run_x=0;
	}
      }
    }

    run_x=run_y=0;

    for(int i=0;i<=n2;i++) {
      
      if(i<n2) {
	curr_x=M2.grades[i].first_index;
	curr_y=M2.grades[i].second_index;
      } else {
	curr_x=M2.num_grades_x-1;
	curr_y=M2.num_grades_y-1;
      }
      
      //std::cout << "M2: Curr: " << curr_x << " " << curr_y << std::endl;

	//assert(run_x<=curr_x);
	//assert(run_y<=curr_y);
      
      while(run_y<curr_y || (run_y==curr_y && run_x <= curr_x)) {
	//std::cout << "M2: Set index of " << run_x << " " << run_y << " to " << i << std::endl;
	M2.start_index_of_pair[std::make_pair(run_x,run_y)]=i;
	run_x++;
	if(run_x==M2.num_grades_x) {
	  run_y++;
	  run_x=0;
	}
      }
    }

  }
  
  template<typename Instream, typename Matrix>
    void create_matrix_from_firep(Instream& instr, 
				  Matrix& matrix1, 
				  Matrix& matrix2) {
    
    std::vector<pre_column> pre_matrix1, pre_matrix2;

    std::string next;
    std::getline(instr,next);
    if(next!="firep") {
      std::cerr << "Keyword 'firep' expected" << std::endl;
      std::exit(1);
    }
    std::string line;
    // Read over label info
    std::getline(instr,line);
    //std::cout << line << std::endl;
    std::getline(instr,line);
    //std::cout << line << std::endl;

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
      //std::cout << "Read " << indices.size() << " indices" << std::endl;

      pre_matrix1.push_back(pre_column(i,grade,indices));
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
      pre_matrix2.push_back(pre_column(i,grade,indices));
    }
    Sort_pre_column<pre_column> sort_pre_column;
    std::sort(pre_matrix1.begin(),pre_matrix1.end(),sort_pre_column);
    std::sort(pre_matrix2.begin(),pre_matrix2.end(),sort_pre_column);

    std::vector<index> re_index;
    re_index.resize(pre_matrix2.size());
    
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

    {
      int n = pre_matrix1.size();
      matrix1.set_num_cols(n);
      
      for(int i=0;i<n;i++) {
	pre_column& pcol = pre_matrix1[i];
        matrix1.grades.push_back(pcol.grade);
	matrix1.set_col(i,pcol.boundary);
      }
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
      matrix2.assign_slave_matrix();
      matrix2.num_rows=r;
      matrix2.assign_pivots();
    }
    assign_grade_indices(matrix1,matrix2);

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
	index start_xy_of_grade = M.start_index_of_pair[std::make_pair(x,y)];
	//std::cout << "Min gens for " << x << " " << y << " traverses through index range " << start_xy << " (" << start_xy_of_grade << ") " << end_xy << std::endl;
	assert(start_xy<=start_xy_of_grade);
	assert(start_xy_of_grade<=end_xy);
	for(index i = start_xy;i<end_xy;i++) {
	  M.reduce_column(i);
	  if(!M.is_empty(i) && i>=start_xy_of_grade) {
	    std::vector<index> col;
	    M.get_col(i,col);
	    //std::cout << "NEW MIN GENERATOR Count" << count << " index " << i << " grade " << x << " " << y << std::endl; 
	    result.set_num_cols(count+1);
	    result.set_col(count++,col);
	    result.grades.push_back(M.grades[i]);
	    result.num_rows=M.num_rows;
	  }
	}
      }
    }
  }

  template<typename GradedMatrix>
    void ker_basis(GradedMatrix& M, GradedMatrix& result) {
    
    index count=0;

    std::set<index> indices_in_kernel;

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
	//std::cout << "Ker basis for " << x << " " << y << " traverses through index range " << start_xy << " " << end_xy << std::endl;
	index start_xy_of_grade = M.start_index_of_pair[std::make_pair(x,y)];
	assert(start_xy<=start_xy_of_grade);
	assert(start_xy_of_grade<=end_xy);
	for(index i = start_xy;i<end_xy;i++) {
	  M.reduce_column(i,true);
	  if(M.is_empty(i) && indices_in_kernel.count(i)==0) {
	    //std::cout << "NEW KERNEL ELEMENT " << i << " Count: " << count << " Grade " << x << " " << y << std::endl;
	    std::vector<index> col;
	    M.slave.get_col(i,col);
	    result.set_num_cols(count+1);
	    result.set_col(count++,col);
	    result.grades.push_back(Grade(x,y,M.x_vals[x],M.y_vals[y]));
	    indices_in_kernel.insert(i);
	  }
	}
      }
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
    
    /*
    std::cout << "Pivots:" << std::endl;
    for(index i=0;i<ker.num_rows;i++) {
      std::cout << i << " -> " << ker.pivots[i] << std::endl;
    }
    */


    for(index i=0;i<cols.get_num_cols();i++) {
      //std::cout << "index" << i << std::endl;
      std::vector<index> col;
      cols.get_col(i,col);
      ker.set_col(ker_cols+i,col);
      //std::cout << "reduce" << std::endl;
      ker.reduce_column(ker_cols+i,true);
      //std::cout << "done" << std::endl;
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
    return std::binary_search(col.begin(),col.end(),m);
  }

  template<typename GradedMatrix_>
    void minimize(GradedMatrix_& M,GradedMatrix_& result) {
    
    GradedMatrix<phat::vector_vector> VVM;
    convert_to_vec_vec(M,VVM);

    typedef std::vector<index> Column;

    std::set<index> rows_to_delete;

    std::vector<Column> cols_to_keep;
    std::vector<Grade> col_grades;

    for(index i=0;i<VVM.get_num_cols();i++) {
      if(VVM.is_empty(i)) {
	continue;
      }
      index col_grade_x = VVM.grades[i].first_index;
      index col_grade_y = VVM.grades[i].second_index;
      index p = M.get_max_index(i);
      index row_grade_x = VVM.row_grades[p].first_index;
      index row_grade_y = VVM.row_grades[p].second_index;
      if(col_grade_x==row_grade_x && col_grade_y==row_grade_y) {
	//std::cout << "Found removable pair " << p << " " << i << std::endl;
	rows_to_delete.insert(p);
	for(index j=i+1;j<VVM.get_num_cols();j++) {
	  if(contains(VVM,j,p)) {
	    VVM.add_to(i,j);
	  }
	}
      } else {
	Column col;
	VVM.get_col(i,col);
	cols_to_keep.push_back(col);
	col_grades.push_back(VVM.grades[i]);
      }
      
    }
    
    index nr = VVM.num_rows;
    //std::cout << "Number of rows of VVM: " << nr << std::endl;
    //std::cout << "Have removed " << rows_to_delete.size() << std::endl;
    //std::cout << "Kept " << cols_to_keep.size() << " columns" << std::endl;
    index count=0;
    std::map<index,index> index_map;
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
      Column& col = cols_to_keep[i];
      for(int j=0;j<col.size();j++) {
	//std::cout << "Lookup " << col[j] << ", count=" << index_map.count(col[j]) << std::endl;
	assert(index_map.count(col[j]));
	col[j]=index_map[col[j]];
      }
    }
    // Now build up the result
    result.set_num_cols(cols_to_keep.size());
    result.num_rows=index_map.size();
    result.grades = col_grades;
    result.row_grades=res_row_grades;
    for(int i=0;i<cols_to_keep.size();i++) {
      result.set_col(i,cols_to_keep[i]);
    }
  }
  
  
}//of namespace phat
