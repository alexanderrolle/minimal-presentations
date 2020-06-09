#include<phat/boundary_matrix.h>
#include<algorithm>

namespace phat {

  typedef std::pair<index,index> index_pair;

  typedef std::priority_queue<index,std::vector<index>,std::greater<index>> PQ;

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

#if SMART_REDUCTION
    std::vector<PQ> pq_row;
#endif

    boundary_matrix<Representation> slave;

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

    void reduce_column(index i, bool use_slave=false, bool notify_pq=false) {
      
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
#if SMART_REDUCTION
      if(notify_pq && p!=-1 && pivots[p]>i) {

	index j = pivots[p];
	index gr_y_index = this->grades[j].second_index;
	//std::cout << "SCHEDULING COLUMN FOR LATER " << i << " " << j << " " << this->grades[j].first_index<< " " << gr_y_index << std::endl;
	this->pq_row[gr_y_index].push(j);
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


  // Assumes that columns are already added to matrix in lex order
  template<typename GradedMatrix>
    void assign_grade_indices(GradedMatrix& M1, GradedMatrix& M2, bool indices_already_assigned=false) {
    
    int n1=M1.get_num_cols();
    int n2=M2.get_num_cols();
      
    if(!indices_already_assigned) {

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
      //std::cout << "Read " << indices.size() << " indices" << std::endl;

      pre_matrix1.push_back(pre_column(i,grade,indices));
    }
    for(int i=0;i<s;i++) {
      std::getline(instr,line);
      std::stringstream sstr(line);
      double x,y;
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
    for(index i=0;i<matrix1.num_rows;i++) {
      matrix1.row_grades.push_back(matrix2.grades[i]);
    }
#if SMART_REDUCTION
    matrix1.pq_row.resize(matrix1.num_grades_y);
    matrix2.pq_row.resize(matrix2.num_grades_y);
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
    for(index i=0;i<M1.get_num_cols();i++) {
      while(M1.is_local(i)) {
	index p = M1.get_max_index(i);
	index j = local_pivots[p];
	if(j!=-1) {
	  M1.add_to(j,i);
	} else {
	  local_pivots[p]=i;
	  break;
	}
      }
    }
    std::cout << "Sparsification" << std::endl;
    for(index i=0;i<M1.get_num_cols();i++) {
      if(M1.is_empty(i) || M1.is_local(i)) {
	continue;
      }
      std::vector<index> col;
      //std::cout << "i=" << i << std::endl;
      while(!M1.is_empty(i)) {
	index p = M1.get_max_index(i);

	index j = local_pivots[p];
	//std::cout << "Pivot=" << p << " " << j << std::endl;
	assert(j==-1 || M1.get_max_index(j)==p);
	if(j!=-1) {
	  M1.add_to(j,i);
	} else {
	  col.push_back(p);
	  M1.remove_max(i);
	}
      }
      std::sort(col.begin(),col.end());
      M1.set_col(i,col);
    }

    std::cout << "Build up smaller matrices" << std::endl;
    
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
    
    result1.set_num_cols(col_count);
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
    result2.set_num_cols(row_count);
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
#if 0
    assign_grade_indices(result1,result2);
#else
    result1.num_grades_x=M1.num_grades_x;
    result1.num_grades_y=M1.num_grades_y;
    std::copy(M1.x_vals.begin(),M1.x_vals.end(),std::back_inserter(result1.x_vals));
    std::copy(M1.y_vals.begin(),M1.y_vals.end(),std::back_inserter(result1.y_vals));
    

    result2.num_grades_x=M2.num_grades_x;
    result2.num_grades_y=M2.num_grades_y;
    std::copy(M2.x_vals.begin(),M2.x_vals.end(),std::back_inserter(result2.x_vals));
    std::copy(M2.y_vals.begin(),M2.y_vals.end(),std::back_inserter(result2.y_vals));

    assign_grade_indices(result1,result2,true);
#endif

#if SMART_REDUCTION
    result1.pq_row.resize(result1.num_grades_y);
    result2.pq_row.resize(result2.num_grades_y);
#endif    

    //result1.print(true,true);
    //result2.print(false,true);

  }



#if SMART_REDUCTION
  template<typename GradedMatrix>
    void min_gens(GradedMatrix& M, GradedMatrix& result) {
    
    index count=0;

    for(index x = 0; x < M.num_grades_x;x++) {
      for(index y = 0; y < M.num_grades_y;y++) {
	
	//std::cout << "Min gens " << x << " " << y << std::endl;
	
	PQ& pq = M.pq_row[y];

	index end_xy;
	if(x<M.num_grades_x-1) {
	  end_xy = M.start_index_of_pair[std::make_pair(x+1,y)];
	} else {
	  end_xy = M.start_index_of_pair[std::make_pair(0,y+1)];
	}
	if(x==M.num_grades_x-1 && y==M.num_grades_y-1) {
	  end_xy=M.get_num_cols();
	}
	index start_xy = M.start_index_of_pair[std::make_pair(x,y)];
	//std::cout << "Min gens for " << x << " " << y << " traverses through index range " << start_xy << " " << end_xy << std::endl;
	assert(start_xy<=end_xy);
	//std::cout << "Before adding, pq of row has size " << pq.size() << std::endl;
	for(index i = start_xy;i<end_xy;i++) {
	  pq.push(i);
	}
	//std::cout << "After adding, pq of row has size " << pq.size() << std::endl;
	while(!pq.empty()) {
	  index i = pq.top();
	  //std::cout << "Index " << i << std::endl;
	  // Remove duplicates
	  while(!pq.empty() && i==pq.top()) {
	    pq.pop();
	  }
	  assert(M.grades[i].first_index<=x);
	  assert(M.grades[i].second_index==y);
	  M.reduce_column(i,false,true);
	  if(!M.is_empty(i) && i>=start_xy&& i<end_xy) {
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
#else

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
	  if(!M.is_empty(i) && i>=start_xy_of_grade&& i<end_xy) {
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

#endif

#if SMART_REDUCTION

  template<typename GradedMatrix>
    void ker_basis(GradedMatrix& M, GradedMatrix& result) {
    
    index count=0;

    std::set<index> indices_in_kernel;

    for(index x = 0; x < M.num_grades_x;x++) {
      for(index y = 0; y < M.num_grades_y;y++) {

	PQ& pq = M.pq_row[y];
	
	index start_xy = M.start_index_of_pair[std::make_pair(x,y)];
	index end_xy;
	if(x<M.num_grades_x-1) {
	  end_xy = M.start_index_of_pair[std::make_pair(x+1,y)];
	} else {
	  end_xy = M.start_index_of_pair[std::make_pair(0,y+1)];
	}
	if(x==M.num_grades_x-1 && y==M.num_grades_y-1) {
	  end_xy=M.get_num_cols();
	}
	assert(start_xy<=end_xy);

	for(index i = start_xy;i<end_xy;i++) {
	  pq.push(i);
	}
	//std::cout << "After adding, pq of row has size " << pq.size() << std::endl;
	test_timer2.resume();
	while(!pq.empty()) {
	  index i = pq.top();
	  // Remove duplicates
	  while(!pq.empty() && i==pq.top()) {
	    pq.pop();
	  }
	  assert(M.grades[i].first_index<=x);
	  assert(M.grades[i].second_index==y);
	  M.reduce_column(i,true,true);
	  if(M.is_empty(i) && indices_in_kernel.count(i)==0) {
	    //std::cout << "NEW KERNEL ELEMENT " << i << " Count: " << count << " Grade " << x << " " << y << std::endl;
	    {
	      std::vector<index> col;
	      M.slave.get_col(i,col);
	      result.set_num_cols(count+1);
	      result.set_col(count++,col);
	      result.grades.push_back(Grade(x,y,M.x_vals[x],M.y_vals[y]));
	      indices_in_kernel.insert(i);
	    }
	  }
	}
	test_timer2.stop();	
      }
    }
    test_timer1.start();
    result.num_rows=M.get_num_cols();
    result.slave.set_num_cols(result.get_num_cols());
    result.assign_pivots();
    for(index i=0;i<result.get_num_cols();i++) {
      result.pivots[result.get_max_index(i)]=i;
      std::vector<index> slave_col;
      slave_col.push_back(i);
      result.slave.set_col(i,slave_col);
    }
    test_timer1.stop();

  }


#else

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


#endif

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


#if LAZY_MINIMIZATION
    GradedMatrix_& VVM=M;
#else    
    GradedMatrix<phat::vector_vector> VVM;
    convert_to_vec_vec(M,VVM);
#endif

    typedef std::vector<index> Column;


    VVM.assign_pivots();

    std::set<index> rows_to_delete;

    std::vector<index> cols_to_keep;
    std::vector<Grade> col_grades;

    for(index i=0;i<VVM.get_num_cols();i++) {
      //std::cout << "i=" << i << std::endl;
      while(! VVM.is_empty(i)) {
	
	index col_grade_x = VVM.grades[i].first_index;
	index col_grade_y = VVM.grades[i].second_index;
	index p = VVM.get_max_index(i);
	index row_grade_x = VVM.row_grades[p].first_index;
	index row_grade_y = VVM.row_grades[p].second_index;
	//std::cout << "Pivot is " <<  p << std::endl;
	if(col_grade_x!=row_grade_x || col_grade_y!=row_grade_y) {
	  cols_to_keep.push_back(i);
	  col_grades.push_back(VVM.grades[i]);
	  break;
	}
	if(VVM.pivots[p]==-1) {
	  //std::cout << "Found removable pair " << p << " " << i << std::endl;
	  rows_to_delete.insert(p);
#if !LAZY_MINIMIZATION
	  for(index j=i+1;j<VVM.get_num_cols();j++) {
	    if(contains(VVM,j,p)) {
	      VVM.add_to(i,j);
	    }
	  }
#endif
	  VVM.pivots[p]=i;
	  break;
	} else {
	  //std::cout << "add a column" << VVM.pivots[p] << std::endl;
	  VVM.add_to(VVM.pivots[p],i);
	}
      }
    }

    std::vector<Column> new_cols;
    new_cols.resize(cols_to_keep.size());

#if LAZY_MINIMIZATION
    //std::cout << "HERE I AM " << std::endl;
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
	}
      }
      std::reverse(col.begin(),col.end());
    }
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
      Column& col = new_cols[i];
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
      result.set_col(i,new_cols[i]);
    }
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
    
    std::map<index,index> re_index_rows;

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

    /*
    out << "x-grades" << std::endl;
    for(index i=0;i<M.num_grades_x;i++) {
      out << M.x_vals[i] << std::endl;
    }
    out << std::endl;
    out << "y-grades" << std::endl;
    for(int i=0;i<M.num_grades_x;i++) {
      out << M.y_vals[i] << std::endl;
    }
    out << std::endl;
    */
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
