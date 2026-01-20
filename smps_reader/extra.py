from smps_reader import read

def two_stage_matvecs(path_to_smps_file, core_file=None,\
  time_file=None, stoch_file=None, strict=True):
    '''construct the matrices and vectors from a 2 stage
    stochastic program (w/ stoch file in scenario form)
    in .smps format such that the problem is:
    min_{x, y_{s}} c^{T} x + \sum_{s \in S} p_{s} q_{s}^{T} y_{s}
    s.t. A x = b
         T_{s} x + Q_{s} y_{s} = r_{s},  s \in S
    Any integer or binary restrictions are ignored.
    First stage and second stage constraints can also 
    be inequalities (<=). Separate variables, ineq_A and
    T_ineq, specify the rows of the A and T, respectively
    that correspond to inequality constraints.
    '''
    prob_data = read(path_to_smps_file, core_file=core_file,\
      time_file=time_file, stoch_file=stoch_file, strict=strict)
    stages = [stage for stage in prob_data['time'].keys() if\
      stage != 'prob_name'] #all keys are stages except for prob_name
    assert len(stages) == 2, "Problem is not 2 stage"
    assert prob_data['scenarios'],\
      "Problem uncertainty is not represented as scenarios"
    #all keys are scenarios except for prob_name and the flags
    #that indicate if the dict contains scenarios or distributions
    scenarios = [scenario for scenario in prob_data['stoch'].keys()\
      if scenarios not in ['scenario', 'distributions', 'prob_name']]
        #extract variables from parsed_file_dict for notational convenience
    
    #in the smps spec it says that the objective must always
    #belong to the first period. But the very first example I found
    #breaks this; it lists the objective first and then constraint
    #rows. The time file then specifies that the first period
    #STARTS with the first constraint row. I will try to work
    #with this, i.e. I will attempt to make it so that it
    #doesn't matter what period you assign the objective to
    #(https://web.archive.org/web/20050618080243/http://www.mgmt.dal.ca/sba/profs/hgassmann/SMPS2.htm#TimeFile)
    assert 'row_start' in prob_data['time'][stages[0]],\
      """Periods are in explicit format.
      Only implicit supported (for now).)"""
    first_stg_1st_col = prob_data['time'][stages[0]]['col_start']
    first_stg_1st_row = prob_data['time'][stages[0]]['row_start']
    second_stg_1st_col = prob_data['time'][stages[1]]['col_start']
    second_stg_1st_row = prob_data['time'][stages[1]]['row_start']
    first_stg_vars = []
    second_stg_vars = []
    hit_2nd_stg = False
    for var in prob_data['core']['columns'].keys():  
        #just grab first elt then break
        assert var==first_stg_1st_col,\
          "First variable not in first stage. Aborting..."
        break
    for var in prob_data['core']['columns'].keys():
        if var==second_stg_1st_col:
            hit_2nd_stg=True
        if hit_2nd_stg:
            second_stg_vars.append(var)
        else:
            first_stg_vars.append(var)
    x = np.empty(len(first_stg_vars))
    c = np.empty(len(first_stg_vars))
    ys = [np.empty(len(second_stg_vars))]
    qs = [np.empty(len(second_stg_vars))]
     

    #COPY PASTED FROM MPS_READER
    #note that obj_shift and prob_name are not used
    rows = parsed_file_dict['rows']
    columns = parsed_file_dict['columns']
    rhs = parsed_file_dict['rhs']
    ranges = parsed_file_dict['ranges']
    bounds = parsed_file_dict['bounds']
    #now we actually build the vectors and matrices
    col_to_ind = dict(zip(columns.keys(), range(len(columns.keys()))))
    n = len(col_to_ind)
    rows_in_objs = [row for row in rows.keys() if rows[row]=='N']
    assert len(rows_in_objs) == 1, "More than 1 objective specified"
    #number of inequality constraints is number of rows w/ type E
    rows_in_eq_ind = [row for row in rows.keys() if rows[row]=='E']
    #number of inequality constraints is number of rows w/ type L or G + 
    #number of additional L and G constraints implied by RANGE
    #the range constraints have an emoji added to their labels
    #which allows referencing them later as unique rows
    #by adding an emoji we guarantee no conflict with another
    #row name since mps only permits ascii characters
    #surprisingly, range on a E row yields two inequality constraints
    rows_in_ub_ind = [row for row in rows.keys() if rows[row] in ['L','G']] +\
        [row+"\U0001f600" for range_ in ranges.values()\
            for (row, value) in range_ if rows[row] in ['L','G']]+\
        [row+"\U0001f600" for range_ in ranges.values()\
            for (row, value) in range_ if rows[row]=='E']+\
        [row+"\U0001f606" for range_ in ranges.values()\
            for (row, value) in range_ if rows[row]=='E']
    m1 = len(rows_in_eq_ind)
    m2 = len(rows_in_ub_ind)
    row_to_eq_ind = dict(zip(rows_in_eq_ind, range(m1)))
    row_to_ub_ind = dict(zip(rows_in_ub_ind, range(m2)))
    #number of fixed variables. We'll treat these separately
    num_fixed = len([kind for bnds in bounds.values() for (kind, col, val) in bnds if kind=='FX'])
    fixed_inds = np.empty(num_fixed, dtype=np.int64)
    fixed_vals = np.empty(num_fixed, dtype=np.float64)
    fixed_itr = 0 #how many we've seen as we loop through them later

    c = np.zeros(n)
    l = np.zeros(n)
    u = np.inf*np.ones(n)
    b_eq = np.zeros(m1)
    b_ub = np.zeros(m2)
    A_eq = scipy.sparse.dok_matrix((m1, n), dtype=np.float64)
    A_ub = scipy.sparse.dok_matrix((m2, n), dtype=np.float64)

    #loop through column section to build c vector and A matrices
    #loop through column section to build c vector and A matrices
    for column in columns.keys():
        rows_for_this_column = columns[column]
        for (row, value) in rows_for_this_column:
            col_ind = col_to_ind[column]
            if rows[row]=='N': #objective
                c[col_ind] = value
            elif rows[row]=='L': #lower bound. negate b/c we only keep track of A_ub @ x <= b_ub
                A_ub[row_to_ub_ind[row], col_ind] = -1.0*float(value)
            elif rows[row]=='G': #upper bound. don't have to negate
                A_ub[row_to_ub_ind[row], col_ind] = float(value)
            elif rows[row]=='E': #equality constraint
                A_eq[row_to_eq_ind[row], col_ind] = float(value)
            else:
                raise ValueError("Row kind " + rows[row] + " not recognized")
    #loop through rhs to build b vectors
    for this_rhs_name in rhs.keys():
        for (row, value) in rhs[this_rhs_name]:
            if rows[row]=='L': #lower bound. negate b/c we only keep track of A_ub @ x <= b_ub
                b_ub[row_to_ub_ind[row]] = -1.0*float(value)
            elif rows[row]=='G': #upper bound. don't have to negate
                b_ub[row_to_ub_ind[row]] = float(value)
            elif rows[row]=='E': #equality constraint
                b_eq[row_to_eq_ind[row]] = float(value)
            else:
                raise ValueError("Row kind " + rows[row] + " not recognized")
    for range_ in ranges.keys():
        for (row, value) in ranges[range_]:
            #The range section has a different meaning depending on whether
            #the row it references is of kind G, L, or E. What I do here
            #is in page 164 of Advanced Linear Programming by Murtagh
            #note that per the mps rules ROWS must come before RANGE (if RANGE exists)
            if rows[row]=='L': #Range adds an upper bound of b[i] + abs(r[i])
                A_ub[row_to_ub_ind[row+"\U0001f600"],:] = A_ub[row_to_ub_ind[row],:]
                b_ub[row_to_ub_ind[row+"\U0001f600"]] = b_ub[row_to_ub_ind[row]]+abs(float(value))
            elif rows[row]=='G': #Range adds a lower bound of b[i] - abs(r[i])
                A_ub[row_to_ub_ind[row+"\U0001f600"],:] = -A_ub[row_to_ub_ind[row],:]
                b_ub[row_to_ub_ind[row+"\U0001f600"]] = -1.*(b_ub[row_to_ub_ind[row]]-abs(float(value)))
            elif rows[row]=='E': #equality constraint. Adds an upper and a lower bound
                #where value depends on the sign
                sign_val = float(value) >= 0
                if sign_val: #(b, b+|r|) constraint
                    A_ub[row_to_ub_ind[row+"\U0001f600"],:] = -A_ub[row_to_ub_ind[row],:]
                    b_ub[row_to_ub_ind[row+"\U0001f600"]] = -1.*b_ub[row_to_ub_ind[row]]
                    A_ub[row_to_ub_ind[row+"\U0001f606"],:] = A_ub[row_to_ub_ind[row],:]
                    b_ub[row_to_ub_ind[row+"\U0001f606"]] = b_ub[row_to_ub_ind[row]]+abs(float(value))
                else: #(b-|r|, b) constraint
                    A_ub[row_to_ub_ind[row+"\U0001f600"],:] = -A_ub[row_to_ub_ind[row],:]
                    b_ub[row_to_ub_ind[row+"\U0001f600"]] = -1.*(b_ub[row_to_ub_ind[row]]-\
                                                                abs(float(value)))
                    A_ub[row_to_ub_ind[row+"\U0001f606"],:] = A_ub[row_to_ub_ind[row],:]
                    b_ub[row_to_ub_ind[row+"\U0001f606"]] = b_ub[row_to_ub_ind[row]]
            else:
                raise ValueError("Row kind " + rows[row] + " not recognized")
    #loop through bounds to build l and u vectors
    for bnd in bounds.keys():
        for (kind, column, value) in bounds[bnd]:
            if kind == 'UP':
                u[col_to_ind[column]] = value
            elif kind == 'LO':
                l[col_to_ind[column]] = value
            elif kind == 'FX': #why do these variables even exist?
                #they should be added to the b terms instead
                fixed_inds[fixed_itr] = col_to_ind[column]
                fixed_vals[fixed_itr] = float(value)
                fixed_itr += 1
            elif kind == 'MI': #x \in (-\infty, 0)
                u[col_to_ind[column]] = 0.0
                l[col_to_ind[column]] = -np.inf
            elif kind == "PL": #x \in (0, \infty)
                continue #this is the default bound. Nothing to do
            else:
                raise ValueError("Bound kind " + kind + " not recognized")
    #convert the completed matrices from dok to csr
    A_eq = scipy.sparse.csr_matrix(A_eq)
    A_ub = scipy.sparse.csr_matrix(A_ub)
    return {'c':c, 'A_eq':A_eq, 'A_ub':A_ub, 'b_eq':b_eq, 'b_ub':b_ub, 'l':l, 'u':u,\
            'fixed_inds':fixed_inds, 'fixed_vals':fixed_vals}

    



    

