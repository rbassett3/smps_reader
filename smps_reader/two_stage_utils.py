import numpy as np
import scipy.sparse
import mps_reader

def extract_matrix_data(parsed_file_dicts):
    '''construct_vecs_and_mats takes a dictionary from
    smps_reader and returns the matrix data defining 
    a two stage problem with discrete scenarios'''
    #extract the dictionaries for each file for further use
    core = parsed_file_dicts['core']
    time = parsed_file_dicts['time']
    stoch = parsed_file_dicts['stoch']
    assert len(time['periods']) == 2, "This problem is not 2 stage"
    assert stoch['scenarios_flag'], "This problem does not list scenarios"

    #build the matrices for the .core file from the extracted dictionary 
    matrix_data = mps_reader.extract_matrix_data(core)

    #construct partition of variables into 1st and 2nd stages
    if time['format'] == 'implicit':
        stage1vars = []
        stage2vars = []
        othervars = []
        stage = 0
        stages_firstvars = \
            [time['periods'][period]['col_start'] for period in time['periods'].keys()]
        stage1_firstvar, stage2_firstvar = stages_firstvars[0], stages_firstvars[1]
        col2ind = {}
        for i, var in enumerate(matrix_data['col_labels']):
            col2ind[var] = i
            if var == stage1_firstvar:
                stage=1
            if var == stage2_firstvar:
                stage=2
            if stage==1:
                stage1vars.append(var)
            if stage==2:
                stage2vars.append(var)
            if stage==0:
                othervars.append(var)
        if len(othervars) != 0:
            print("Some variables aren't assigned to a stage! They are:")
            print(othervars)
            assert False
      
        #this is pretty good but I need to find a way to index into each array  
        bool_stage1vars = np.zeros(len(matrix_data['col_labels']), dtype=np.bool_)
        bool_stage1vars[[col2ind[var] for var in stage1vars]] = True
        ZeroWBlock = matrix_data['A'][:, ~bool_stage1vars]
        ATBlock = matrix_data['A'][:,bool_stage1vars]
        A_block_inds = np.array([ZeroWBlock[i,:].nnz==0 for i in range(ZeroWBlock.shape[0])])
        A = ATBlock[A_block_inds,:]
        T = ATBlock[~A_block_inds,:]
        W = ZeroWBlock[~A_block_inds,:]
        b = matrix_data['b'][A_block_inds]
        r = matrix_data['b'][~A_block_inds]
        l1 = matrix_data['l'][bool_stage1vars]
        u1 = matrix_data['u'][bool_stage1vars]
        l2 = matrix_data['l'][~bool_stage1vars]
        u2 = matrix_data['u'][~bool_stage1vars]
        c = matrix_data['c'][bool_stage1vars]
        q = matrix_data['c'][~bool_stage1vars]
        ineq_b = matrix_data['ineq_b'][A_block_inds]
        ineq_r = matrix_data['ineq_b'][~A_block_inds]

        var2ATind = dict(zip(stage1vars, range(len(stage1vars)))) #index in A or T
        var2Wind = dict(zip(stage2vars, range(len(stage2vars)))) #index in W
        row2Aind = dict(zip(
            [lab for i, lab in enumerate(matrix_data['row_labels']) 
                if A_block_inds[i]],
            range(A_block_inds.sum())))
        row2WTind = dict(zip(
            [lab for i, lab in enumerate(matrix_data['row_labels']) 
                if ~A_block_inds[i]],
            range((~A_block_inds).sum())))

        #guaranteed that there is only one objective row from
        #mps_reader
        obj_row = [row for row in core['rows'].keys() \
                    if core['rows'][row]=='N'][0]

        prob_data = {'A':A, 'b':b, 'c':c, 'l1':l1, 'u1':u1, 'l2':l2, 'u2':u2,\
            'T_root':T, 'W_root':W, 'r_root':r, 'q_root':q, 'ineq_b':ineq_b,\
            'ineq_r':ineq_r, 'scenarios':{}}
     
        for scen in stoch['scenarios'].keys():
            this_scen = {}
            this_scen['prob'] = stoch['scenarios'][scen]['prob']
            parent = stoch['scenarios'][scen]['parent']
            #"'Root'" is required, but common typo which
            #we gracefully deal with
            if parent == 'ROOT' or "'ROOT'":
                this_scen['T'] = prob_data['T_root'].copy()
                this_scen['W'] = prob_data['W_root'].copy()
                this_scen['q'] = prob_data['q_root'].copy()
                this_scen['r'] = prob_data['r_root'].copy()
                #todo: support bounds by having l2 and u2
                #depend on scenario. File type supports this
            else:
                par_scen = prob_data[scenarios][parent]
                this_scen['T'] = par_scen['T'].copy()
                this_scen['W'] = par_scen['W'].copy()
                this_scen['q'] = par_scen['q'].copy()
                this_scen['r'] = par_scen['r'].copy()
                #todo: support bounds by having l2 and u2
                #depend on scenario. File type supports this

            #next we loop through the data for this scenario
            #there are 3 types of data here
            #stoch['scenarios'][scen] is a tuple containing
            #(type, bound, column, value) in the case of bound updates
            #('', range, row, value) in the case of range updates
            #('', col, row, value) in the case of coefficient updates
            #('', rhs, row, value) in the case of rhs updates
            #each of these correspond to the formatting of the mps file
            for data in stoch['scenarios'][scen]['data']:
                isbound = (data[0] != '')
                if isbound: #it's a bound update
                    print("It's a bound update!")
                    assert False, "Not supported yet"
                elif data[2] == obj_row:
                    #It's an objective update
                    this_scen['q'][var2Wind[data[1]]] = data[3]
                    #I am confused. The docs from haussman's website makes it clear
                    #that rhs side updates should just look like a rhs data field.
                    #but I keep seeing files that use RHS as the first entry as the data
                    #field. I will work around it below
                elif data[1]=='RHS' or data[1] in core['rhs'].keys():
                    #It's a rhs update
                    this_scen['r'][row2WTind[data[2]]] = data[3]
                elif data[1] in core['ranges']:
                    print("It's a range update!")
                    assert False, "Not supported yet"
                elif data[1] in var2Wind.keys():
                    #It's a W update
                    this_scen['W'][row2WTind[data[2]], var2Wind[data[1]]]\
                        = data[3]
                elif data[1] in var2ATind.keys():
                    #It's a T update!
                    this_scen['T'][row2WTind[data[2]], var2ATind[data[1]]]\
                        = data[3]
                else:
                    print("data is", data)
                    assert False, "not a recognized update!"
            prob_data['scenarios'][scen] = this_scen
    else: #explicit scenarios
        assert False, "Only implicit scenarios have been implemented"
    return prob_data
