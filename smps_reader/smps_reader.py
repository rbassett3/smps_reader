from mps_reader import parse_mps_file as parse_core_file
from mps_reader import get_fields_
from mps_reader import reset_flags_to_false

def read(path_to_smps_file, core_file=None, time_file=None,
    stoch_file=None, strict=True):
    '''read takes a path to an smps file as input.
    Problems written in smps format have 3 files, a
    core file, and time file, and a stochastics (stochs)
    file. These are assumed to have extensions .cor, .tim,
    and .sto, respectively.

    The smps file provided to read can be any one of these
    files. The others will be found by replacing the file
    extension of the provided path.

    Otherwise, you can provide the files directly via
    keyword arguments'''

    base_path = path_to_smps_file.split('.')[0]
    if core_file is None:
        core_file = base_path + ".cor"
    if time_file is None:
        time_file = base_path + ".tim"
    if stoch_file is None:
        stoch_file = base_path + ".sto"
        
    core_dict = parse_core_file(core_file, strict=strict)
    time_dict = parse_time_file(time_file, strict=strict)
    stoch_dict = parse_stoch_file(stoch_file, strict=strict)
    assert core_dict['prob_name'] == time_dict['prob_name'] \
      and time_dict['prob_name'] == stoch_dict['prob_name'],\
      "Problem name inconsistent across files"
    return {'core':core_dict, 'time':time_dict, 'stoch':stoch_dict}
    
def parse_time_file(path_to_time_file, strict=True):
    #define the strict (or not) reader for data records
    get_fields = lambda x: get_fields_(x, strict=strict)
    periods = {}
    flags = {'in_periods':False, 'explicit':False, 
        'in_rows':False, 'in_columns':False}
    with open(path_to_time_file, 'r') as f:
        for line in f:
            if line[0] != ' ': #this is a section
                reset_flags_to_false(flags)
                sec_name = line.split()[0]
                if sec_name == "TIME":
                    prob_name = line.split()[1]
                elif sec_name == "PERIODS":
                    #assume that rows & columns are ordered
                    #by the period they are in, otherwise
                    #EXPLICIT needs to be specified
                    if len(line.split()) == 1 \
                      or line.split()[1] == 'EXPLICIT':
                        flags['explicit'] = True
                    flags['in_periods'] = True
                elif sec_name == "ROWS":
                    flags['in_rows'] = True
                elif sec_name == "COLUMNS":
                    flags['in_columns'] = True
                elif sec_name == "ENDATA":
                    break #end of file
                else:
                    raise ValueError("SMPS time file has unrecognized section " + sec_name)
            elif flags['in_periods']:
                #ignore 'core' field
                _, col, row, period, _, _ = get_fields(line)
                if flags['explicit']: #col & row are not used
                    periods[period]['rows'] = []
                    periods[period]['col'] = []
                else: 
                    periods[period] = {'row_start':row,\
                      'col_start':col}
            elif flags['in_rows']:
                #format says this only exists if 
                #PERIOD is explicit
                _, row, period, _, _, _ = get_fields(line)
                periods[period]['rows'].append(row)
            elif flags['in_columns']:
                #format says this only exists if 
                #PERIOD is explicit
                col, period, _, _, _, _ = get_fields(line)
                periods[period]['cols'].append(col)
    periods['prob_name'] = prob_name
    return periods

def parse_stoch_file(path_to_stoch_file, strict=True):
    '''Doesn't support the following sections:
    SIMPLE
    ROBUST
    PLINQUAD
    CHANCE
    ICC
    NODE
    DISTRIB
    '''
    get_fields = lambda x: get_fields_(x, strict=strict)
    flags = {'in_scenarios':False, 'in_indep':False,
        'in_blocks':False}
    return_scenarios = False
    return_discrete = False
    return_dict = {}
    with open(path_to_stoch_file, 'r') as f:
        for line in f:
            if line[0] != ' ': #this is a section
                reset_flags_to_false(flags)
                sec_name = line.split()[0]
                #print("In section ", sec_name) 
                if sec_name == "STOCH":
                    prob_name = line.split()[1]
                elif sec_name == "SCENARIOS":
                    flags['in_scenarios'] = True
                    return_scenarios = True
                elif sec_name == "INDEP":
                    rv_type = line.split()[1]
                    assert rv_type == 'DISCRETE',\
                        "Only DISCRETE supported at this time"
                    #assert rv_type in ['DISCRETE', 'UNIFORM'\
                    #    'NORMAL', 'GAMMA', 'BETA', 'LOGNORM']
                    if len(line.split()) == 3:
                        modify_type = line.split()[2]
                        assert modify_type == 'REPLACE',\
                            "Only REPLACE supported at this time."
                        #assert modify_type in ['ADD', 'MULTIPLY',\
                        #    'REPLACE']
                    else: #default is REPLACE
                        modify_type = 'REPLACE'
                    flags['in_indep'] = True
                    return_discrete = True
                elif sec_name == "BLOCKS":
                    rv_type = line.split()[1]
                    assert rv_type == 'DISCRETE',\
                        "Only DISCRETE supported at this time"
                    #assert rv_type in ['DISCRETE', 'MVNORMAL']
                    if len(line.split()) == 3:
                        modify_type = line.split()[2]
                        assert modify_type == 'REPLACE',\
                            "Only REPLACE supported at this time."
                        #assert modify_type in ['ADD', 'MULTIPLY',\
                        #    'REPLACE']
                    else: #default is REPLACE
                        modify_type = 'REPLACE'

                    flags['in_blocks'] = True
                    return_discrete = True
                elif sec_name == "ENDATA":
                    break #end of file
                else:
                    raise ValueError("SMPS time file has unrecognized section " + sec_name)
            elif flags['in_scenarios']:
                field1, field2, field3, field4, field5, field6 = \
                    get_fields(line)
                if field1 == "SC": #this is a new scenario
                    this_scen = field2
                    #pdb.set_trace()
                    return_dict[this_scen] = {}
                    return_dict[this_scen]['parent'] = field3
                    return_dict[this_scen]['prob'] = float(field4)
                    #seems like field 6 is unused?
                    return_dict[this_scen]['period'] = field6
                    return_dict[this_scen]['data'] = []
                else:
                    #more information on current scenario
                    #either
                    #type, bound, column, value
                    #or 
                    #"", column, row, value (for matrices in constraints)
                    #or
                    #"", rhs, column, value (for rhs in constraints)
                    return_dict[this_scen]['data'].append(
                        (field1, field2, field3, float(field4)))

            elif flags['in_indep']:
                field1, field2, field3, field4, field5, field6 = \
                    get_fields(line)
                #note. period is often ''
                col, row, period = field2, field3, field5
                if ((col, row), period) in return_dict.keys():
                    return_dict[(col, row, period)]\
                      ['values'].append(float(field4))
                    return_dict[(col, row, period)]\
                      ['probs'].append(float(field6)) 
                else:
                    return_dict[((col, row), period)] =\
                    {'values':[float(field4),], 'probs':[float(field6),]}
            elif flags['in_blocks']:
                field1, field2, field3, field4, field5, field6 = \
                    get_fields(line)
                if field1 == 'BL': #new block. 
                    this_block = field2
                    period = field3
                    if (this_block, period) in\
                      return_dict.keys():
                        return_dict[(this_block, period)] \
                          = [{'prob':float(field4), 'col/row':[],\
                             'value':[]}]
                    else:
                        return_dict[(this_block, period)]\
                        .append({'prob':float(field4), 'col/row':[],\
                             'value':[]})
                return_dict[(this_block, period)]['col/row'].append((field2, field3))
                return_dict[(this_block, period)]['value'].append(field4)

        return_dict['scenarios'] = return_scenarios
        return_dict['distributions'] = return_discrete
        return_dict['prob_name'] = prob_name
        assert return_scenarios or return_discrete, "Neither distribution nor scenario representation"
        return return_dict #returns a dictionary contaning scenarious or discrete distributions on elements.
        #the scenarios and distributions keys tell these cases apart.

                    
                
                        
            
    
               

