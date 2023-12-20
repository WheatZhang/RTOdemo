import pyomo
import os
import re
from pyomo.environ import value, Set, RangeSet, Var
import numpy

def get_indexList_name(model, variable_name):
    '''
    Example:
    model.A=Set(initialize=[1,2,3])
    model.B=Set(initialize=[1,2,3])
    model.time = ContinuousSet(initialize=range(0,10))
    model.x=Var(model.time,model.B)
    model.y=Var(model.A)
    print(get_indexList_name(model, 'x'))
    print(get_indexList_name(model, 'y'))
    >>
    ['time', 'B']
    ['A']
    :param model: A pyomo model.
    :param variable_name: String
    :return: A list
    '''
    index = getattr(model, variable_name)._index
    if isinstance(index, pyomo.core.base.set.SetProduct_FiniteSet):
        result = []
        for s in index.subsets():
            result.append(s._name)
        return result
    else:
        return [index._name]


def write_one_var_to_tpl(f, model, attr):
    if isinstance(getattr(model, attr), pyomo.core.base.var.ScalarVar):
        f.write(attr + '\t')
        f.write(str(getattr(model, attr)._value))
        f.write('\n\n')
    elif isinstance(getattr(model, attr), pyomo.core.base.var.IndexedVar):
        index_list = get_indexList_name(model, attr)
        index_line = ""
        first_flag = True
        for i in index_list:
            if first_flag:
                first_flag = False
            else:
                index_line += ","
            index_line += i
        f.write(attr + '\n' + index_line + '\n')
        for key, value in getattr(model, attr)._data.items():
            f.write(str(key) + '\t')
            f.write(str(value._value) + '\n')
        f.write('\n')
    elif isinstance(getattr(model, attr), pyomo.dae.DerivativeVar):
        index_list = get_indexList_name(model, attr)
        index_line = ""
        first_flag = True
        for i in index_list:
            if first_flag:
                first_flag = False
            else:
                index_line += ","
            index_line += i
        f.write(attr + '\n' + index_line + '\n')
        for key, value in getattr(model, attr)._data.items():
            f.write(str(key) + '\t')
            f.write(str(value._value) + '\n')
        f.write('\n')

def to_template(model, file_name):
    all_attr = dir(model)
    f = open(file_name, "w")
    for attr in all_attr:
        write_one_var_to_tpl(f, model, attr)
    f.close()


def load_init_from_tpl_free_collocation(model, file_name,ignore_init_mismatch=False):
    for name, value, index in load_general_txt_init(file_name):
        set_initials_free_collocation(model, name, value, ignore_init_mismatch)


def load_init_from_template(model, file_name,ignore_init_mismatch=False):
    # print(file_name)
    for name, value, index in load_general_txt_init(file_name):
        set_initials(model, name, value, ignore_init_mismatch)


def load_general_txt_init(file_name):
    status = "Ready"
    if not os.path.exists(file_name):
        raise Exception("File not exist: %s."%file_name)
    file = open(file_name)
    for line in file:
        line_split = line.split('\t')
        # print(line_split)
        for i in range(len(line_split)):
            line_split[i] = line_split[i].strip()
        if len(line_split) >= 1 and line_split[0].startswith('#'):
            continue
        if status == "Ready":
            if len(line_split) == 0:
                continue
            elif len(line_split) == 1:
                if line_split[0] == "":
                    continue
                else:
                    status = "WaitIndexLine"
                    var_name = line_split[0]
                    dict = {}
                    continue
            elif len(line_split) == 2:
                if line_split[0] == "":
                    continue
                elif line_split[1] == "":
                    status = "WaitIndexLine"
                    var_name = line_split[0]
                    dict = {}
                    continue
                else:
                    yield line_split[0], to_float_incase_None(line_split[1]),None
                    continue
            else:
                raise Exception("Wrong template format.")
        elif status == 'WaitIndexLine':
            if len(line_split) == 0:
                continue
            elif len(line_split) == 1:
                if line_split[0] == "":
                    continue
                else:
                    status = "ReadIndexedData"
                    index = line_split[0]
                    continue
            elif len(line_split) == 2:
                if line_split[0] == "":
                    continue
                elif line_split[1] == "":
                    status = "ReadIndexedData"
                    index = line_split[0]
                    continue
                else:
                    raise Exception("Wrong template format. Missing index line for the IndexedVar.In line: "+line)
            else:
                raise Exception("Wrong template format.")
        elif status == 'ReadIndexedData':
            if len(line_split) == 0:
                yield var_name, dict, index
                status = "Ready"
                continue
            elif len(line_split) == 1:
                if line_split[0] == "":
                    yield var_name, dict, index
                    status = "Ready"
                    continue
                else:
                    yield var_name, dict, index
                    status = "WaitIndexLine"
                    var_name = line_split[0]
                    continue
            elif len(line_split) == 2:
                if line_split[0] == "":
                    yield var_name, dict, index
                    status = "Ready"
                    continue
                else:
                    # if line_split[1] == "None":
                    #     continue
                    str = line_split[0]
                    if str[0] == "\"" or str[0] == "\'":
                        str = str[1:-1]
                    if str[0] == '(':
                        dict[eval(str)] = to_float_incase_None(line_split[1])
                    else:
                        if str.isdigit():
                            dict[int(str)] = to_float_incase_None(line_split[1])
                        elif re.match(r'[a-zA-Z]',str) == None:
                            if type(eval(str)) == float:
                                dict[float(str)] = to_float_incase_None(line_split[1])
                            else:
                                raise Exception("Wrong index format.")
                        else:
                            try:
                                dict[str] = to_float_incase_None(line_split[1])
                            except ValueError:
                                print(line_split)
                                raise ValueError("could not convert string to float")

            else:
                raise Exception("Wrong template format.")
    if status == 'ReadIndexedData':
        yield var_name, dict, index
    file.close()

def to_float_incase_None(value):
    # value should be a string
    if value == "None":
        return None
    else:
        return float(value)


def set_initials(model, var_name, initial_value,ignore_init_mismatch=False):
    all_attr = dir(model)
    if var_name in all_attr:
        if isinstance(initial_value, dict):
            if (not isinstance(getattr(model, var_name),pyomo.core.base.var.IndexedVar))\
                and (not isinstance(getattr(model, var_name),pyomo.dae.DerivativeVar)):
                if not ignore_init_mismatch:
                    raise Exception("%s is not a IndexedVar or a DerivativeVar" % var_name)
                else:
                    return 0
            for key, value in initial_value.items():
                getattr(model, var_name)._data[key]._value = value
            getattr(model, var_name)._value_init_value = initial_value
        else:
            if not isinstance(getattr(model, var_name),pyomo.core.base.var.ScalarVar):
                if not ignore_init_mismatch:
                    raise Exception("%s is not a ScalarVar" % var_name)
                else:
                    return 0
            getattr(model, var_name)._value = initial_value
            getattr(model, var_name)._value_init_value = initial_value
    else:
        if not ignore_init_mismatch:
            raise Exception("There is no variable called:%s"%var_name)


def find_by_nearest_key(initial_value, key):
    nearest_value=None
    distance=10000
    for k, v in initial_value.items():
        if isinstance(k, tuple):
            flag = True
            for i in range(len(k)-1):
                if k[i+1] != key[i+1]:
                    flag = False
                    break
            if not flag:
                continue
            this_distance = numpy.abs(k[0] - key[0])
        else:
            this_distance = numpy.abs(k - key)
        if this_distance < distance:
            nearest_value = v
            distance = this_distance
    return nearest_value


def set_initials_free_collocation(model, var_name, initial_value,ignore_init_mismatch=False):
    all_attr = dir(model)
    if var_name in all_attr:
        if isinstance(initial_value, dict):
            if (not isinstance(getattr(model, var_name),pyomo.core.base.var.IndexedVar))\
                and (not isinstance(getattr(model, var_name),pyomo.dae.DerivativeVar)):
                if not ignore_init_mismatch:
                    raise Exception("%s is not a IndexedVar or a DerivativeVar" % var_name)
                else:
                    return 0
            initial_value_new_collocation={}
            for key, value in getattr(model, var_name)._data.items():
                initial_value_new_collocation[key]=find_by_nearest_key(initial_value, key)
            for key, value in initial_value_new_collocation.items():
                getattr(model, var_name)._data[key]._value = value
            getattr(model, var_name)._value_init_value = initial_value_new_collocation
        else:
            if not isinstance(getattr(model, var_name),pyomo.core.base.var.ScalarVar):
                if not ignore_init_mismatch:
                    raise Exception("%s is not a ScalarVar" % var_name)
                else:
                    return 0
            getattr(model, var_name)._value = initial_value
            getattr(model, var_name)._value_init_value = initial_value
    else:
        if not ignore_init_mismatch:
            raise Exception("There is no variable called:%s"%var_name)


def unfixed_variables_generator(block):
    """
    Generator which returns all unfixed Var components in a model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all unfixed Var components block
    """
    for v in block.component_data_objects(
            ctype=Var, active=True, descend_into=True):
        if not v.fixed:
            yield v

def back_up_var_value(block,filename):
    with open(filename, 'w') as fp:
        for v in unfixed_variables_generator(block):
            fp.write(str(v) + '\t')
            fp.write(str(value(v)))
            fp.write('\n')

def load_backup_var_value(block, filename):
    with open(filename, 'r') as fp:
        line = fp.readline()
        while line:
            line_split = line.split(sep='\t')
            name_before_revised = line_split[0]

            def replace_rule(matchobj):
                return "\'"+matchobj.group(1)+"\'"
            name_after_revised = re.sub(r'(?<=[\,\[])\s*?([a-zA-Z_][a-zA-Z0-9]*?)\s*?(?=[\,\]])', replace_rule,
                                            name_before_revised)

            eval("block."+name_after_revised)._value = float(line_split[1])
            line = fp.readline()