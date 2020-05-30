#!/usr/local/bin/python3

import yaml
import os
import argparse

def run():
  
  with open(args.config, 'r') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
  with open(args.template, 'r') as file:
    cpp = file.read()

  opt_vars = config['optimization_variables']
  derived_vars = config['derived_variables']
  all_vars = list(opt_vars.keys()) + list(derived_vars.keys())
    
  # Class name
  name = os.path.basename(args.config).replace('.yaml', '')
  cpp = cpp.replace('OBJECTIVE_TEMPLATE', name)

  # Variable declarations
  declarations = ''
  for var in all_vars:
    declarations += '  double {};\n'.format(var)
  cpp = cpp.replace('  // DECLARE_VARS', declarations[:-1])

  # Variable assignments
  assignments = ''
  for idx, var in enumerate(opt_vars.keys()):
    assignments += '    {} = vars[{}];\n'.format(var, idx)
  for idx, var in enumerate(derived_vars.keys()):
    assignments += '    {} = {};\n'.format(var, derived_vars[var])
  cpp = cpp.replace('    // SET_VARS', assignments[:-1])  

  # Constraints
  constraints = ''
  for constraint in config['constraints']:
    lhs = constraint.split(' == ')[0]
    rhs = constraint.split(' == ')[1]
    constraints += '    equalityPenalty({}, {});\n'.format(lhs, rhs)
  cpp = cpp.replace('    // SET_CONSTRAINTS', constraints[:-1])

  # Status
  status = ''
  for var in all_vars:
    status += '    oss << prefix << "{} = " << {} << endl;\n'.format(var, var)
  cpp = cpp.replace('    // PRINT_STATUS', status[:-1])

  # Lower bounds
  lower = '    VectorXd lower({});\n'.format(len(opt_vars))
  for idx, var in enumerate(opt_vars.keys()):
    lower += '    lower[{}] = {};\n'.format(idx, opt_vars[var][0])
  cpp = cpp.replace('    // LOWER_BOUNDS', lower[:-1])    

  # Upper bounds
  upper = '    VectorXd upper({});\n'.format(len(opt_vars))
  for idx, var in enumerate(opt_vars.keys()):
    upper += '    upper[{}] = {};\n'.format(idx, opt_vars[var][1])
  cpp = cpp.replace('    // UPPER_BOUNDS', upper[:-1])    
  
  print(cpp)
  

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("config", help="yaml file specifying the optimization problem")
  parser.add_argument("-t", "--template", type=str, default='bds_template.cpp', help='')
  args = parser.parse_args()

  run()

