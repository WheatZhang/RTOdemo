import os
from do_test_lib import generate_noise_file, do_test_MA, \
    do_test_PE, do_test_MAy,do_test_GPE,do_test_ISOPE

#------------------------------------
noise_filename="noise/noise0.txt"
solver_executable=r"F:\Research\RTOdemo\external\bin\ipopt.exe"
result_filename_folder="data/0/"

#------------------------------------
profit_noise_level=0.6
composition_noise_level=1e-4
starting_point={
        "Fb": 4,
        "Tr": 75,
    }

#------------------------------------
perturbation_stepsize={
        "Fb": 0.2,
        "Tr": 2,
    }
filtering_factor=0.5
ka_relative_uncertainty=0.01
kb_relative_uncertainty=0.1
factor_n=10

#------------------------------------
print_iter_data=False
max_iter=20

#------------------------------------
if not os.path.exists(noise_filename):
    generate_noise_file(profit_noise_level, composition_noise_level, noise_filename)
if profit_noise_level <= 0:
    profit_noise_level=0.6
if composition_noise_level <= 0:
    composition_noise_level=1e-4

#------------------------------------
print("\nTesting MA")
result_filename_header=result_filename_folder+"MA_"
do_test_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header)

print("\nTesting MAy")
result_filename_header=result_filename_folder+"MAy_"
do_test_MAy(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header)

print("\nTesting PE")
result_filename_header=result_filename_folder+"PE_"
do_test_PE(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header, composition_noise_level, \
           ka_relative_uncertainty, kb_relative_uncertainty)

print("\nTesting GPE")
result_filename_header=result_filename_folder+"GPE_"
do_test_GPE(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header,profit_noise_level,composition_noise_level,\
                ka_relative_uncertainty,kb_relative_uncertainty,factor_n)

print("\nTesting ISOPE")
result_filename_header=result_filename_folder+"ISOPE_"
do_test_ISOPE(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header, composition_noise_level)