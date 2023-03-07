from example.TR.tr_quadr_eg import original_MA_with_max_tr, algo2_TR_MA, algo1_TR_MA
import os
import matplotlib.pyplot as plt
from example.TR.draw_lib import plot_contour
import pandas
import numpy

global_linewidth = 1
global_point_size= 4
pic_constant = 0.39370
font_factor = numpy.sqrt(1/pic_constant)
font_legend = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 17/font_factor,
         }
font_title = {'family': 'Times New Roman',
              'size': 17/font_factor,
              'weight': 'normal'
              }
font_label = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 17/font_factor,
             }
global_tick_size=17/font_factor
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def do_test_global_convergence():
    u1_0 = 2 #1.1
    u2_0 = -2 #0.1
    r_0 = 1
    sigma = 10
    xi_N = 0.5
    r_max = 2
    filtering_factor=0.1
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/escape/global_conv/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    starting_point = {
        "u1": u1_0,
        "u2": u2_0,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "u1": 0.01,
        "u2": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = r_0
    sigma = sigma
    xi_N = xi_N


    # ------------------------------------
    for model_name in ["norminal", "obj_partial_wrong_curv", "con_wrong_curv2"]:
        print("\nTesting Original_MA_with_Rmax")
        result_filename_header = result_filename_folder + model_name+"_MA_"
        original_MA_with_max_tr(model_name, perturbation_stepsize, starting_point, r_max, filtering_factor, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

        # ------------------------------------
        print("\nTesting CompoStep_TR_MA")
        result_filename_header = result_filename_folder + model_name+"_CompoStep_"
        algo1_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius,r_max, xi_N, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

def do_test_penalty():
    f_u1_0 = 2 #1.1
    f_u2_0 = -2 #0.1
    inf_u1_0 = 1.5
    inf_u2_0 = 1.5
    r_0 = 1
    xi_N = 0.5
    r_max = 2
    filtering_factor=0.1
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/escape/penalty_coeff/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    feasible_starting_point = {
        "u1": f_u1_0,
        "u2": f_u2_0,
    }
    infeasible_starting_point = {
        "u1": inf_u1_0,
        "u2": inf_u2_0,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "u1": 0.01,
        "u2": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = r_0
    xi_N = xi_N
    model_name = "norminal"

    # ------------------------------------
    name=["LP","MP","SP"]
    for no,sigma in enumerate([1,0.1,0.01]):
        print("\nTesting Penalty_MA")
        result_filename_header = result_filename_folder + name[no]+"F_Penalty_"
        algo2_TR_MA(model_name, perturbation_stepsize, feasible_starting_point, sigma, initial_trust_radius, r_max, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

        print("\nTesting CompoStep_TR_MA")
        result_filename_header = result_filename_folder + name[no]+"F_CompoStep_"
        algo1_TR_MA(model_name, perturbation_stepsize, feasible_starting_point, sigma, initial_trust_radius,r_max, xi_N, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)
        # ------------------------------------
        # print("\nTesting Penalty_MA")
        # result_filename_header = result_filename_folder + name[no] + "INF_Penalty_"
        # algo2_TR_MA(model_name, perturbation_stepsize, infeasible_starting_point, sigma, initial_trust_radius, r_max, \
        #             noise_filename, solver_executable, print_iter_data, max_iter, \
        #             result_filename_header)
        #
        # print("\nTesting CompoStep_TR_MA")
        # result_filename_header = result_filename_folder + name[no] +"INF_CompoStep_"
        # algo1_TR_MA(model_name, perturbation_stepsize, infeasible_starting_point, sigma, initial_trust_radius, r_max, xi_N, \
        #             noise_filename, solver_executable, print_iter_data, max_iter, \
        #             result_filename_header)

# def draw_contour_global_convergence():
#     global_marker_size = 2
#     line_width = 1
#     max_iter = 20
#     fig = plt.figure(figsize=(9,6))
#     plt.subplot(231)
#     plot_contour("data/contour_data/res40_cost_plant.txt", "plant cost")
#     plt.title("")
#     plt.subplot(232)
#     plot_contour("data/contour_data/res40_cost_plant.txt", "plant cost")
#     plt.ylabel("")
#     plt.title("")
#     plt.subplot(233)
#     plot_contour("data/contour_data/res40_cost_plant.txt", "plant cost")
#     plt.title("")
#     plt.ylabel("")
#     plt.subplot(234)
#     plot_contour("data/contour_data/res40_con_plant.txt", "plant constraint")
#     plt.title("")
#     plt.subplot(235)
#     plot_contour("data/contour_data/res40_con_plant.txt", "plant constraint")
#     plt.ylabel("")
#     plt.title("")
#     plt.subplot(236)
#     plot_contour("data/contour_data/res40_con_plant.txt", "plant constraint")
#     plt.ylabel("")
#     plt.title("")
#
#     # for ax_no in[231,232,233]:
#     #     plt.subplot(ax_no)
#     #     plt.ylim([-2,2])
#     # for ax_no in[234,235,236]:
#     #     plt.subplot(ax_no)
#     #     plt.ylim([-2,2])
#
#     line_style={
#         "norminal":"-",
#         "obj_partial_wrong_curv":'-',
#         "con_wrong_curv2":"-",
#     }
#     for no,prefix in enumerate(["LPF","SPF","LPIF","SPIF"]):
#         compo_step_data = pandas.read_csv("data/escape/global_conv/"+ prefix+"_CompoStep_input_data.txt", \
#                                           index_col=0, header=0, sep='\t')
#         origin_ma = pandas.read_csv("data/escape/global_conv/"+ prefix+"_MA_input_data.txt", \
#                                        index_col=0, header=0, sep='\t')
#
#         plt.subplot(2,3,no+1)
#         for iter in range(max_iter - 1):
#             x = [compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter + 1, 'u1']]
#             y = [compo_step_data.loc[iter, 'u2'], compo_step_data.loc[iter + 1, 'u2']]
#             plt.plot(x, y, color="black", linewidth=line_width, linestyle=line_style[model_name])
#
#         for iter in range(max_iter):
#             plt.plot(compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter, 'u2'], \
#                      marker='o', c='black', markersize=global_marker_size)
#
#         for iter in range(max_iter - 1):
#             x = [origin_ma.loc[iter, 'u1'], origin_ma.loc[iter + 1, 'u1']]
#             y = [origin_ma.loc[iter, 'u2'], origin_ma.loc[iter + 1, 'u2']]
#             plt.plot(x, y, color="red", linewidth=line_width/3, linestyle=line_style[model_name])
#
#         for iter in range(max_iter):
#             plt.plot(origin_ma.loc[iter, 'u1'], origin_ma.loc[iter, 'u2'], \
#                      marker='o', c='red', markersize=global_marker_size/3)
#
#         plt.subplot(2,3,no+4)
#         for iter in range(max_iter - 1):
#             x = [compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter + 1, 'u1']]
#             y = [compo_step_data.loc[iter, 'u2'], compo_step_data.loc[iter + 1, 'u2']]
#             plt.plot(x, y, color="black", linewidth=line_width/3, linestyle=line_style[model_name])
#
#         for iter in range(max_iter):
#             plt.plot(compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter, 'u2'], \
#                      marker='o', c='black', markersize=global_marker_size/3)
#
#         for iter in range(max_iter - 1):
#             x = [origin_ma.loc[iter, 'u1'], origin_ma.loc[iter + 1, 'u1']]
#             y = [origin_ma.loc[iter, 'u2'], origin_ma.loc[iter + 1, 'u2']]
#             plt.plot(x, y, color="red", linewidth=line_width/3, linestyle=line_style[model_name])
#
#         for iter in range(max_iter):
#             plt.plot(origin_ma.loc[iter, 'u1'], origin_ma.loc[iter, 'u2'], \
#                      marker='o', c='red', markersize=global_marker_size/3)
#
#     plt.savefig("pic/escape/global_convengence", dpi=600)
#     plt.close()

def draw_profile_global_convergence():
    max_iter = 20
    global_marker_size = 2
    linewidth = 1
    fig = plt.figure(figsize=(21 * pic_constant, 10 * pic_constant))
    for ax_no in [231, 232, 233]:
        plt.subplot(ax_no)
        optimal = 1.453659e-01 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-1, 5])
        plt.grid(ls='--', axis='y')
    for ax_no in [234, 235, 236]:
        plt.subplot(ax_no)
        optimal = 0 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-2, 2])
        plt.grid(ls='--', axis='y')
    for no,prefix in enumerate(["norminal", "obj_partial_wrong_curv", "con_wrong_curv2"]):
        compo_step_data = pandas.read_csv("data/escape/global_conv/" + prefix + "_CompoStep_plant_data.txt", \
                                          index_col=0, header=0, sep='\t')
        origin_ma_data = pandas.read_csv("data/escape/global_conv/" + prefix + "_MA_plant_data.txt", \
                                          index_col=0, header=0, sep='\t')
        plt.subplot(2, 3, no + 1)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth, label="Composite-step")
        plt.plot(range(1, max_iter + 1), origin_ma_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='red', markersize=global_marker_size / 3, linewidth=linewidth / 3, label="Modifier adaptation")
        if no == 0:
            plt.ylabel("cost", font_label)
            plt.title(r"Model with correct curvatures", fontdict=font_title)
        elif no == 1:
            plt.title(r"Model with wrong objective curvature", fontdict=font_title)
        elif no == 2:
            plt.title(r"Model with wrong constraint curvature", fontdict=font_title)
            plt.legend(loc='upper right', prop=font_legend)
        plt.xticks([0, 5, 10, 15, 20])
        plt.yticks([-1, 1, 3, 5])
        plt.tick_params(labelsize=global_tick_size)

        plt.subplot(2, 3, no + 4)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'con'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
        plt.plot(range(1, max_iter + 1), origin_ma_data.loc[1:max_iter, 'con'], \
                 marker='o', c='red', markersize=global_marker_size / 3, linewidth=linewidth / 3)
        if no == 0:
            plt.ylabel("constraint", font_label)
        plt.xlabel("RTO iteration", font_label)
        plt.xticks([0, 5, 10, 15, 20])
        plt.tick_params(labelsize=global_tick_size)

    for i in range(1, 7):
        plt.subplot(2, 3, i)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    plt.savefig("pic/escape/global_convengence_profile", dpi=600)
    plt.close()

def draw_profile_penalty():
    max_iter = 20
    global_marker_size = 2
    linewidth = 1
    fig = plt.figure(figsize=(21 * pic_constant, 10 * pic_constant))
    for ax_no in [231,232,233]:
        plt.subplot(ax_no)
        optimal = 1.453659e-01 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-1,5])
        plt.grid(ls='--',axis='y')
    for ax_no in [234,235,236]:
        plt.subplot(ax_no)
        optimal = 0 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-2, 2])
        plt.grid(ls='--',axis='y')
    for no,prefix in enumerate(["LPF","MPF","SPF"]):
        compo_step_data = pandas.read_csv("data/escape/penalty_coeff/"+ prefix+"_CompoStep_plant_data.txt", \
                                          index_col=0, header=0, sep='\t')
        penalty_tr_data = pandas.read_csv("data/escape/penalty_coeff/"+ prefix+"_Penalty_plant_data.txt", \
                                       index_col=0, header=0, sep='\t')
        plt.subplot(2,3,no+1)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth,label="Composite-step")
        plt.plot(range(1, max_iter + 1), penalty_tr_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='red', markersize=global_marker_size/3, linewidth=linewidth/3,label="Penalty")
        if no ==0:
            plt.ylabel("cost",font_label)
            plt.title(r"$\sigma=1$", fontdict=font_title)
        elif no == 1:
            plt.title(r"$\sigma=0.1$", fontdict=font_title)
        elif no == 2:
            plt.title(r"$\sigma=0.01$", fontdict=font_title)
            plt.legend(loc='upper right', prop=font_legend)
        plt.xticks([0, 5, 10, 15, 20])
        plt.yticks([-1,1,3,5])
        plt.tick_params(labelsize=global_tick_size)

        plt.subplot(2,3, no + 4)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'con'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
        plt.plot(range(1, max_iter + 1), penalty_tr_data.loc[1:max_iter, 'con'], \
                 marker='o', c='red', markersize=global_marker_size/3, linewidth=linewidth/3)
        if no == 0:
            plt.ylabel("constraint",font_label)
        plt.xlabel("RTO iteration",font_label)
        plt.xticks([0, 5, 10, 15, 20])
        plt.tick_params(labelsize=global_tick_size)

    for i in range(1,7):
        plt.subplot(2,3,i)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("pic/escape/penalty_comparison", dpi=600)
    plt.close()

if __name__ == "__main__":
    # do_test_global_convergence()
    draw_profile_global_convergence()
    # do_test_penalty()
    draw_profile_penalty()