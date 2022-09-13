import matplotlib.pyplot as plt
from draw_lib import generate_all_contour_data, plot_contour
import pandas
import os
from algo_tr_ma import do_all_batches_for_all_model

def draw_nominal_contour(resolution):
    datafile_folder = "data/contour_data/res%d" % resolution
    data_files = {"cost_plant": datafile_folder + "_cost_plant.txt",
                  "con_plant": datafile_folder + "_con_plant.txt",
                  "cost_model": datafile_folder + "_cost_model.txt",
                  "con_model": datafile_folder + "_con_model.txt",
                  }
    for title,data_file in data_files.items():
        fig = plt.figure()
        plot_contour(data_file, title)
        plt.savefig("pic/contour/"+title+"_%d"%resolution, dpi=600)
        plt.close()

def draw_all_contour():
    datafile_folder = "data/contour_data"
    for root, dirs, files in os.walk(datafile_folder):
        for filename in files:
            if filename.endswith(".txt"):
                fig = plt.figure()
                title = filename.split(".")[0]
                plot_contour(os.path.join(root, filename), title)
                plt.savefig("pic/contour/" + title + ".png", dpi=600)
                plt.close()

def draw_algo_12_input_on_contour(data_prefix):
    global_marker_size = 2
    line_width = 1
    max_iter = 10
    compo_step_data = pandas.read_csv("data/batch12/"+data_prefix+"/CompoStep_TR_MA_input_data.txt", \
                                      index_col=0, header=0, sep='\t')
    penalty_data = pandas.read_csv("data/batch12/"+data_prefix+"/Penalty_TR_MA_input_data.txt", \
                                      index_col=0, header=0, sep='\t')
    fig = plt.figure()
    plt.subplot(121)
    plot_contour("data/contour_data/res40_cost_plant.txt", "plant cost")

    for iter in range(max_iter-1):
        x = [compo_step_data.loc[iter,'u1'], compo_step_data.loc[iter+1,'u1']]
        y = [compo_step_data.loc[iter,'u2'], compo_step_data.loc[iter+1,'u2']]
        plt.plot(x, y, color="black", linewidth=line_width)

    for iter in range(max_iter):
        plt.plot(compo_step_data.loc[iter,'u1'], compo_step_data.loc[iter,'u2'],\
                 marker='o', c='black', markersize=global_marker_size)

    for iter in range(max_iter-1):
        x = [penalty_data.loc[iter, 'u1'], penalty_data.loc[iter + 1, 'u1']]
        y = [penalty_data.loc[iter, 'u2'], penalty_data.loc[iter + 1, 'u2']]
        plt.plot(x, y, color="red", linewidth=line_width/3)

    for iter in range(max_iter):
        plt.plot(penalty_data.loc[iter, 'u1'], penalty_data.loc[iter, 'u2'], \
                 marker='o', c='red', markersize=global_marker_size/3)

    plt.subplot(122)
    plot_contour("data/contour_data/res40_con_plant.txt", "plant constraint")

    for iter in range(max_iter - 1):
        x = [compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter + 1, 'u1']]
        y = [compo_step_data.loc[iter, 'u2'], compo_step_data.loc[iter + 1, 'u2']]
        plt.plot(x, y, color="black", linewidth=line_width)

    for iter in range(max_iter):
        plt.plot(compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter, 'u2'], \
                 marker='o', c='black', markersize=global_marker_size)

    for iter in range(max_iter - 1):
        x = [penalty_data.loc[iter, 'u1'], penalty_data.loc[iter + 1, 'u1']]
        y = [penalty_data.loc[iter, 'u2'], penalty_data.loc[iter + 1, 'u2']]
        plt.plot(x, y, color="red", linewidth=line_width / 3)

    for iter in range(max_iter):
        plt.plot(penalty_data.loc[iter, 'u1'], penalty_data.loc[iter, 'u2'], \
                 marker='o', c='red', markersize=global_marker_size / 3)
    plt.ylabel("")

    plt.savefig("pic/batch12/"+data_prefix+"_algo12_input_iterates", dpi=600)
    plt.close()

def draw_algo_13_input_on_contour(data_prefix):
    global_marker_size = 2
    line_width = 1
    max_iter = 10
    compo_step_data = pandas.read_csv("data/batch13/"+data_prefix+"/CompoStep_TR_MA_input_data.txt", \
                                      index_col=0, header=0, sep='\t')
    inf_averse_data = pandas.read_csv("data/batch13/"+data_prefix+"/InfAverse_TR_MA_input_data.txt", \
                                      index_col=0, header=0, sep='\t')
    fig = plt.figure()
    plt.subplot(121)
    plot_contour("data/contour_data/res40_cost_plant.txt", "plant cost")

    for iter in range(max_iter-1):
        x = [compo_step_data.loc[iter,'u1'], compo_step_data.loc[iter+1,'u1']]
        y = [compo_step_data.loc[iter,'u2'], compo_step_data.loc[iter+1,'u2']]
        plt.plot(x, y, color="black", linewidth=line_width)

    for iter in range(max_iter):
        plt.plot(compo_step_data.loc[iter,'u1'], compo_step_data.loc[iter,'u2'],\
                 marker='o', c='black', markersize=global_marker_size)

    for iter in range(max_iter-1):
        x = [inf_averse_data.loc[iter, 'u1'], inf_averse_data.loc[iter + 1, 'u1']]
        y = [inf_averse_data.loc[iter, 'u2'], inf_averse_data.loc[iter + 1, 'u2']]
        plt.plot(x, y, color="darkorange", linewidth=line_width/3)

    for iter in range(max_iter):
        plt.plot(inf_averse_data.loc[iter, 'u1'], inf_averse_data.loc[iter, 'u2'], \
                 marker='o', c='darkorange', markersize=global_marker_size/3)

    plt.subplot(122)
    plot_contour("data/contour_data/res40_con_plant.txt", "plant constraint")

    for iter in range(max_iter - 1):
        x = [compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter + 1, 'u1']]
        y = [compo_step_data.loc[iter, 'u2'], compo_step_data.loc[iter + 1, 'u2']]
        plt.plot(x, y, color="black", linewidth=line_width)

    for iter in range(max_iter):
        plt.plot(compo_step_data.loc[iter, 'u1'], compo_step_data.loc[iter, 'u2'], \
                 marker='o', c='black', markersize=global_marker_size)

    for iter in range(max_iter - 1):
        x = [inf_averse_data.loc[iter, 'u1'], inf_averse_data.loc[iter + 1, 'u1']]
        y = [inf_averse_data.loc[iter, 'u2'], inf_averse_data.loc[iter + 1, 'u2']]
        plt.plot(x, y, color="darkorange", linewidth=line_width / 3)

    for iter in range(max_iter):
        plt.plot(inf_averse_data.loc[iter, 'u1'], inf_averse_data.loc[iter, 'u2'], \
                 marker='o', c='darkorange', markersize=global_marker_size / 3)
    plt.ylabel("")

    plt.savefig("pic/batch13/"+data_prefix+"_algo13_input_iterates", dpi=600)
    plt.close()

def plot_algo12_profile(data_prefix):
    max_iter=20
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data = pandas.read_csv("data/batch12/"+data_prefix+"/CompoStep_TR_MA_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    penalty_plant_data = pandas.read_csv("data/batch12/"+data_prefix+"/Penalty_TR_MA_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')
    compo_step_model_data = pandas.read_csv("data/batch12/"+data_prefix+"/CompoStep_TR_MA_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    penalty_model_data = pandas.read_csv("data/batch12/"+data_prefix+"/Penalty_TR_MA_model_data.txt", \
                                   index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,9))
    plt.subplot(411)
    plt.plot(range(1,max_iter+1), compo_step_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), penalty_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant cost")
    plt.subplot(412)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), penalty_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(413)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), penalty_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$\rho$")
    plt.subplot(414)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), penalty_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("radius")
    plt.xlabel("#iteration")

    plt.savefig("pic/batch12/"+data_prefix+"_algo12_profile", dpi=600)
    plt.close()

def plot_algo13_profile(data_prefix):
    max_iter=20
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data = pandas.read_csv("data/batch13/"+data_prefix+"/CompoStep_TR_MA_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    inf_averse_plant_data = pandas.read_csv("data/batch13/"+data_prefix+"/InfAverse_TR_MA_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')
    compo_step_model_data = pandas.read_csv("data/batch13/"+data_prefix+"/CompoStep_TR_MA_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    inf_averse_model_data = pandas.read_csv("data/batch13/"+data_prefix+"/InfAverse_TR_MA_model_data.txt", \
                                   index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,9))
    plt.subplot(411)
    plt.plot(range(1,max_iter+1), compo_step_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='darkorange', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant cost")
    plt.subplot(412)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='darkorange', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(413)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='darkorange', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$\rho$")
    plt.subplot(414)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='darkorange', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("radius")
    plt.xlabel("#iteration")

    plt.savefig("pic/batch13/"+data_prefix+"_algo13_profile", dpi=600)
    plt.close()


def draw_batch_pic_algo12():
    batch_prefixes = ["N_", "O1_", "O2_", "C_", "OC_"]
    batch_names = [
        "U0","U1","U2","U3","U4","U5","U6","U7",
        "R0","R1","R2","R3","R4",
        "FP0","FP1","FP2","FP3","FP4",
        "IFP0","IFP1","IFP2","IFP3","IFP4",
        "XINM0","XINM1","XINM2",
        "XISR0","XISR1","XISR2",
    ]
    for batch_prefix in batch_prefixes:
        for data_prefix in batch_names:
            draw_algo_12_input_on_contour(batch_prefix+data_prefix)
            plot_algo12_profile(batch_prefix+data_prefix)

def draw_batch_pic_algo13():
    batch_prefixes = ["N_", "O1_", "O2_", "C_", "OC_"]
    batch_names = [
        "U0","U1","U2","U3","U4",
    ]
    for batch_prefix in batch_prefixes:
        for data_prefix in batch_names:
            draw_algo_13_input_on_contour(batch_prefix+data_prefix)
            plot_algo13_profile(batch_prefix+data_prefix)


if __name__ == "__main__":
    # resolution=40
    # generate_all_contour_data(resolution)
    # draw_nominal_contour(resolution)
    # draw_all_contour()
    #
    # do_all_batches_for_all_model()
    # draw_batch_pic_algo12()
    draw_batch_pic_algo13()