import pandas
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import MultipleLocator

def draw_overall_pic(folder, pic_filename):
    # pic_filename="pic/overall_pic0-0.png"
    # =======================================
    #            Load Data
    # =======================================
    # folder="data/0-0/"
    MA_plant_data = pandas.read_csv(folder + "MA_plant_data.txt", sep='\t', index_col=0, header=0)
    MA_input_data = pandas.read_csv(folder + "MA_input_data.txt", sep='\t', index_col=0, header=0)
    MAy_plant_data = pandas.read_csv(folder + "MAy_plant_data.txt", sep='\t', index_col=0, header=0)
    MAy_input_data = pandas.read_csv(folder + "MAy_input_data.txt", sep='\t', index_col=0, header=0)
    PE_plant_data = pandas.read_csv(folder + "PE_plant_data.txt", sep='\t', index_col=0, header=0)
    PE_input_data = pandas.read_csv(folder + "PE_input_data.txt", sep='\t', index_col=0, header=0)
    GPE_plant_data = pandas.read_csv(folder + "GPE_plant_data.txt", sep='\t', index_col=0, header=0)
    GPE_input_data = pandas.read_csv(folder + "GPE_input_data.txt", sep='\t', index_col=0, header=0)
    ISOPE_plant_data = pandas.read_csv(folder + "ISOPE_plant_data.txt", sep='\t', index_col=0, header=0)
    ISOPE_input_data = pandas.read_csv(folder + "ISOPE_input_data.txt", sep='\t', index_col=0, header=0)

    max_iter=20
    # #----------Optimal-------------
    ori_plant_profit = 190.98 * np.ones(max_iter+1)
    ori_FB = 4.7874 * np.ones(max_iter+1)
    ori_TR = 89.704 * np.ones(max_iter+1)
    # =======================================
    #            Draw Pictures
    # =======================================
    global_font_size = 20
    global_tick_size = 20
    global_linewidth = 1
    global_point_size = 4
    pic_constant = 0.39370
    font_factor = np.sqrt(1 / pic_constant)
    global_font_size = global_font_size / font_factor
    global_tick_size = global_tick_size / font_factor

    fig = plt.figure(figsize=(13 * 0.39370, 26 * 0.39370))
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_tick_size,
             }
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font3 = {'family': 'Times New Roman',
             'weight': 'bold',
             'size': global_font_size,
             }
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plt.subplot(3, 1, 1)
    Optimal, = plt.plot(range(max_iter+1), ori_plant_profit, linewidth=global_linewidth, label='Optimal', color='gray',
                        linestyle='--')
    MA, = plt.plot(range(max_iter), MA_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='MA', linestyle='-',
                   marker='o', markersize=global_point_size)
    MAy, = plt.plot(range(max_iter), MAy_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='MAy', linestyle='-',
                    marker='o', markersize=global_point_size)
    PA, = plt.plot(range(max_iter), PE_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='PE',
                   linestyle='-', marker='o', markersize=global_point_size)
    ISOPE, = plt.plot(range(max_iter), ISOPE_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='PE+MA',
                      linestyle='-', marker='o', markersize=global_point_size)
    GPE, = plt.plot(range(max_iter), GPE_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='GPE',
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.legend(handles=[MA, MAy, PA, ISOPE, GPE, Optimal], prop=font1)
    plt.legend(loc='lower right',
               fancybox=False,
               edgecolor='k',
               fontsize=global_point_size,
               shadow=False,
               facecolor='w',
               framealpha=1.0,
               prop={'family': 'Times New Roman', 'size': global_tick_size})
    plt.ylabel('(a) Plant Profit', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(145, 195)

    plt.subplot(3, 1, 2)
    OptimalFB, = plt.plot(range(max_iter + 1), ori_FB, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='MA', linestyle='-',
                   marker='o', markersize=global_point_size)
    MAy, = plt.plot(range(max_iter), MAy_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='MAy', linestyle='-',
                    marker='o', markersize=global_point_size)
    PA, = plt.plot(range(max_iter), PE_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='PE',
                   linestyle='-', marker='o', markersize=global_point_size)
    ISOPE, = plt.plot(range(max_iter), ISOPE_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='PE+MA',
                      linestyle='-',
                      marker='o', markersize=global_point_size)

    GPE, = plt.plot(range(max_iter), GPE_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='GPE',
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.ylabel('(b) Flowrates of B, $F_B$ (kg/s)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.yticks([4, 4.5, 5, 5.5])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(3.9, 5.6)

    plt.subplot(3, 1, 3)
    OptimalTR, = plt.plot(range(max_iter + 1), ori_TR, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='MA', linestyle='-',
                   marker='o', markersize=global_point_size)
    MAy, = plt.plot(range(max_iter), MAy_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='MAy', linestyle='-',
                    marker='o', markersize=global_point_size)
    PA, = plt.plot(range(max_iter), PE_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='PE',
                   linestyle='-', marker='o', markersize=global_point_size)
    ISOPE, = plt.plot(range(max_iter), ISOPE_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='PE+MA',
                      linestyle='-',
                      marker='o', markersize=global_point_size)
    GPE, = plt.plot(range(max_iter), GPE_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='GPE',
                      linestyle='-', marker='o', markersize=global_point_size)

    # plt.legend(handles=[OptimalTR, MA, MAy, PA, GPE, ISOPE], prop=font1)
    plt.xlabel('RTO Iterations', font2)
    plt.ylabel('(c) CSTR Temperature, $T_R$ (℃)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax = plt.gca()
    plt.tick_params(labelsize=global_tick_size)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(74, 101)
    y_major_locator = MultipleLocator(5)
    ax = plt.gca()
    ax.yaxis.set_major_locator(y_major_locator)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_filename, dpi=600)
    plt.close()


def draw_GPE_different_weight():
    pic_filename="pic/gpe-different-weight.png"
    # =======================================
    #            Load Data
    # =======================================
    weights=[0.1,1,10,100,1000]
    plant_data = []
    input_data = []
    for i, n in enumerate(weights):
        plant_data.append(pandas.read_csv("data/batch_GPE_%s/GPE_plant_data.txt"%str(n), \
                                               sep='\t', index_col=0, header=0))
        input_data.append(pandas.read_csv("data/batch_GPE_%s/GPE_input_data.txt" % str(n), \
                                               sep='\t', index_col=0, header=0))

    max_iter=20
    # #----------Optimal-------------
    ori_plant_profit = 190.98 * np.ones(max_iter+1)
    ori_FB = 4.7874 * np.ones(max_iter+1)
    ori_TR = 89.704 * np.ones(max_iter+1)
    # =======================================
    #            Draw Pictures
    # =======================================
    global_font_size = 20
    global_tick_size = 20
    global_linewidth = 1
    global_point_size = 4
    pic_constant = 0.39370
    font_factor = np.sqrt(1 / pic_constant)
    global_font_size = global_font_size / font_factor
    global_tick_size = global_tick_size / font_factor

    fig = plt.figure(figsize=(13 * 0.39370, 26 * 0.39370))
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_tick_size,
             }
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font3 = {'family': 'Times New Roman',
             'weight': 'bold',
             'size': global_font_size,
             }
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plt.subplot(3, 1, 1)
    Optimal, = plt.plot(range(max_iter+1), ori_plant_profit, linewidth=global_linewidth, label='Optimal', color='gray',
                        linestyle='--')
    GPE1, = plt.plot(range(max_iter), plant_data[0].loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label=str(weights[0]), linestyle='-',
                   marker='o', markersize=global_point_size)
    GPE2, = plt.plot(range(max_iter), plant_data[1].loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label=str(weights[1]), linestyle='-',
                    marker='o', markersize=global_point_size)
    GPE3, = plt.plot(range(max_iter), plant_data[2].loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label=str(weights[2]),
                   linestyle='-', marker='o', markersize=global_point_size)
    GPE4, = plt.plot(range(max_iter), plant_data[3].loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label=str(weights[3]),
                      linestyle='-', marker='o', markersize=global_point_size)
    GPE5, = plt.plot(range(max_iter), plant_data[4].loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label=str(weights[4]),
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.legend(handles=[GPE1, GPE2, GPE3, GPE4, GPE5, Optimal], prop=font1)
    plt.legend(loc='lower right',
               fancybox=False,
               edgecolor='k',
               fontsize=global_point_size,
               shadow=False,
               facecolor='w',
               framealpha=1.0,
               prop={'family': 'Times New Roman', 'size': global_tick_size})
    plt.ylabel('(a) Plant Profit', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(145, 195)

    plt.subplot(3, 1, 2)
    Optimal, = plt.plot(range(max_iter + 1), ori_FB, linewidth=global_linewidth, label='Optimal',
                        color='gray',
                        linestyle='--')
    GPE1, = plt.plot(range(max_iter), input_data[0].loc[1:(max_iter + 1), "Fb"],
                     linewidth=global_linewidth, label=str(weights[0]), linestyle='-',
                     marker='o', markersize=global_point_size)
    GPE2, = plt.plot(range(max_iter), input_data[1].loc[1:(max_iter + 1), "Fb"],
                     linewidth=global_linewidth, label=str(weights[1]), linestyle='-',
                     marker='o', markersize=global_point_size)
    GPE3, = plt.plot(range(max_iter), input_data[2].loc[1:(max_iter + 1), "Fb"],
                     linewidth=global_linewidth, label=str(weights[2]),
                     linestyle='-', marker='o', markersize=global_point_size)
    GPE4, = plt.plot(range(max_iter), input_data[3].loc[1:(max_iter + 1), "Fb"],
                     linewidth=global_linewidth, label=str(weights[3]),
                     linestyle='-', marker='o', markersize=global_point_size)
    GPE5, = plt.plot(range(max_iter), input_data[4].loc[1:(max_iter + 1), "Fb"],
                     linewidth=global_linewidth, label=str(weights[4]),
                     linestyle='-', marker='o', markersize=global_point_size)

    plt.ylabel('(b) Flowrates of B, $F_B$ (kg/s)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.yticks([4, 4.5, 5, 5.5])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(3.9, 5.6)

    plt.subplot(3, 1, 3)
    OptimalTR, = plt.plot(range(max_iter + 1), ori_TR, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    GPE1, = plt.plot(range(max_iter), input_data[0].loc[1:(max_iter + 1), "Tr"],
                     linewidth=global_linewidth, label=str(weights[0]), linestyle='-',
                     marker='o', markersize=global_point_size)
    GPE2, = plt.plot(range(max_iter), input_data[1].loc[1:(max_iter + 1), "Tr"],
                     linewidth=global_linewidth, label=str(weights[1]), linestyle='-',
                     marker='o', markersize=global_point_size)
    GPE3, = plt.plot(range(max_iter), input_data[2].loc[1:(max_iter + 1), "Tr"],
                     linewidth=global_linewidth, label=str(weights[2]),
                     linestyle='-', marker='o', markersize=global_point_size)
    GPE4, = plt.plot(range(max_iter), input_data[3].loc[1:(max_iter + 1), "Tr"],
                     linewidth=global_linewidth, label=str(weights[3]),
                     linestyle='-', marker='o', markersize=global_point_size)
    GPE5, = plt.plot(range(max_iter), input_data[4].loc[1:(max_iter + 1), "Tr"],
                     linewidth=global_linewidth, label=str(weights[4]),
                     linestyle='-', marker='o', markersize=global_point_size)

    plt.xlabel('RTO Iterations', font2)
    plt.ylabel('(c) CSTR Temperature, $T_R$ (℃)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax = plt.gca()
    plt.tick_params(labelsize=global_tick_size)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(74, 101)
    y_major_locator = MultipleLocator(5)
    ax = plt.gca()
    ax.yaxis.set_major_locator(y_major_locator)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_filename, dpi=600)
    plt.close()

def draw_para_set_comp_pic(folder, pic_filename):
    # pic_filename="pic/overall_pic0-0.png"
    # =======================================
    #            Load Data
    # =======================================
    # folder="data/0-0/"
    PE_plant_data = pandas.read_csv(folder + "PE_plant_data.txt", sep='\t', index_col=0, header=0)
    PE_input_data = pandas.read_csv(folder + "PE_input_data.txt", sep='\t', index_col=0, header=0)
    GPE_plant_data = pandas.read_csv(folder + "GPE_plant_data.txt", sep='\t', index_col=0, header=0)
    GPE_input_data = pandas.read_csv(folder + "GPE_input_data.txt", sep='\t', index_col=0, header=0)
    ISOPE_plant_data = pandas.read_csv(folder + "ISOPE_plant_data.txt", sep='\t', index_col=0, header=0)
    ISOPE_input_data = pandas.read_csv(folder + "ISOPE_input_data.txt", sep='\t', index_col=0, header=0)

    max_iter=20
    # #----------Optimal-------------
    ori_plant_profit = 190.98 * np.ones(max_iter+1)
    ori_FB = 4.7874 * np.ones(max_iter+1)
    ori_TR = 89.704 * np.ones(max_iter+1)
    # =======================================
    #            Draw Pictures
    # =======================================
    global_font_size = 20
    global_tick_size = 20
    global_linewidth = 1
    global_point_size = 4
    pic_constant = 0.39370
    font_factor = np.sqrt(1 / pic_constant)
    global_font_size = global_font_size / font_factor
    global_tick_size = global_tick_size / font_factor

    fig = plt.figure(figsize=(13 * 0.39370, 26 * 0.39370))
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_tick_size,
             }
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font3 = {'family': 'Times New Roman',
             'weight': 'bold',
             'size': global_font_size,
             }
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plt.subplot(3, 1, 1)
    Optimal, = plt.plot(range(max_iter+1), ori_plant_profit, linewidth=global_linewidth, label='Optimal', color='gray',
                        linestyle='--')
    PA, = plt.plot(range(max_iter), PE_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='PE',
                   linestyle='-', marker='o', markersize=global_point_size)
    ISOPE, = plt.plot(range(max_iter), ISOPE_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='PE+MA',
                      linestyle='-', marker='o', markersize=global_point_size)
    GPE, = plt.plot(range(max_iter), GPE_plant_data.loc[1:(max_iter+1),"profit"].multiply(-1), linewidth=global_linewidth, label='GPE',
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.legend(handles=[PA, ISOPE, GPE, Optimal], prop=font1)
    plt.legend(loc='lower right',
               fancybox=False,
               edgecolor='k',
               fontsize=global_point_size,
               shadow=False,
               facecolor='w',
               framealpha=1.0,
               prop={'family': 'Times New Roman', 'size': global_tick_size})
    plt.ylabel('(a) Plant Profit', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(145, 195)

    plt.subplot(3, 1, 2)
    OptimalFB, = plt.plot(range(max_iter + 1), ori_FB, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    PA, = plt.plot(range(max_iter), PE_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='PE',
                   linestyle='-', marker='o', markersize=global_point_size)
    ISOPE, = plt.plot(range(max_iter), ISOPE_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='PE+MA',
                      linestyle='-',
                      marker='o', markersize=global_point_size)

    GPE, = plt.plot(range(max_iter), GPE_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='GPE',
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.ylabel('(b) Flowrates of B, $F_B$ (kg/s)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.yticks([4, 4.5, 5, 5.5])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(3.9, 5.6)

    plt.subplot(3, 1, 3)
    OptimalTR, = plt.plot(range(max_iter + 1), ori_TR, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    PA, = plt.plot(range(max_iter), PE_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='PE',
                   linestyle='-', marker='o', markersize=global_point_size)
    ISOPE, = plt.plot(range(max_iter), ISOPE_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='PE+MA',
                      linestyle='-',
                      marker='o', markersize=global_point_size)
    GPE, = plt.plot(range(max_iter), GPE_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='GPE',
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.xlabel('RTO Iterations', font2)
    plt.ylabel('(c) CSTR Temperature, $T_R$ (℃)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax = plt.gca()
    plt.tick_params(labelsize=global_tick_size)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(74, 101)
    y_major_locator = MultipleLocator(5)
    ax = plt.gca()
    ax.yaxis.set_major_locator(y_major_locator)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_filename, dpi=600)
    plt.close()

def generate_contour_data(x,y,func,filename):
    x_len = x.shape[0]
    y_len = y.shape[0]
    with open(filename, 'w') as fp:
        for i in range(x_len):
            if i != 0:
                fp.write('\n')
            for j in range(y_len):
                if j != 0:
                    fp.write('\t')
                fp.write("%.6e" % func(x[i], y[j]))
        fp.write('\n')


def plot_contour(data_file,title):
    global_font_size = 20
    global_marker_size = 3
    contour_linewidth = 1.5
    global_tick_size = 15
    global_contour_label_size = 18
    font_factor = np.sqrt(1 / 0.39370)
    global_font_size = global_font_size / font_factor
    contour_linewidth = contour_linewidth / font_factor
    global_contour_label_size = global_contour_label_size / font_factor

    font_title = {'family': 'Times New Roman',
             'size': global_font_size,
             }

    data_stored = pandas.read_csv(data_file, index_col=None, header=None, sep='\t')
    x_len = data_stored.shape[0]
    y_len = data_stored.shape[1]
    Z = np.zeros(shape=(x_len, y_len))

    for i in range(x_len):
        for j in range(y_len):
            Z[j, i] = data_stored.iloc[i,j]

    a=np.linspace(start=3.01,stop=6.01,num=x_len,endpoint=True)
    b=np.linspace(start=70, stop=100, num=y_len,endpoint=True)
    N = np.array([20,60,100,140,160,180,190,200,210])  # 用来指明等高线对应的值为多少时才出图线
    CS = plt.contour(a,b,Z, N, linewidths=contour_linewidth, cmap=plt.get_cmap('jet'))  # 画出等高线图，cmap表示颜色的图层。
    plt.tick_params(top= ['on', 'true'], right= ['on', 'true'], which='both')
    plt.clabel(CS, inline=True, fmt='%d', fontsize=global_contour_label_size)  # 在等高线图里面加入每条线对应的值
    plt.title(title,fontdict=font_title)


def plot_model_maintenance(algo_results_folder, datafile_folder, pic_filename):
    from rtolib.core.solve import PyomoOptimizer
    from rtolib.model.wo_reactor import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor, default_WOR_description
    from pyomo.environ import SolverFactory
    data_files={"Plant":datafile_folder+"_Plant.txt",
                "Model":datafile_folder+"_Model.txt",
                "PE":datafile_folder+"_PE.txt",
                "ISOPE":datafile_folder+"_ISOPE.txt",
                "GPE":datafile_folder+"_GPE.txt",}

    plant_optimizer = PyomoOptimizer(RTO_Plant_WO_reactor())
    model_optimizer = PyomoOptimizer(RTO_Mismatched_WO_reactor())
    plant_optimizer.build(default_WOR_description)
    model_optimizer.build(default_WOR_description)
    default_options = {'max_iter': 100}
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    plant_optimizer.set_solver(solver, tee=False, default_options=default_options)
    model_optimizer.set_solver(solver, tee=False, default_options=default_options)
    #--------read data--------------------
    PE_model_data = pandas.read_csv(algo_results_folder + "PE_model_data.txt", sep='\t', index_col=0, header=0)
    GPE_model_data = pandas.read_csv(algo_results_folder + "GPE_model_data.txt", sep='\t', index_col=0, header=0)
    ISOPE_model_data = pandas.read_csv(algo_results_folder + "ISOPE_model_data.txt", sep='\t', index_col=0, header=0)

    # =======================================
    #            Draw Pictures
    # =======================================

    pic_dpi = 600

    global_font_size = 20
    global_marker_size = 3
    contour_linewidth = 1.5
    global_tick_size = 15
    global_contour_label_size = 18

    pic_constant = 0.39370
    font_factor = np.sqrt(1 / 0.39370)
    global_font_size = global_font_size / font_factor
    global_marker_size = global_marker_size / font_factor
    contour_linewidth = contour_linewidth / font_factor
    global_tick_size = global_tick_size / font_factor
    global_contour_label_size = global_contour_label_size / font_factor
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font_point = {'family': 'Times New Roman',
                  'weight': 'normal',
                  'size': global_contour_label_size,
                  }
    fig = plt.figure(figsize=(26 * pic_constant, 8 * pic_constant))
    print('Plant:')
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.subplot(1, 5, 1)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.xlabel('Flowrate of B (kg/s)', font2)
    plt.ylabel('CSTR Temperature, $T_R$ (℃)', font2)
    plt.tick_params(labelsize=global_tick_size)
    plot_contour(data_files['Plant'],  '(a) Plant')

    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]

    optimized_input, solve_status = plant_optimizer.optimize(input_values={},
                                          param_values=None,
                                          use_homo=True)
    Fb=optimized_input['Fb']
    Tr=optimized_input['Tr']
    print(Fb)
    print(Tr)
    x_pos = Fb
    y_pos = Tr
    plt.plot(x_pos, y_pos, marker='o', c='black', markersize=global_marker_size)
    plt.text(x_pos, y_pos - 0.1, '%.2f, %.1f' % (Fb, Tr),  fontdict=font_point,
             verticalalignment="top",  # ‘center’ | ‘top’ | ‘bottom’
             horizontalalignment="center")  # ‘center’ | ‘right’ | ‘left’
    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    y_major_locator = MultipleLocator(5)
    ax.yaxis.set_major_locator(y_major_locator)

    print('Model:')
    plt.subplot(1, 5, 2)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plot_contour(data_files['Model'],  '(b) Model')
    plt.xlabel('Flowrate of B (kg/s)', font2)
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    optimized_input, solve_status = model_optimizer.optimize(input_values={},
                                                             param_values=None,
                                                             use_homo=True)
    Fb = optimized_input['Fb']
    Tr = optimized_input['Tr']
    print(Fb)
    print(Tr)
    x_pos = Fb
    y_pos = Tr
    plt.plot(x_pos, y_pos, marker='o', c='black', markersize=global_marker_size)
    plt.text(x_pos, y_pos - 0.1, '%.2f, %.1f' % (Fb, Tr),  fontdict=font_point,
             verticalalignment="top",  # ‘center’ | ‘top’ | ‘bottom’
             horizontalalignment="center")  # ‘center’ | ‘right’ | ‘left’
    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    y_major_locator = MultipleLocator(5)
    ax.yaxis.set_major_locator(y_major_locator)

    print('PA:')
    plt.subplot(1, 5, 3)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plot_contour(data_files['PE'],  '(c) PE')
    plt.xlabel('Flowrate of B (kg/s)', font2)
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    param_values={}
    for k in model_optimizer.pyomo_model.parameters.keys():
        param_values[k]=PE_model_data.loc[20, k]
    # print(param_values)
    optimized_input, solve_status = model_optimizer.optimize(input_values={},
                                                             param_values=param_values,
                                                             use_homo=True)
    Fb = optimized_input['Fb']
    Tr = optimized_input['Tr']
    print(Fb)
    print(Tr)
    x_pos = Fb
    y_pos = Tr
    plt.plot(x_pos, y_pos, marker='o', c='black', markersize=global_marker_size)
    plt.text(x_pos, y_pos - 0.1, '%.2f, %.1f' % (Fb, Tr),  fontdict=font_point,
             verticalalignment="top",  # ‘center’ | ‘top’ | ‘bottom’
             horizontalalignment="center")  # ‘center’ | ‘right’ | ‘left’
    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    y_major_locator = MultipleLocator(5)
    ax.yaxis.set_major_locator(y_major_locator)

    print('ISOPE:')
    plt.subplot(1, 5, 4)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plot_contour(data_files['ISOPE'],  '(d) PE+MA')
    plt.xlabel('Flowrate of B (kg/s)', font2)
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    param_values = {}
    for k in model_optimizer.pyomo_model.parameters.keys():
        param_values[k] = ISOPE_model_data.loc[20, k]
    optimized_input, solve_status = model_optimizer.optimize(input_values={},
                                                             param_values=param_values,
                                                             use_homo=True)
    Fb = optimized_input['Fb']
    Tr = optimized_input['Tr']
    print(Fb)
    print(Tr)
    x_pos = Fb
    y_pos = Tr
    plt.plot(x_pos, y_pos, marker='o', c='black', markersize=global_marker_size)
    plt.text(x_pos, y_pos - 0.1, '%.2f, %.1f' % (Fb, Tr),  fontdict=font_point,
             verticalalignment="top",  # ‘center’ | ‘top’ | ‘bottom’
             horizontalalignment="center")  # ‘center’ | ‘right’ | ‘left’
    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    y_major_locator = MultipleLocator(5)
    ax.yaxis.set_major_locator(y_major_locator)

    print('IPAMA:')
    plt.subplot(1, 5, 5)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plot_contour(data_files['GPE'],  '(e) GPE')
    plt.xlabel('Flowrate of B (kg/s)', font2)
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    param_values = {}
    for k in model_optimizer.pyomo_model.parameters.keys():
        param_values[k] = GPE_model_data.loc[20, k]
    optimized_input, solve_status = model_optimizer.optimize(input_values={},
                                                             param_values=param_values,
                                                             use_homo=True)
    Fb = optimized_input['Fb']
    Tr = optimized_input['Tr']
    print(Fb)
    print(Tr)
    x_pos = Fb
    y_pos = Tr
    plt.plot(x_pos, y_pos, marker='o', c='black', markersize=global_marker_size)
    plt.text(x_pos, y_pos - 0.1, '%.2f, %.1f' % (Fb, Tr),  fontdict=font_point,
             verticalalignment="top",  # ‘center’ | ‘top’ | ‘bottom’
             horizontalalignment="center")  # ‘center’ | ‘right’ | ‘left’
    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    y_major_locator = MultipleLocator(5)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.tight_layout(pad=0.1, w_pad=0, h_pad=0.1)
    plt.savefig(pic_filename, dpi=pic_dpi)
    plt.close()

if __name__ == "__main__":
    draw_overall_pic("data/1/", "pic/overall_pic1.png")