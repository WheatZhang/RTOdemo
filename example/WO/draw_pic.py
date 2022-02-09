import pandas
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import MultipleLocator

def draw_overall_pic():
    pic_filename="pic/overall_pic0.png"
    # =======================================
    #            Load Data
    # =======================================
    folder="data/0/"
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



if __name__ == "__main__":
    draw_overall_pic()