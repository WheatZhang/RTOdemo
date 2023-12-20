import numpy as np
import matplotlib.pyplot as plt
import pandas
from matplotlib.pyplot import MultipleLocator


def draw_compare_unconstrained_wo_adchem(folder, pic_filename):
    # pic_filename="pic/overall_pic0-0.png"
    # =======================================
    #            Load Data
    # =======================================
    # folder="data/0-0/"
    MA_plant_data = pandas.read_csv(folder + "MA_plant_data.txt", sep='\t', index_col=0, header=0)
    MA_input_data = pandas.read_csv(folder + "MA_input_data.txt", sep='\t', index_col=0, header=0)
    CPWL_MA_plant_data = pandas.read_csv(folder + "CPWL_MA_plant_data.txt", sep='\t', index_col=0, header=0)
    CPWL_MA_input_data = pandas.read_csv(folder + "CPWL_MA_input_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_Subgrad_plant_data = pandas.read_csv(folder + "QCPWL_Subgrad_plant_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_Subgrad_input_data = pandas.read_csv(folder + "QCPWL_Subgrad_input_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_MINLP_plant_data = pandas.read_csv(folder + "QCPWL_MINLP_plant_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_MINLP_input_data = pandas.read_csv(folder + "QCPWL_MINLP_input_data.txt", sep='\t', index_col=0, header=0)

    max_iter=20
    # #----------Optimal-------------
    ori_plant_cost = -190.98 * np.ones(max_iter+1)
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

    fig = plt.figure(figsize=(13 * 0.39370, 18 * 0.39370))
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
    Optimal, = plt.plot(range(max_iter+1), ori_plant_cost, linewidth=global_linewidth, label='Optimal', color='gray',
                        linestyle='--')
    MA, = plt.plot(range(max_iter), MA_plant_data.loc[1:(max_iter+1),"cost"].multiply(1), linewidth=global_linewidth, label='Nonlinear model', linestyle='dashdot',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_plant_data.loc[1:(max_iter+1),"cost"].multiply(1), linewidth=global_linewidth, label='DCA+CPWL', linestyle='dotted',
                    marker='o', markersize=global_point_size)
    QCPWL_Subgrad, = plt.plot(range(max_iter), QCPWL_Subgrad_plant_data.loc[1:(max_iter+1), "cost"].multiply(1),
                              linewidth=global_linewidth,
                              label='DCA+QCPWL',
                              linestyle='-', marker='o', markersize=global_point_size)

    plt.legend(handles=[MA, CPWL_MA, QCPWL_Subgrad, Optimal], prop=font1)
    plt.legend(loc='upper right',
               fancybox=False,
               edgecolor='k',
               fontsize=global_point_size,
               shadow=False,
               facecolor='w',
               framealpha=1.0,
               prop={'family': 'Times New Roman', 'size': global_tick_size})
    plt.ylabel('(a) Plant cost', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(-195, -145)

    plt.subplot(3, 1, 2)
    OptimalFB, = plt.plot(range(max_iter + 1), ori_FB, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='Nonlinear model', linestyle='dashdot',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='DCA+CPWL', linestyle='dotted',
                    marker='o', markersize=global_point_size)
    QCPWL_Subgrad, = plt.plot(range(max_iter), QCPWL_Subgrad_input_data.loc[range(max_iter), "Fb"],
                              linewidth=global_linewidth,
                              label='DCA+QPWL',
                              linestyle='-', marker='o', markersize=global_point_size)

    plt.ylabel('(b) Flowrate of B $F_B$ (kg/s)', font2)
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
    plt.ylim(3.9, 5)

    plt.subplot(3, 1, 3)
    OptimalTR, = plt.plot(range(max_iter + 1), ori_TR, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='Nonlinear model', linestyle='dashdot',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='DCA+CPWL', linestyle='dashed',
                    marker='o', markersize=global_point_size)
    QCPWL_Subgrad, = plt.plot(range(max_iter), QCPWL_Subgrad_input_data.loc[range(max_iter), "Tr"], linewidth=global_linewidth,
                         label='DCA+QPWL',
                         linestyle='-', marker='o', markersize=global_point_size)

    # plt.legend(handles=[OptimalTR, MA, CPWL_MA, PA, QCPWL_MA, ISOQCPWL_MA], prop=font1)
    plt.xlabel('RTO Iterations', font2)
    plt.ylabel('(c) Temperature $T_R$ (℃)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax = plt.gca()
    plt.tick_params(labelsize=global_tick_size)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(74, 95)
    y_major_locator = MultipleLocator(5)
    ax = plt.gca()
    ax.yaxis.set_major_locator(y_major_locator)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_filename, dpi=300)
    plt.close()
def draw_compare_unconstrained_wo(folder, pic_filename):
    # pic_filename="pic/overall_pic0-0.png"
    # =======================================
    #            Load Data
    # =======================================
    # folder="data/0-0/"
    MA_plant_data = pandas.read_csv(folder + "MA_plant_data.txt", sep='\t', index_col=0, header=0)
    MA_input_data = pandas.read_csv(folder + "MA_input_data.txt", sep='\t', index_col=0, header=0)
    CPWL_MA_plant_data = pandas.read_csv(folder + "CPWL_MA_plant_data.txt", sep='\t', index_col=0, header=0)
    CPWL_MA_input_data = pandas.read_csv(folder + "CPWL_MA_input_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_Subgrad_plant_data = pandas.read_csv(folder + "QCPWL_Subgrad_plant_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_Subgrad_input_data = pandas.read_csv(folder + "QCPWL_Subgrad_input_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_MINLP_plant_data = pandas.read_csv(folder + "QCPWL_MINLP_plant_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_MINLP_input_data = pandas.read_csv(folder + "QCPWL_MINLP_input_data.txt", sep='\t', index_col=0, header=0)

    max_iter=20
    # #----------Optimal-------------
    ori_plant_cost = 190.98 * np.ones(max_iter+1)
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
    Optimal, = plt.plot(range(max_iter+1), ori_plant_cost, linewidth=global_linewidth, label='Optimal', color='gray',
                        linestyle='--')
    MA, = plt.plot(range(max_iter), MA_plant_data.loc[1:(max_iter+1),"cost"].multiply(-1), linewidth=global_linewidth, label='Nonlinear model', linestyle='dashdot',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_plant_data.loc[1:(max_iter+1),"cost"].multiply(-1), linewidth=global_linewidth, label='DCA+CPWL', linestyle='dotted',
                    marker='o', markersize=global_point_size)
    QCPWL_Subgrad, = plt.plot(range(max_iter), QCPWL_Subgrad_plant_data.loc[1:(max_iter+1), "cost"].multiply(-1),
                              linewidth=global_linewidth,
                              label='DCA+QPWL',
                              linestyle='-', marker='o', markersize=global_point_size)
    QCPWL_MINLP, = plt.plot(range(max_iter), QCPWL_MINLP_plant_data.loc[1:(max_iter+1), "cost"].multiply(-1),
                            linewidth=global_linewidth,
                            label='MINLP+QPWL',
                            linestyle='dashed', marker='o', markersize=global_point_size/2)

    plt.legend(handles=[MA, CPWL_MA, QCPWL_Subgrad, QCPWL_MINLP, Optimal], prop=font1)
    plt.legend(loc='lower right',
               fancybox=False,
               edgecolor='k',
               fontsize=global_point_size,
               shadow=False,
               facecolor='w',
               framealpha=1.0,
               prop={'family': 'Times New Roman', 'size': global_tick_size})
    plt.ylabel('(a) Plant cost', font2)
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
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='Nonlinear model', linestyle='dashdot',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='DCA+CPWL', linestyle='dotted',
                    marker='o', markersize=global_point_size)
    QCPWL_Subgrad, = plt.plot(range(max_iter), QCPWL_Subgrad_input_data.loc[range(max_iter), "Fb"],
                              linewidth=global_linewidth,
                              label='DCA+QPWL',
                              linestyle='-', marker='o', markersize=global_point_size)
    QCPWL_MINLP, = plt.plot(range(max_iter), QCPWL_MINLP_input_data.loc[range(max_iter), "Fb"],
                            linewidth=global_linewidth,
                            label='MINLP+QPWL',
                            linestyle='dashed', marker='o', markersize=global_point_size/2)

    plt.ylabel('(b) $F_B$ (kg/s)', font2)
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
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='Nonlinear model', linestyle='dashdot',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='DCA+CPWL', linestyle='dashed',
                    marker='o', markersize=global_point_size)
    QCPWL_Subgrad, = plt.plot(range(max_iter), QCPWL_Subgrad_input_data.loc[range(max_iter), "Tr"], linewidth=global_linewidth,
                         label='DCA+QPWL',
                         linestyle='-', marker='o', markersize=global_point_size)
    QCPWL_MINLP, = plt.plot(range(max_iter), QCPWL_MINLP_input_data.loc[range(max_iter), "Tr"], linewidth=global_linewidth,
                         label='MINLP+QPWL',
                         linestyle='dashed', marker='o', markersize=global_point_size/2)

    # plt.legend(handles=[OptimalTR, MA, CPWL_MA, PA, QCPWL_MA, ISOQCPWL_MA], prop=font1)
    plt.xlabel('RTO Iterations', font2)
    plt.ylabel('(c) $T_R$ (℃)', font2)
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

def HPC_plant_optimum(iter):
    if iter % 30 >= 20:
        return {'TA_Flow': 271.9211899185857, 'MA_Flow': 1687.6282683022116, 'GAN_Flow': 750, 'LN_Flow': 100,
         'Purity_GAN': 0.9989999897973798, 'T_GAN': -178.05452423466437, 'T_30thTray': -177.26340112601463,
         'WN_Flow': 1099.549456249657, 'Drain_Flow': 9.999998719778556, 'obj': 3992.5635796531133,
         'feed_con1': -0.07596624389862891, 'feed_con2': -0.5118985933139635, 'purity_con': 1.0202620193133782e-08,
         'drain_con': 1.2802214435225778e-06}

    elif iter % 30 >= 10:
        return {'TA_Flow': 277.6176452038008, 'MA_Flow': 1791.9833231204288, 'GAN_Flow': 800, 'LN_Flow': 100,
         'Purity_GAN': 0.9989999898170692, 'T_GAN': -178.05452488930501, 'T_30thTray': -177.26131328659508,
         'WN_Flow': 1159.6009664557073, 'Drain_Flow': 9.999998778993419, 'obj': 4194.195296505531,
         'feed_con1': -0.07065754818279839, 'feed_con2': -0.5502227420732334, 'purity_con': 1.0182930831881265e-08,
         'drain_con': 1.2210065811757431e-06}
    else:
        return {'TA_Flow': 275.3390329075402, 'MA_Flow': 1750.240047108621, 'GAN_Flow': 780, 'LN_Flow': 100,
                'Purity_GAN': 0.9989999898095535, 'T_GAN': -178.05452463566908, 'T_30thTray': -177.2621236594108,
                'WN_Flow': 1135.579078107841, 'Drain_Flow': 9.999998756038286, 'obj': 4113.540342439814,
                'feed_con1': -0.0727811247137334, 'feed_con2': -0.5348925990450646, 'purity_con': 1.0190446486646465e-08,
                'drain_con': 1.2439617140813652e-06}

def draw_compare_HPC3(pic_name):
    max_iter_no=30
    pic_constant = 0.39370
    global_font_size = 15
    global_tick_size = 15
    global_linewidth = 1.5
    global_point_size = 4
    font_factor = np.sqrt(1 / pic_constant)
    fig = plt.figure(figsize=(13 * pic_constant, 30 * pic_constant))
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font3 = {'family': 'Times New Roman',
             'weight': 'bold',
             'size': global_font_size,
             }
    font_legend = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 16/font_factor,
             }
    font_title = {'family': 'Times New Roman',
                  'size': global_font_size,
                  'weight': 'bold'
                  }
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # algorighms = ['QCPWL_MA', 'MA','QCPWL_Subgrad','QCPWL_Subgrad_MINLP']
    algorighms = ['MA','QCPWL_Subgrad','QCPWL_Subgrad_MINLP']
    algo_name_dict = {
        'MA':"Nonlinear model",
        'QCPWL_Subgrad':'DCA+QPWL',
        'QCPWL_Subgrad_MINLP':'MINLP+QPWL',
    }
    algo_color = {
        'MA':"blue",
        'QCPWL_Subgrad':"red",
        'QCPWL_Subgrad_MINLP':"black"
    }
    output_measurements = ['obj','Purity_GAN','Drain_Flow']
    inputs = ['TA_Flow', 'MA_Flow']
    mv_label_dict = {
        'TA_Flow': '$F_{TA}$ (kmol/h)',
        'MA_Flow': '$F_{MA}$ (kmol/h)',
    }
    y_label_dict = {
        'obj':'$\phi_p$',
        'Purity_GAN':'$c_{GAN,N_2}$',
        'Drain_Flow':'$F_{Drain}$ (kmol/h)',
    }
    linestyle_dict={
        'MA': "dashed",
        'QCPWL_Subgrad': "dashdot",
        'QCPWL_Subgrad_MINLP': "dotted"
    }
    marker_size={
        'MA': 4,
        'QCPWL_Subgrad': 4,
        'QCPWL_Subgrad_MINLP': 2
    }
    for op_index, op in enumerate(output_measurements):
        plt.subplot(5, 1, op_index + 1)
        time=np.zeros((max_iter_no,))
        for iter in range(1,max_iter_no):
            time[iter]=iter+1
        optimal=np.zeros((max_iter_no,))
        for index,iter in enumerate(time):
            optimal[index]=HPC_plant_optimum(iter)[op]
        optimal[-1]=HPC_plant_optimum(time[-2])[op]
        plt.plot(time, optimal, linewidth=global_linewidth/2,
                 label='Optimal',
                 color='darkgreen',
                 linestyle='-')
        for algo_index, algo in enumerate(algorighms):
            model_data = pandas.read_csv("data/hpc3d/" + algo + '_model_data' + ".txt", sep='\t', index_col=0, header=0)
            plant_data = pandas.read_csv("data/hpc3d/" + algo + '_plant_data' + ".txt", sep='\t', index_col=0, header=0)
            # a stupid remedy for the bug of the plant data
            if op_index == 1:
                y=list(plant_data.iloc[0:(max_iter_no-1), :].loc[:, op])
                y.insert(0,0.9999)
                time = model_data.index
                plt.plot(time[0:max_iter_no], y,
                         linewidth=global_linewidth, label=algo_name_dict[algo],
                         color=algo_color[algo],
                         linestyle=linestyle_dict[algo], marker='o', markersize=marker_size[algo])
            else:
                plant_data.iloc[max_iter_no-1] = plant_data.iloc[max_iter_no - 2]
                time = model_data.index
                plt.plot(time[0:max_iter_no], plant_data.iloc[0:(max_iter_no),:].loc[:,op], linewidth=global_linewidth, label=algo_name_dict[algo],
                                           color=algo_color[algo],
                                           linestyle=linestyle_dict[algo], marker='o', markersize=marker_size[algo])
        if op_index == 0:
            plt.legend(prop=font_legend)
        plt.ylabel(y_label_dict[op], font2)
        plt.xticks([0,10,20,30])
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.tick_params(labelsize=global_tick_size)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        plt.tick_params(top= ['on', 'true'], right= ['on', 'true'], which='both')
        plt.grid(ls='--',axis='y')

    for op_index, op in enumerate(inputs):
        plt.subplot(5, 1, op_index + 4)
        time = np.zeros((max_iter_no,))
        for iter in range(1, max_iter_no):
            time[iter] = iter + 1
        optimal = np.zeros((max_iter_no,))
        for index, iter in enumerate(time):
            optimal[index] = HPC_plant_optimum(iter)[op]
        optimal[-1] = HPC_plant_optimum(time[-2])[op]
        plt.plot(time, optimal, linewidth=global_linewidth / 2,
                 label='Optimal',
                 color='darkgreen',
                 linestyle='-')
        for algo_index, algo in enumerate(algorighms):
            model_data = pandas.read_csv("data/hpc3d/" + algo + '_model_data' + ".txt", sep='\t', index_col=0, header=0)
            input_data = pandas.read_csv("data/hpc3d/" + algo + '_input_data' + ".txt", sep='\t', index_col=0, header=0)
            input_data.iloc[max_iter_no] = input_data.iloc[max_iter_no-1]
            time = model_data.index
            plt.plot(time[0:max_iter_no], input_data.iloc[0:(max_iter_no),:].loc[:,op], linewidth=global_linewidth, label=algo,
                                       color=algo_color[algo],
                                           linestyle=linestyle_dict[algo], marker='o', markersize=marker_size[algo])

        plt.ylabel(mv_label_dict[op], font2)
        if op_index == 1:
            plt.xlabel("RTO iteration", font2)
        plt.xticks([0,10,20,30])
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.tick_params(labelsize=global_tick_size)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        plt.tick_params(top= ['on', 'true'], right= ['on', 'true'], which='both')
        plt.grid(ls='--',axis='y')

    plt.tight_layout(pad=0.5,w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_name,dpi=600)
    plt.close()


def draw_compare_HPC3_CPCC(pic_name):
    max_iter_no=30
    pic_constant = 0.39370
    global_font_size = 15
    global_tick_size = 15
    global_linewidth = 1.5
    global_point_size = 4
    font_factor = np.sqrt(1 / pic_constant)
    fig = plt.figure(figsize=(16 * pic_constant, 16 * pic_constant))
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font3 = {'family': 'Times New Roman',
             'weight': 'bold',
             'size': global_font_size,
             }
    font_legend = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 20/font_factor,
             }
    font_title = {'family': 'Times New Roman',
                  'size': global_font_size,
                  'weight': 'bold'
                  }
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # algorighms = ['QCPWL_MA', 'MA','QCPWL_Subgrad','QCPWL_Subgrad_MINLP']
    algorighms = ['MA','QCPWL_Subgrad']
    algo_name_dict = {
        'MA':"Modifier adaptation",
        'QCPWL_Subgrad':'Proposed',
    }
    algo_color = {
        'MA':"blue",
        'QCPWL_Subgrad':"red",
    }
    output_measurements = ['obj','Purity_GAN','Drain_Flow']
    inputs = ['TA_Flow', 'MA_Flow']
    mv_label_dict = {
        'TA_Flow': 'Turbine air flowrate (kmol/h)',
        'MA_Flow': 'Main air flowrate (kmol/h)',
    }
    y_label_dict = {
        'obj':'Cost',
        'Purity_GAN':'Purity of GAN',
        'Drain_Flow':'Drain flowrate (kmol/h)',
    }
    linestyle_dict={
        'MA': "dashed",
        'QCPWL_Subgrad': "dashdot",
    }
    marker_size={
        'MA': 2,
        'QCPWL_Subgrad': 2,
    }
    for op_index, op in enumerate(output_measurements):
        plt.subplot(3, 1, op_index + 1)
        time=np.zeros((max_iter_no,))
        for iter in range(1,max_iter_no):
            time[iter]=iter+1
        optimal=np.zeros((max_iter_no,))
        for index,iter in enumerate(time):
            optimal[index]=HPC_plant_optimum(iter)[op]
        optimal[-1]=HPC_plant_optimum(time[-2])[op]
        plt.plot(time, optimal, linewidth=global_linewidth/2,
                 label='Optimal',
                 color='black',
                 linestyle='-')
        for algo_index, algo in enumerate(algorighms):
            model_data = pandas.read_csv("data/hpc3d/" + algo + '_model_data' + ".txt", sep='\t', index_col=0, header=0)
            plant_data = pandas.read_csv("data/hpc3d/" + algo + '_plant_data' + ".txt", sep='\t', index_col=0, header=0)
            plant_data.iloc[max_iter_no-1] = plant_data.iloc[max_iter_no - 2]
            time = model_data.index
            plt.plot(time[0:max_iter_no], plant_data.iloc[0:(max_iter_no),:].loc[:,op], linewidth=global_linewidth, label=algo_name_dict[algo],
                                       color=algo_color[algo],
                                       linestyle=linestyle_dict[algo], marker='o', markersize=marker_size[algo])
        if op_index == 0:
            plt.legend(prop=font_legend)
        plt.ylabel(y_label_dict[op], font2)
        plt.xticks([0,10,20,30])
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.tick_params(labelsize=global_tick_size)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        plt.tick_params(top= ['on', 'true'], right= ['on', 'true'], which='both')
        plt.grid(ls='--',axis='y')

    plt.subplot(3, 1, 3)
    plt.xlabel("RTO iteration", font2)

    plt.tight_layout(pad=0.5,w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_name,dpi=600)
    plt.close()

if __name__ == "__main__":
    # draw_compare_HPC3("pic/paper/compare_HPC3.png")
    # draw_compare_HPC3_CPCC("pic/paper/compare_HPC3_CPCC.png")
    # draw_compare_unconstrained_wo("data/uncon-wo/", "pic/paper/compare_unconstrained_wo.png")
    draw_compare_unconstrained_wo_adchem("data/uncon-wo/",
                                  "pic/paper/compare_unconstrained_wo_adchem.png")