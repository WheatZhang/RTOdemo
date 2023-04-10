import pandas
import numpy

def calculate_indices(data, optimal):
    devi = 0
    for i in range(30):
        devi += numpy.abs(data.loc[i+1]-optimal)

    devi10 = 0
    for i in range(10):
        devi10 += numpy.abs(data.loc[i+1]-optimal)

    devi25 = 0
    for i in range(5):
        devi25 += numpy.abs(data.loc[i + 26] - optimal)

    aver5 = 0
    for i in range(5):
        aver5 += data.loc[i+6]-optimal
    aver5 /= 5

    aver25 = 0
    for i in range(5):
        aver25 += data.loc[i + 26] - optimal
    aver25 /= 5

    for i in range(30):
        if numpy.abs(data.loc[i+1]-optimal) < max(abs(optimal),0.1)/10:
            break
    speed = i

    monoto = 0
    for i in range(29):
        monoto += numpy.abs(data.loc[i+2]-data.loc[i+1])
    monoto -= numpy.abs(data.loc[1]-(aver25+optimal))
    monoto = max(monoto, 0)

    return devi, devi10, aver5, devi25, aver25, speed, monoto


def compo_step2_test_report():
    optimal = {
        "cost" : 1.453659e-01,
        "con" : 0,
        "u1" : 3.684054e-01,
        "u2" : -3.929466e-01,
    }

    algo_report = pandas.DataFrame()
    batch_prefixes = ["N_", "O1_", "O2_", "C_", "OC_"]
    batch_names = [
        "U0", "U1", "U2", "U3", "U4", "U5", "U6", "U7",
        "R0", "R1", "R2", "R3", "R4",
        "P0", "P1", "P2", "P3", "P4", "P5",
        "PF0", "PF1", "PF2", "PF3", "PF4", "PF5",
        "PIF0", "PIF1", "PIF2", "PIF3", "PIF4", "PIF5",
        "X0", "X1", "X2", "X3", "X4",
        "XF0", "XF1", "XF2", "XF3", "XF4",
        "XIF0", "XIF1", "XIF2", "XIF3", "XIF4",
        "XR0", "XR1", "XR2", "XR3", "XR4",
        "XRF0", "XRF1", "XRF2", "XRF3", "XRF4",
        "XRIF0", "XRIF1", "XRIF2", "XRIF3", "XRIF4",
    ]
    for batch_prefix in batch_prefixes:
        for batch_name in batch_names:
            data_prefix = batch_prefix + batch_name
            compo_step_plant_data = pandas.read_csv("data/batch12/"+data_prefix+"/CompoStep_TR_MA_plant_data.txt", \
                                              index_col=0, header=0, sep='\t')
            penalty_plant_data = pandas.read_csv("data/batch12/"+data_prefix+"/Penalty_TR_MA_plant_data.txt", \
                                           index_col=0, header=0, sep='\t')
            compo_step_model_data = pandas.read_csv("data/batch12/"+data_prefix+"/CompoStep_TR_MA_model_data.txt", \
                                              index_col=0, header=0, sep='\t')
            penalty_model_data = pandas.read_csv("data/batch12/"+data_prefix+"/Penalty_TR_MA_model_data.txt", \
                                           index_col=0, header=0, sep='\t')
            compo_step_input_data = pandas.read_csv("data/batch12/" + data_prefix + "/CompoStep_TR_MA_input_data.txt", \
                                                    index_col=0, header=0, sep='\t')
            penalty_input_data = pandas.read_csv("data/batch12/" + data_prefix + "/Penalty_TR_MA_input_data.txt", \
                                                    index_col=0, header=0, sep='\t')
            test_algo="compo"
            this_test_data = {}
            for test_var in ["cost", "con"]:
                devi, devi10, aver5, devi25, aver25, speed, monoto = \
                    calculate_indices(compo_step_plant_data.loc[:,test_var], optimal[test_var])
                this_test_data[test_algo+"_"+test_var+"_"+"devi"] = devi
                this_test_data[test_algo+"_"+test_var+"_"+"devi10"] = devi10
                this_test_data[test_algo+"_"+test_var+"_"+"aver5"] = aver5
                this_test_data[test_algo+"_"+test_var+"_"+"devi25"] = devi25
                this_test_data[test_algo+"_"+test_var+"_"+"aver25"] = aver25
                this_test_data[test_algo+"_"+test_var+"_"+"speed"] = speed
                this_test_data[test_algo+"_"+test_var+"_"+"monoto"] = monoto

            for test_var in ["u1", "u2"]:
                devi, devi10, aver5, devi25, aver25, speed, monoto = \
                    calculate_indices(compo_step_input_data.loc[:,test_var], optimal[test_var])
                this_test_data[test_algo+"_"+test_var+"_"+"devi"] = devi
                this_test_data[test_algo+"_"+test_var+"_"+"devi10"] = devi10
                this_test_data[test_algo+"_"+test_var+"_"+"aver5"] = aver5
                this_test_data[test_algo+"_"+test_var+"_"+"devi25"] = devi25
                this_test_data[test_algo+"_"+test_var+"_"+"aver25"] = aver25
                this_test_data[test_algo+"_"+test_var+"_"+"speed"] = speed
                this_test_data[test_algo+"_"+test_var+"_"+"monoto"] = monoto

            this_test_data[test_algo+"_last_r"] = compo_step_model_data.loc[26:31,"tr"].mean()

            test_algo = "penalty"
            for test_var in ["cost", "con"]:
                devi, devi10, aver5, devi25, aver25, speed, monoto = \
                    calculate_indices(penalty_plant_data.loc[:, test_var], optimal[test_var])
                this_test_data[test_algo + "_" + test_var + "_" + "devi"] = devi
                this_test_data[test_algo + "_" + test_var + "_" + "devi10"] = devi10
                this_test_data[test_algo + "_" + test_var + "_" + "aver5"] = aver5
                this_test_data[test_algo + "_" + test_var + "_" + "devi25"] = devi25
                this_test_data[test_algo + "_" + test_var + "_" + "aver25"] = aver25
                this_test_data[test_algo + "_" + test_var + "_" + "speed"] = speed
                this_test_data[test_algo + "_" + test_var + "_" + "monoto"] = monoto

            for test_var in ["u1", "u2"]:
                devi, devi10, aver5, devi25, aver25, speed, monoto = \
                    calculate_indices(penalty_input_data.loc[:, test_var], optimal[test_var])
                this_test_data[test_algo + "_" + test_var + "_" + "devi"] = devi
                this_test_data[test_algo + "_" + test_var + "_" + "devi10"] = devi10
                this_test_data[test_algo + "_" + test_var + "_" + "aver5"] = aver5
                this_test_data[test_algo + "_" + test_var + "_" + "devi25"] = devi25
                this_test_data[test_algo + "_" + test_var + "_" + "aver25"] = aver25
                this_test_data[test_algo + "_" + test_var + "_" + "speed"] = speed
                this_test_data[test_algo + "_" + test_var + "_" + "monoto"] = monoto

            this_test_data[test_algo+"_last_r"] = penalty_model_data.loc[26:31, "tr"].mean()

            algo_report = algo_report.append(pandas.DataFrame(data=this_test_data, index=[data_prefix]))
    algo_report.to_csv("data/report/compo_step2_report.txt", sep='\t', header=True, index=True)


def synthesize_compo_step2_report():
    original_report = pandas.read_csv("data/report/compo_step2_report.txt", sep='\t', header=0, index_col=0)
    synthesized_report = pandas.DataFrame(index=original_report.index, \
                                          columns=["Converged", "Monotone", "Algo12Diff", "LastTR",\
                                                   "ObjSpeed", "ConSpeed"])

    for index, row in original_report.iterrows():
        algo_1_converged = not(row.loc["compo_cost_speed"] >= 25 or row.loc["compo_con_speed"] >= 25 or \
                               row.loc["compo_cost_devi25"] > 1e-1 or row.loc["compo_con_devi25"] > 1e-1)
        algo_2_converged = not(row.loc["penalty_cost_speed"] >= 25 or row.loc["penalty_con_speed"] >= 25 or\
                               row.loc["penalty_cost_devi25"] > 1e-1 or row.loc["penalty_con_devi25"] > 1e-1)
        if algo_1_converged and algo_2_converged:
            synthesized_report.loc[index, "Converged"] = "Both"
        elif algo_1_converged and not algo_2_converged:
            synthesized_report.loc[index, "Converged"] = "Alg1"
        elif not algo_1_converged and algo_2_converged:
            synthesized_report.loc[index, "Converged"] = "Alg2"
        else:
            synthesized_report.loc[index, "Converged"] = "Fail"

        algo_1_monotone = row.loc["compo_cost_monoto"] < 1 and row.loc["compo_con_monoto"] < 1
        algo_2_monotone = row.loc["penalty_cost_monoto"] < 1 and row.loc["penalty_con_monoto"] < 1
        if algo_1_monotone and algo_2_monotone:
            synthesized_report.loc[index, "Monotone"] = "Mono"
        elif algo_1_monotone and not algo_2_monotone:
            synthesized_report.loc[index, "Monotone"] = "Alg1"
        elif not algo_1_monotone and algo_2_monotone:
            synthesized_report.loc[index, "Monotone"] = "Alg2"
        else:
            synthesized_report.loc[index, "Monotone"] = "Osc"

        algo_12_diff = max(abs(row.loc["compo_cost_devi"]-row.loc["penalty_cost_devi"]),\
                           abs(row.loc["compo_con_devi"]-row.loc["penalty_con_devi"]))
        if algo_12_diff < 1e-1:
            synthesized_report.loc[index, "Algo12Diff"] = "Identical"
        elif algo_12_diff < 1:
            synthesized_report.loc[index, "Algo12Diff"] = "Likewise"
        elif algo_12_diff < 5:
            synthesized_report.loc[index, "Algo12Diff"] = "Different"
        else:
            synthesized_report.loc[index, "Algo12Diff"] = "VeryDiff"

        if row.loc["compo_last_r"] < 0.1 and row.loc["penalty_last_r"] < 0.1:
            synthesized_report.loc[index, "LastTR"] = "Small"
        elif row.loc["compo_last_r"] < 0.1 and row.loc["penalty_last_r"] < 1:
            synthesized_report.loc[index, "LastTR"] = "1S2M"
        elif row.loc["compo_last_r"] < 0.1:
            synthesized_report.loc[index, "LastTR"] = "1S2L"
        elif row.loc["compo_last_r"] < 1 and row.loc["penalty_last_r"] < 0.1:
            synthesized_report.loc[index, "LastTR"] = "1M2S"
        elif row.loc["compo_last_r"] < 1 and row.loc["penalty_last_r"] < 1:
            synthesized_report.loc[index, "LastTR"] = "Medium"
        elif row.loc["compo_last_r"] < 1:
            synthesized_report.loc[index, "LastTR"] = "1M2L"
        elif row.loc["penalty_last_r"] < 0.1:
            synthesized_report.loc[index, "LastTR"] = "1L2S"
        elif row.loc["penalty_last_r"] < 1:
            synthesized_report.loc[index, "LastTR"] = "1L2M"
        else:
            synthesized_report.loc[index, "LastTR"] = "Large"

        speed_diff = row.loc["compo_cost_devi10"]-row.loc["penalty_cost_devi10"]
        if speed_diff < -1e-2:
            synthesized_report.loc[index, "ObjSpeed"] = "Alg1Faster"
        elif algo_12_diff > 1e-2:
            synthesized_report.loc[index, "ObjSpeed"] = "Alg2Faster"
        else:
            synthesized_report.loc[index, "ObjSpeed"] = "Similar"

        speed_diff = row.loc["compo_con_devi10"] - row.loc["penalty_con_devi10"]
        if speed_diff < -1e-2:
            synthesized_report.loc[index, "ConSpeed"] = "Alg1Faster"
        elif algo_12_diff > 1e-2:
            synthesized_report.loc[index, "ConSpeed"] = "Alg2Faster"
        else:
            synthesized_report.loc[index, "ConSpeed"] = "Similar"

    synthesized_report.index.name = "Case"
    synthesized_report.to_csv("data/report/compo_step2_synthesized_report.txt", sep='\t', header=True, index=True)

if __name__ == "__main__":
    compo_step2_test_report()
    synthesize_compo_step2_report()