import pandas
import do_test_lib
import numpy

do_test_lib.test_compare_modifiers()

folder="data/comp_modifiers/"
MA_model_data = pandas.read_csv(folder + "MA_model_data.txt", sep='\t', index_col=0, header=0)
GPE10_model_data = pandas.read_csv(folder + "GPE10_model_data.txt", sep='\t', index_col=0, header=0)
GPE1_model_data = pandas.read_csv(folder + "GPE1_model_data.txt", sep='\t', index_col=0, header=0)
ISOPE_model_data = pandas.read_csv(folder + "ISOPE_model_data.txt", sep='\t', index_col=0, header=0)

with open("data/comp_modifiers/synthesis.txt", 'w') as fp:
    fp.write("method\tepsilon_data\tlambda_Fb_data\tlambda_Tr_data\n")
    fp.write("MA\t%.2f\t%.2f\t%.2f\n"%(numpy.log10(MA_model_data.loc[:,'cost_eps'].abs().mean()),\
                                     numpy.log10(MA_model_data.loc[:,'Fb_cost_lam'].abs().mean()),\
                                     numpy.log10(MA_model_data.loc[:,'Tr_cost_lam'].abs().mean())))
    fp.write("PE+MA\t%.2f\t%.2f\t%.2f\n" % (numpy.log10(ISOPE_model_data.loc[:, 'cost_eps'].abs().mean()), \
                                         numpy.log10(ISOPE_model_data.loc[:, 'Fb_cost_lam'].abs().mean()), \
                                         numpy.log10(ISOPE_model_data.loc[:, 'Tr_cost_lam'].abs().mean())))
    fp.write("GPE(n=10)\t%.2f\t%.2f\t%.2f\n" % (numpy.log10(GPE10_model_data.loc[:, 'cost_eps'].abs().mean()), \
                                         numpy.log10(GPE10_model_data.loc[:, 'Fb_cost_lam'].abs().mean()), \
                                         numpy.log10(GPE10_model_data.loc[:, 'Tr_cost_lam'].abs().mean())))
    fp.write("GPE(n=1)\t%.2f\t%.2f\t%.2f\n" % (numpy.log10(GPE1_model_data.loc[:, 'cost_eps'].abs().mean()), \
                                         numpy.log10(GPE1_model_data.loc[:, 'Fb_cost_lam'].abs().mean()), \
                                         numpy.log10(GPE1_model_data.loc[:, 'Tr_cost_lam'].abs().mean())))