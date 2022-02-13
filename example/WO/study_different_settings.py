from do_test_lib import batch_test_GPE, batch_test_all_algo, different_para_set_test
from draw_pic_lib import draw_para_set_comp_pic,draw_overall_pic,draw_GPE_different_weight
import os

for n in [0, 0.5, 1, 2]:
    print("n=%f" % n)
    batch_test_all_algo(delta_u_Fb=0.2, filtering_factor=0.5, noise_level_coeff=n)

for K in [0.25, 0.5, 0.75]:
    print("K=%f" % K)
    batch_test_all_algo(delta_u_Fb=0.2, filtering_factor=K, noise_level_coeff=1)

for delta_u in [0.1, 0.2, 0.4]:
    print("delta_u=%f" % delta_u)
    batch_test_all_algo(delta_u_Fb=delta_u, filtering_factor=0.5, noise_level_coeff=1)

for K in [0.25, 0.5, 0.75]:
    print("K=%f" % K)
    batch_test_all_algo(delta_u_Fb=0.2, filtering_factor=K, noise_level_coeff=0)

for delta_u in [0.1, 0.2, 0.4]:
    print("delta_u=%f" % delta_u)
    batch_test_all_algo(delta_u_Fb=delta_u, filtering_factor=0.5, noise_level_coeff=0)

for factor_n in [0.1, 1, 10, 100, 1000]:
    print("factor_n=%f" % factor_n)
    batch_test_GPE(factor_n=factor_n)

different_para_set_test(ka_relative_uncertainty=1, kb_relative_uncertainty=1)
different_para_set_test(ka_relative_uncertainty=0.1, kb_relative_uncertainty=0.1)
different_para_set_test(ka_relative_uncertainty=0.01, kb_relative_uncertainty=0.01)
different_para_set_test(ka_relative_uncertainty=0.1, kb_relative_uncertainty=-1)
different_para_set_test(ka_relative_uncertainty=-1, kb_relative_uncertainty=0.1)

for fname in os.listdir("data/", ):
    if os.path.isdir("data/" + fname):
        if fname.startswith("batch_all"):
            print(fname)
            draw_overall_pic("data/" + fname + "/", "pic/overall_pic_" + fname + ".png")

draw_GPE_different_weight()

for fname in os.listdir("data/", ):
    if os.path.isdir("data/" + fname):
        if fname.startswith("batch_para"):
            print(fname)
            draw_para_set_comp_pic("data/" + fname + "/", "pic/overall_pic_" + fname + ".png")

