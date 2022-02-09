import pandas
import matplotlib.pyplot as plt

def draw_overview_pic():
    files = ['model_data', 'plant_data']

    for filename in files:
        temp = pandas.read_csv("result/"+filename+".txt", sep='\t', index_col=0, header=0)
        time = temp.index
        vars = temp.columns.to_list()
        for v in vars:
            fig = plt.figure()
            temp = pandas.read_csv("result/" + filename + ".txt", sep='\t', index_col=0, header=0)
            plt.plot(time, temp.loc[:, v])
            plt.title(v)
            plt.savefig("pic/"+v+"_"+filename+".png")
            plt.close(fig)


if __name__ == "__main__":
    draw_overview_pic()