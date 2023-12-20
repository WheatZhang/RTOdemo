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


def draw_x2():
    # plot y=x^2
    # x in [-2.5, 2.5]
    x = [i/10 for i in range(-25, 26)]
    y = [i**2 for i in x]
    plt.plot(x, y)
    plt.show()
if __name__ == "__main__":
    # draw_overview_pic()
    draw_x2()