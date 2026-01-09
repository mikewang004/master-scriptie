import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 

class make_plot():
    def __init__(self):
        return

    def scatter_plot(xdata: list, ydata: list, labels: list = None, xlabel = None, ylabel = None, title = None, save_string: str = None, show_plot: bool = False, marker = "o"):
        """xdata, ydata lists of numpy arrays"""

        #print(xdata, ydata)
        if isinstance(xdata, list):
            if len(xdata) == len(ydata):
                for i in range(len(xdata)):
                    if labels != None:
                        print(xdata[i], ydata[i], labels[i])
                        #print(xdata[i].shape, ydata[i].shape, labels[i].shape)
                        plt.scatter(xdata[i],ydata[i],label = labels[i], marker = marker)
                        plt.legend()
                    else:
                        plt.scatter(xdata[i],ydata[i], marker = marker)
            else:
                raise Exception("x,y need have same len")
        else:
            if isinstance(xdata, np.ndarray):
                if isinstance(ydata, np.ndarray):
                    if isinstance(labels, str):
                        plt.scatter(xdata, ydata, label = labels, marker = marker)
                    else:
                        plt.scatter(xdata, ydata)
                else: 
                    raise Exception("x,y need have same len")


        if isinstance(xlabel, str):
            plt.xlabel(xlabel)
        if isinstance(ylabel, str):
            plt.ylabel(ylabel)
        if isinstance(title, str):
            plt.title(title)

        if isinstance(save_string, str):
            plt.savefig(save_string)
        if show_plot == True:
            plt.show()


    