import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt


class DAneal:
    def __init__(self, k, Beta, BetaMax, alpha):
        self.k = k
        self.Beta = Beta
        self.BetaMax = BetaMax
        self.alpha = alpha

    def Anneal(self, data):
        #step 1
        self.data = data
        x = data[:,0]
        y = data[:,1]
        N = x.shape[0];

        #Initialization
        centers = np.ones((self.k, 2))
        Pcenters = np.ones((self.k, N))

        #PERTURB = 0.001*ones((self.k,2))
        PERTURB = 0.001
        STOP = 1e-4

        for i in range(self.k):
            centers[i,0] = np.sum(x)/N
            centers[i,1] = np.sum(y)/N

        print("The initial centers are")
        print(centers)

        while self.Beta <= self.BetaMax:
            #update the Pcenters for all the points
            # itr = 1  #Put the F old here
            # while itr < 20:
            for point in range(N):
                for center in range(self.k):
                    curr_dist = (x[point]-centers[center,0])**2 + (y[point]-centers[center,1])**2
                    num = np.exp(-self.Beta*curr_dist)

                    denum = 0.0
                    for i in range(self.k):
                        dist = (x[point]-centers[i,0])**2 + (y[point]-centers[i,1])**2
                        denum += np.exp(-self.Beta*dist)
                        #print("the denumorator is")
                        #print(denum)
                    Pcenters[center,point] = num / denum


            #update all centers
            for j in range(self.k):
                Xcoord = 0.0
                Ycoord = 0.0
                denum = 0.0
                for i in range(N):
                    Xcoord += Pcenters[j,i] * x[i]
                    Ycoord += Pcenters[j,i] * y[i]
                    denum += Pcenters[j,i]
                centers[j,0] = Xcoord/denum + PERTURB*random.random()
                centers[j,1] = Ycoord/denum + PERTURB*random.random()
                # itr += 1

            self.Beta *= self.alpha
        Pcenters = np.around(Pcenters)
        return (centers, Pcenters)



#Import the data set
rawData = pd.read_csv("ipl.csv")
rawData = rawData[['one', 'two']]
data_arr = rawData.values

#Create the object and pass the data in
DA = DAneal(3, 0.001, 1000, 1.1)
centers, Pcenters = DA.Anneal(data_arr)
print(Pcenters)

colors = 10*["r", "g", "c", "b", "k"]
for centroid in range(DA.k):
    plt.scatter(centers[centroid][0], centers[centroid][1], s = 130, marker = "x")

for cluster in range(DA.k):
    color = colors[cluster]
    valid_X = data_arr[:,0][np.where(Pcenters[cluster] == 1)]
    valid_Y = data_arr[:,1][np.where(Pcenters[cluster] == 1)]
    plt.scatter(valid_X, valid_Y, color = color, s = 30)
