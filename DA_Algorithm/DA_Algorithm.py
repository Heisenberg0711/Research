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

        while self.Beta <= self.BetaMax:
            #update the Pcenters for all the points
             itr = 1  #Put the F old here+
             while itr < 20:
                 #Updating the Pcenters
                    # for point in range(N):
                    #     for center in range(self.k):
                    #         curr_dist = (x[point]-centers[center,0])**2 + (y[point]-centers[center,1])**2
                    #         num = np.exp(-self.Beta*curr_dist)
                    #
                    #         denum = 0.0
                    #         for i in range(self.k):
                    #             dist = (x[point]-centers[i,0])**2 + (y[point]-centers[i,1])**2
                    #             denum += np.exp(-self.Beta*dist)
                    #
                    #         Pcenters[center,point] = num / denum

                    Pcenters_new = np.zeros([self.k, N])
                    for centroid in range(self.k):
                        dist = data - np.tile(centers[centroid,:], (N, 1))
                        dist = (dist[:,0]**2 + dist[:,1]**2).reshape(N, 1)
                        dist = np.exp(-self.Beta * dist).reshape(1,N)
                        Pcenters_new[centroid,:] = dist

                    denum = np.sum(Pcenters)
                    for i in range(self.k):
                        Pcenters[i,:] = Pcenters_new[i,:] / denum


                    #update all centers
                    # for j in range(self.k):
                    #     Xcoord = 0.0
                    #     Ycoord = 0.0
                    #     denum = 0.0
                    #     for i in range(N):
                    #         Xcoord += Pcenters[j,i] * x[i]
                    #         Ycoord += Pcenters[j,i] * y[i]
                    #         denum += Pcenters[j,i]
                    #     centers[j,0] = Xcoord/denum + PERTURB*random.random()
                    #     centers[j,1] = Ycoord/denum + PERTURB*random.random()


                    for p in range(self.k):
                        curr = Pcenters[p,:]
                        pvalues = np.tile(curr, (2,1))
                        multiply = pvalues * data.T
                        centers[p,:] = np.sum(multiply, axis=1) / np.sum(curr)
                    centers = centers + PERTURB*random.random()
                    print(centers)
                    itr += 1
             self.Beta *= self.alpha
        Pcenters = np.around(Pcenters)
        return (centers, Pcenters)



#Import the data set
df1 = pd.read_csv("Datasets/ipl.csv")
df1 = df1[['one', 'two']]
data_arr = rawData.values

df2 = pd.read_table("Datasets/s1.txt", delim_whitespace=True, header=None).astype(float)
df2 = df2 * 1e-4
data_arr = df2.values

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
