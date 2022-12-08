import numpy as np

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

def generateMST(twoPinList):
    pinList = []
    for i in range(len(twoPinList)):
        pinList.append(twoPinList[i][0])
        pinList.append(twoPinList[i][1])

    pinList = list(set(pinList))
    # print(pinList)

    # Set Manhattan distance as weights of tree
    recDist = np.zeros((len(pinList),len(pinList)))
    for i in range(len(pinList)):
        for j in range(len(pinList)):
            recDist[i,j] = int(np.abs(pinList[i][0]-pinList[j][0])+ np.abs(pinList[i][1]-pinList[j][1])\
               +np.abs(pinList[i][2] - pinList[j][2]))

    X = csr_matrix(recDist)
    Tcsr = minimum_spanning_tree(X)
    Tree = Tcsr.toarray().astype(int)
    # print(Tree)

    twoPinListSorted = []
    for i in range(Tree.shape[0]):
        for j in range(Tree.shape[1]):
            if Tree[i,j] != 0:
                twoPinListSorted.append([pinList[i],pinList[j]])
    #print ('Sorted Two pin list: ',twoPinListSorted)
    return twoPinListSorted


if __name__ == '__main__':

    # X = csr_matrix([[0, 8, 0, 3],
    #                 [0, 0, 2, 5],
    #                 [0, 0, 0, 6],
    #                 [0, 0, 0, 0]])
    #
    # Tcsr = minimum_spanning_tree(X)
    # print(Tcsr.toarray().astype(int))


    # Test Generate MST
    twoPinList = [[(0,0,0),(2,0,0)],[(2,0,0),(2,0,1)],[(2,0,1),(2,2,0)],[(2,2,0),(0,0,1)]]
    MST = generateMST(twoPinList)
    print(MST)
