import numpy as np

# retorna os indexes referentes a lista de nós de parâmetro
def getDuplicatedNodes(nodeList: list):
    duplicatedNodes = np.array([], dtype=int)

    for i in range(len(nodeList)):
        nodeToBeChecked = nodeList[i]
        nodeList.count(nodeToBeChecked)

        if nodeList.count(nodeToBeChecked) > 1:
            duplicatedNodes = np.append(duplicatedNodes, i)

    return duplicatedNodes

# nodeList = [[0,0],[1,0],[2,0],[2,0],[2,1],[2,2]]
# print(len(getDuplicatedNodes(nodeList)))