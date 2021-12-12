import numpy as np
from timeit import default_timer as timer

# retorna os indexes referentes a lista de nós de parâmetro
def getDuplicatedNodes(nodeList: list):
    duplicatedNodes = np.array([], dtype=int)

    for i in range(len(nodeList)):
        nodeToBeChecked = nodeList[i]
        nodeList.count(nodeToBeChecked)

        if nodeList.count(nodeToBeChecked) > 1:
            duplicatedNodes = np.append(duplicatedNodes, i)

    return duplicatedNodes

# retorn uma lista com os nós referentes à malha de colocação
def generateAuxiliaryMesh(elementsList: list, duplicatedNodes: list, geometricNodes: list):
    auxiliaryMesh = []

    for i in range(len(elementsList)):
        elementAuxiliaryCoordinates = elementsList[i].getAuxiliaryNodesCoordinates(duplicatedNodes, geometricNodes)
        
        for j in range(len(elementAuxiliaryCoordinates)):
            if elementAuxiliaryCoordinates[j] not in auxiliaryMesh:
                auxiliaryMesh.append(elementAuxiliaryCoordinates[j])

    return auxiliaryMesh

def handleAuxiliaryMesh(elements: list, geometricNodes: list):
    start = timer()

    # elementsList = getElementsList(elements)
    duplicatedNodes = getDuplicatedNodes(geometricNodes)
    auxiliaryMesh = generateAuxiliaryMesh(elements, duplicatedNodes, geometricNodes)

    end = timer()
    print("  2 - Preparação da malha e dos elementos: ", "%.5f" % (end - start), " segundos.\n")

    return duplicatedNodes, auxiliaryMesh