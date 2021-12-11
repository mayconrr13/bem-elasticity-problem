# retorn uma lista com os nós referentes à malha de colocação
from Element import Element

def generateAuxiliaryMesh(elementsList: list, duplicatedNodes: list, geometricNodes: list):
    auxiliaryMesh = []

    for i in range(len(elementsList)):
        elementAuxiliaryCoordinates = elementsList[i].getAuxiliaryNodesCoordinates(duplicatedNodes, geometricNodes)
        
        for j in range(len(elementAuxiliaryCoordinates)):
            if elementAuxiliaryCoordinates[j] not in auxiliaryMesh:
                auxiliaryMesh.append(elementAuxiliaryCoordinates[j])

    return auxiliaryMesh
    
# elementList = [Element([0,1,2]), Element([3,4,5])]
# geometricNodes = [[0,0],[1,0],[2,0],[2,0],[2,1],[2,2]]
# print(generateAuxiliaryMesh(elementList, [2,3], geometricNodes))