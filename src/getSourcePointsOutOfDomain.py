from Element import Element
from pointsProperties import getPointProperties

# retorna os pontos fontes localizados fora do domínio
def getSourcePoints(duplicatedNodes, geometricNodes, elementsList):
    sourcePointsList = []

    for i in range(len(elementsList)):
        elementAuxiliaryNodes = elementsList[i].getAuxiliaryNodesCoordinates(duplicatedNodes, geometricNodes)
        dimensionlessPoints = elementsList[i].getDimensionlessPointsBasedOnGeometricCoordinates()
        
        for k in range(len(elementsList[i].nodeList)):
            _, normalVector, _ = getPointProperties(k, elementAuxiliaryNodes, dimensionlessPoints)
                        
            xCoordinate = normalVector[0] * 0.25 + elementAuxiliaryNodes[k][0]
            yCoordinate = normalVector[1] * 0.25 + elementAuxiliaryNodes[k][1]

            sourcePointCoordinates = [xCoordinate, yCoordinate]

            if sourcePointsList.count(sourcePointCoordinates) == 0:
                sourcePointsList.append(sourcePointCoordinates)

    return sourcePointsList

# elementList = [Element([0,1,2]), Element([3,4,5])]
# geometricNodes = [[0,0],[1,0],[2,0],[2,0],[2,1],[2,2]]
# print(getSourcePoints([2,3], geometricNodes,elementList))