# Retorna o rário e seus componentes tomando dois pontos como referência (final e inicial)
def getRadius(sourcePoint: list, integrationPointCoordinates: list):

    xComponent = integrationPointCoordinates[0] - sourcePoint[0]
    yComponent = integrationPointCoordinates[1] - sourcePoint[1]
    radius = (xComponent ** 2 + yComponent ** 2) ** (1/2)

    return [[xComponent, yComponent], radius]

# verifica se o ponto fonte está sobre o elemento
def handleSourcePointOnElement(element, sourcePointCoordinate, duplicatedNodes, geometricNodes):
    auxiliaryNodesCoordinates = element.getAuxiliaryNodesCoordinates(duplicatedNodes, geometricNodes)

    if sourcePointCoordinate in auxiliaryNodesCoordinates:
        return True
    else:
        return False

# retorna o índice de um array
def getIndex(newList, parameter):
    for i in range(len(newList)):
        if parameter == newList[i]:
            return i
