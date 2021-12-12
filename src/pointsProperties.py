from src.diffShapeFunction import *

def getPointProperties(point: float, elementNodes: list, adimentionalPoints: list):
    tangentVector = [0, 0]
    normalVector = [0, 0]    

    for j in range(len(elementNodes)):
        tangentVector[0] += getDiffShapeFunction(point, j, adimentionalPoints) * elementNodes[j][0]
        tangentVector[1] += getDiffShapeFunction(point, j, adimentionalPoints) * elementNodes[j][1]
    
    jacobian = (tangentVector[0] * tangentVector[0] + tangentVector[1] * tangentVector[1]) ** (1 / 2)

    normalVector[0] = tangentVector[1] / jacobian
    normalVector[1] = - tangentVector[0] / jacobian

    return tangentVector, normalVector, jacobian

# teste = getPointProperties(1, [[0,0], [1,1], [2,2]], [-1,0,1])
# print(teste)