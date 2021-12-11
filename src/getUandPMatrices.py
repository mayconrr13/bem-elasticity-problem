from sympy.integrals.quadrature import gauss_legendre
import numpy as np
import math
from getDuplicatedNodes import getDuplicatedNodes
from getElementList import getElementsList
from generateAuxiliaryMesh import generateAuxiliaryMesh
from getIntegrationPointCoordinates import getIntegrationPointCoordinates
from getSourcePointsOutOfDomain import getSourcePoints
from pointsProperties import getPointProperties
from readInputFile import readInputFile
from shapeFunctions import getShapeFunctionValueOnNode

prescribedU, prescribedP, geometricNodes, _, elements = readInputFile("src/ex1_inputFile.txt")

elementsList = getElementsList(elements)

duplicatedNodes = getDuplicatedNodes(geometricNodes)
auxiliaryMesh = generateAuxiliaryMesh(elementsList, duplicatedNodes, geometricNodes)
sourcePoints = getSourcePoints(duplicatedNodes, geometricNodes, elementsList)

dirac = [[1,0],[0,1]]

def getRadius(sourcePoint: list, integrationPointCoordinates: list):

    xComponent = integrationPointCoordinates[0] - sourcePoint[0]
    yComponent = integrationPointCoordinates[1] - sourcePoint[1]
    radius = (xComponent ** 2 + yComponent ** 2) ** (1/2)

    return [[xComponent, yComponent], radius]

def getHandGMatrices():
    HMatrix = np.zeros((2 * len(sourcePoints), 2 * len(sourcePoints)))
    GMatrix = np.zeros((2 * len(sourcePoints), 2 * len(sourcePoints)))

    integrationPoints, weights = gauss_legendre(12, 5)
    G = 0.5 # kN/m2
    poisson = 0

    for sp in range(len(sourcePoints)):
        for el in range(len(elementsList)):
            elementNodes = elementsList[el].getElementNodesRealCoordinates(geometricNodes)
            dimensionlessPoints = elementsList[el].getDimensionlessPointsBasedOnGeometricCoordinates()
            auxiliaryDimensionlessPoints = elementsList[el].getDimensionlessPointsBasedOnElementContinuity(duplicatedNodes)

            DH = np.zeros((2, 2 * len(elementNodes)))
            DG = np.zeros((2, 2 * len(elementNodes)))
                        
            for ip in range(len(integrationPoints)):           
                integrationPointsRealCoordinates = getIntegrationPointCoordinates(integrationPoints[ip], elementNodes, dimensionlessPoints)

                integrationPointRadius = getRadius(sourcePoints[sp], integrationPointsRealCoordinates)
                radius = integrationPointRadius[1]
                radiusDiff = [integrationPointRadius[0][0] / integrationPointRadius[1], integrationPointRadius[0][1] / integrationPointRadius[1]]

                _, normalVector, jacobian = getPointProperties(ip, elementNodes, dimensionlessPoints)
                
                DRDN = radiusDiff[0] * normalVector[0] +  radiusDiff[1] * normalVector[1]
                
                U = np.zeros((2,2))
                P = np.zeros((2,2))

                for i in range(2):
                    for j in range(2):
                        U[i,j] = (1 / (8 * math.pi * G * (1 - poisson))) * ((- 3 + 4 * poisson) * math.log(radius) * dirac[i][j] + radiusDiff[i] * radiusDiff[j])
                        P[i,j] = (-1 / (4 * math.pi * (1 - poisson) * radius)) * (DRDN * ((1 - 2 * poisson) * dirac[i][j] + 2 * radiusDiff[i] * radiusDiff[j]) + (1 - 2 * poisson) * (normalVector[i] * radiusDiff[j] - normalVector[j] * radiusDiff[i]))
                                
                fi = np.zeros((2, 2 * len(elementNodes))) 

                for en in range(len(elementNodes)):
                    shapeFunctionValueOnIP = getShapeFunctionValueOnNode(integrationPoints[ip], en, auxiliaryDimensionlessPoints)
                    
                    fi[0, 2 * en] = shapeFunctionValueOnIP
                    fi[1, 2 * en + 1] = shapeFunctionValueOnIP
                
                DH = DH + np.dot(P, fi) * jacobian * weights[ip]
                DG = DG + np.dot(U, fi) * jacobian * weights[ip]

            # colocar DH e DG na matriz global
            for en in range(len(elementNodes)):
                for i in range(2):
                    HMatrix[2*sp][2 * elementsList[el].nodeList[en] + i] += DH[0][2 * en + i]
                    HMatrix[2*sp + 1][2 * elementsList[el].nodeList[en] + i] += DH[1][2 * en + i]

                    GMatrix[2*sp][2 * elementsList[el].nodeList[en] + i] += DG[0][2 * en + i]
                    GMatrix[2*sp + 1][2 * elementsList[el].nodeList[en] + i] += DG[1][2 * en + i]
                                
    return HMatrix, GMatrix

HMatrix, GMatrix = getHandGMatrices()

def applyBoundaryConditions(HMatrix, GMatrix, prescribedU, prescribedP, sourcePoints, auxiliaryMesh):
    FHMatrix = HMatrix
    FGMatrix = GMatrix
    FVector = np.zeros(2 * len(sourcePoints), dtype=float)

    for j in range(len(prescribedU)):
        index = prescribedU[j][0]

        FHMatrix[:, 2 * index] = - GMatrix[:, 2 * index]
        FHMatrix[:, 2 * index + 1] = - GMatrix[:, 2 * index + 1]

        FGMatrix[:, 2 * index] = - HMatrix[:, 2 * index]
        FGMatrix[:, 2 * index + 1] = - HMatrix[:, 2 * index + 1]

        FVector[2 * index] += prescribedU[j][1][0]
        FVector[2 * index + 1] += prescribedU[j][1][1]

    for k in range(len(prescribedP)):
        index = prescribedP[k][0]

        FVector[2 * index] += prescribedP[k][1][0]
        FVector[2 * index + 1] += prescribedP[k][1][1]
    
    return FHMatrix, FGMatrix, FVector

FHMatrix, FGMatrix, FVector = applyBoundaryConditions(HMatrix, GMatrix, prescribedU, prescribedP, sourcePoints, auxiliaryMesh)

resultsVector = np.linalg.solve(FHMatrix, np.dot(FGMatrix, FVector))

print(resultsVector)