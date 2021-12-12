import numpy as np
import math
from timeit import default_timer as timer

from src.getIntegrationPointCoordinates import getIntegrationPointCoordinates
from src.pointsProperties import getPointProperties
from src.shapeFunctions import getShapeFunctionValueOnNode
from src.auxiliaryFunctions import * 

dirac = [[1,0],[0,1]]

#### Resolução do problema no contorno
# H[] e G[]
def getHandGMatrices(sourcePoints: list, auxiliaryMesh: list, duplicatedNodes: list, elementsList: list, geometricNodes: list, integrationPoints: list, weights: list, poisson: float, G:float):
    start = timer()

    HMatrix = np.zeros((2 * len(sourcePoints), 2 * len(auxiliaryMesh)))
    GMatrix = np.zeros((2 * len(sourcePoints), 2 * len(auxiliaryMesh)))

    # verificação de pontos fontes no contorno
    if sourcePoints == auxiliaryMesh:
        HMatrix += np.identity(2 * len(sourcePoints)) * 1/2

    for sp in range(len(sourcePoints)):
        for el in range(len(elementsList)):
            elementNodes = elementsList[el].getElementNodesRealCoordinates(geometricNodes)
            dimensionlessPoints = elementsList[el].getDimensionlessPointsBasedOnGeometricCoordinates()
            auxiliaryDimensionlessPoints = elementsList[el].getDimensionlessPointsBasedOnElementContinuity(duplicatedNodes)
            auxiliaryCoordinates = elementsList[el].getAuxiliaryNodesCoordinates(duplicatedNodes, geometricNodes)
            nodeIndex = getIndex(list(auxiliaryCoordinates), sourcePoints[sp])

            sourcePointIsOnElement = handleSourcePointOnElement(elementsList[el], sourcePoints[sp], duplicatedNodes, geometricNodes)

            DH = np.zeros((2, 2 * len(elementNodes)))
            DG = np.zeros((2, 2 * len(elementNodes)))
            displacementCPVContribution = np.zeros((2,2 * len(elementNodes)))
            forceCPVContribution = np.zeros((2,2 * len(elementNodes)))

            if sourcePointIsOnElement:
                displacementCPV = np.zeros((2,2))
                forceCPV = np.zeros((2,2))

                sourcePointDimensionlessCoordinate = auxiliaryDimensionlessPoints[nodeIndex]
                fi_ = np.zeros((2, 2 * len(elementNodes)))

                for en in range(len(elementNodes)):
                    shapeFunctionValueOnSource = getShapeFunctionValueOnNode(sourcePointDimensionlessCoordinate, en, auxiliaryDimensionlessPoints)
                    
                    fi_[0, 2 * en] = shapeFunctionValueOnSource
                    fi_[1, 2 * en + 1] = shapeFunctionValueOnSource

                # CPV
                # pontos fonte nas extremidades -1 ou 1
                sourceTangentVector, sourceNormalVector, sourceJacobain = getPointProperties(sourcePointDimensionlessCoordinate, elementNodes, dimensionlessPoints)
                radiusDotDiff = [sourceTangentVector[0] / sourceJacobain, sourceTangentVector[1] / sourceJacobain]

                if sourcePointDimensionlessCoordinate == -1 or sourcePointDimensionlessCoordinate == 1:
                    displacementCPV[0][0] = (((-3 + 4 * poisson) * sourceJacobain) / (8 * math.pi * G * (1 - poisson))) * (2 * math.log(abs(2 * sourceJacobain)) - 2)
                    displacementCPV[1][1] = displacementCPV[0][0]


                    if sourcePointDimensionlessCoordinate == -1:
                        forceCPV[0,0] = (((-1 + 2 * poisson) * (sourceNormalVector[0] * radiusDotDiff[0] - sourceNormalVector[0] * radiusDotDiff[0])) / (4 * math.pi * (1 - poisson))) * math.log(2)
                        forceCPV[0,1] = (((-1 + 2 * poisson) * (sourceNormalVector[0] * radiusDotDiff[1] - sourceNormalVector[1] * radiusDotDiff[0])) / (4 * math.pi * (1 - poisson))) * math.log(2)
                        forceCPV[1,0] = (((-1 + 2 * poisson) * (sourceNormalVector[1] * radiusDotDiff[0] - sourceNormalVector[0] * radiusDotDiff[1])) / (4 * math.pi * (1 - poisson))) * math.log(2)
                        forceCPV[1,1] = (((-1 + 2 * poisson) * (sourceNormalVector[1] * radiusDotDiff[1] - sourceNormalVector[1] * radiusDotDiff[1])) / (4 * math.pi * (1 - poisson))) * math.log(2)
                        
                    else:
                        forceCPV[0,0] = (((-1 + 2 * poisson) * (sourceNormalVector[0] * radiusDotDiff[0] - sourceNormalVector[0] * radiusDotDiff[0])) / (4 * math.pi * (1 - poisson))) * math.log(2) * -1
                        forceCPV[0,1] = (((-1 + 2 * poisson) * (sourceNormalVector[0] * radiusDotDiff[1] - sourceNormalVector[1] * radiusDotDiff[0])) / (4 * math.pi * (1 - poisson))) * math.log(2) * -1
                        forceCPV[1,0] = (((-1 + 2 * poisson) * (sourceNormalVector[1] * radiusDotDiff[0] - sourceNormalVector[0] * radiusDotDiff[1])) / (4 * math.pi * (1 - poisson))) * math.log(2) * -1
                        forceCPV[1,1] = (((-1 + 2 * poisson) * (sourceNormalVector[1] * radiusDotDiff[1] - sourceNormalVector[1] * radiusDotDiff[1])) / (4 * math.pi * (1 - poisson))) * math.log(2) * -1
                                    
                # ponto fonte na parte interna diferente de -1 e 1
                else:
                    partial1 = (1 + sourcePointDimensionlessCoordinate) * math.log(abs(sourceJacobain *  (1 + sourcePointDimensionlessCoordinate)))
                    partial2 = (1 - sourcePointDimensionlessCoordinate) * math.log(abs(sourceJacobain *  (1 - sourcePointDimensionlessCoordinate)))
                    displacementCPV[0][0] = (((-3 + 4 * poisson) * sourceJacobain) / (8 * math.pi * G * (1 - poisson))) * (partial1 + partial2 - 2)
                    displacementCPV[1][1] = displacementCPV[0][0]

                    forceCPV[0,0] = (((-1 + 2 * poisson) * (sourceNormalVector[0] * radiusDotDiff[0] - sourceNormalVector[0] * radiusDotDiff[0])) / (4 * math.pi * (1 - poisson))) * (math.log(1 - sourcePointDimensionlessCoordinate) - math.log(1 + sourcePointDimensionlessCoordinate))
                    forceCPV[0,1] = (((-1 + 2 * poisson) * (sourceNormalVector[0] * radiusDotDiff[1] - sourceNormalVector[1] * radiusDotDiff[0])) / (4 * math.pi * (1 - poisson))) * (math.log(1 - sourcePointDimensionlessCoordinate) - math.log(1 + sourcePointDimensionlessCoordinate))
                    forceCPV[1,0] = (((-1 + 2 * poisson) * (sourceNormalVector[1] * radiusDotDiff[0] - sourceNormalVector[0] * radiusDotDiff[1])) / (4 * math.pi * (1 - poisson))) * (math.log(1 - sourcePointDimensionlessCoordinate) - math.log(1 + sourcePointDimensionlessCoordinate))
                    forceCPV[1,1] = (((-1 + 2 * poisson) * (sourceNormalVector[1] * radiusDotDiff[1] - sourceNormalVector[1] * radiusDotDiff[1])) / (4 * math.pi * (1 - poisson))) * (math.log(1 - sourcePointDimensionlessCoordinate) - math.log(1 + sourcePointDimensionlessCoordinate))
                    
                # print("Pvpc: ", forceCPV)
                displacementCPVContribution = np.dot(displacementCPV, fi_)  
                forceCPVContribution = np.dot(forceCPV, fi_)
                                    
                        
            for ip in range(len(integrationPoints)):           
                integrationPointsRealCoordinates = getIntegrationPointCoordinates(integrationPoints[ip], elementNodes, dimensionlessPoints)

                integrationPointRadius = getRadius(sourcePoints[sp], integrationPointsRealCoordinates)
                radius = integrationPointRadius[1]
                radiusDiff = [integrationPointRadius[0][0] / integrationPointRadius[1], integrationPointRadius[0][1] / integrationPointRadius[1]]

                _, normalVector, jacobian = getPointProperties(ip, elementNodes, dimensionlessPoints)
                
                DRDN = radiusDiff[0] * normalVector[0] +  radiusDiff[1] * normalVector[1]
                
                U = np.zeros((2,2))
                P = np.zeros((2,2))
                _U = np.zeros((2,2))
                _P = np.zeros((2,2))

                for i in range(2):
                    for j in range(2):
                        # remover o termo final de U
                        U[i,j] += (1 / (8 * math.pi * G * (1 - poisson))) * ((- 3 + 4 * poisson) * math.log(radius) * dirac[i][j] + radiusDiff[i] * radiusDiff[j] - (0.5 * (7 - 8 * poisson)) * dirac[i][j])
                        P[i,j] += (-1 / (4 * math.pi * (1 - poisson) * radius)) * (DRDN * ((1 - 2 * poisson) * dirac[i][j] + 2 * radiusDiff[i] * radiusDiff[j]) + (1 - 2 * poisson) * (normalVector[i] * radiusDiff[j] - normalVector[j] * radiusDiff[i]))
                        
                        
                        # contribuir com parcela de regularização
                        if sourcePointIsOnElement:                            
                            radiusDot = jacobian * (integrationPoints[ip] - auxiliaryDimensionlessPoints[nodeIndex])

                            sourceTangentVector, sourceNormalVector, sourceJacobain = getPointProperties(auxiliaryDimensionlessPoints[nodeIndex], elementNodes, dimensionlessPoints)
                            radiusDotDiff = [sourceTangentVector[0] / sourceJacobain, sourceTangentVector[1] / sourceJacobain]
                            
                            _U[i,j] += (-1 / (8 * math.pi * G * (1 - poisson))) * ((- 3 + 4 * poisson) * math.log(abs(radiusDot)) * dirac[i][j])
                            _P[i,j] += (1 / (4 * math.pi * (1 - poisson) * radiusDot)) * (1 - 2 * poisson) * (sourceNormalVector[i] * radiusDotDiff[j] - sourceNormalVector[j] * radiusDotDiff[i])

                fi = np.zeros((2, 2 * len(elementNodes))) 
                _fi = np.zeros((2, 2 * len(elementNodes))) 

                for en in range(len(elementNodes)):
                    shapeFunctionValueOnIP = getShapeFunctionValueOnNode(integrationPoints[ip], en, auxiliaryDimensionlessPoints)
                    
                    fi[0, 2 * en] = shapeFunctionValueOnIP
                    fi[1, 2 * en + 1] = shapeFunctionValueOnIP

                    if sourcePointIsOnElement:
                        shapeFunctionValueOnSource = getShapeFunctionValueOnNode(sourcePointDimensionlessCoordinate, en, auxiliaryDimensionlessPoints)
                        
                        _fi[0, 2 * en] = shapeFunctionValueOnSource
                        _fi[1, 2 * en + 1] = shapeFunctionValueOnSource
                
                DH = DH + np.dot(P, fi) * jacobian * weights[ip] + np.dot(_P, _fi) * jacobian * weights[ip]
                DG = DG + np.dot(U, fi) * jacobian * weights[ip] + np.dot(_U, _fi) * jacobian * weights[ip]
            
            # colocar DH e DG na matriz global
            for en in range(len(elementNodes)):
                for i in range(2):
                    HMatrix[2*sp][2 * elementsList[el].nodeList[en] + i] += DH[0][2 * en + i] + forceCPVContribution[0][2 * en + i]
                    HMatrix[2*sp + 1][2 * elementsList[el].nodeList[en] + i] += DH[1][2 * en + i] + forceCPVContribution[1][2 * en + i]

                    GMatrix[2*sp][2 * elementsList[el].nodeList[en] + i] += DG[0][2 * en + i] + displacementCPVContribution[0][2 * en + i]
                    GMatrix[2*sp + 1][2 * elementsList[el].nodeList[en] + i] += DG[1][2 * en + i] + displacementCPVContribution[1][2 * en + i]

    end = timer()
    print("    -- Geração das matrizes H e G: ", "%.5f" % (end - start), " segundos.")     

    return HMatrix, GMatrix

# Aplicação da condições de contorno
def applyBoundaryConditions(HMatrix, GMatrix, prescribedU, prescribedP, sourcePoints):
    start = timer()

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

    end = timer()
    print("    -- Aplicação da condições de contorno: ", "%.5f" % (end - start), " segundos.")  
    
    return FHMatrix, FGMatrix, FVector

# criação do vertor de forças de superfície e de deslocamento no contorno
def handleResultantVectors(resultsVector, FVector, prescribedU):
    start = timer()

    boundaryDisplacements = np.array(resultsVector)
    boundaryForces = np.array(FVector)

    for i in range(len(prescribedU)):
        index = prescribedU[i][0] 

        boundaryForces[2 * index] = resultsVector[2 * index]
        boundaryForces[2 * index + 1] = resultsVector[2 * index + 1] 

        boundaryDisplacements[2 * index] = prescribedU[i][1][0]
        boundaryDisplacements[2 * index + 1] = prescribedU[i][1][1]

    end = timer()
    print("    -- Obtenção dos vetores resultantes: ", "%.5f" % (end - start), " segundos.")  
    
    return boundaryDisplacements, boundaryForces

# Resolução do problema no contorno
def solveBoundaryProblem(sourcePoints: list, prescribedU: list, prescribedP: list, auxiliaryMesh: list, duplicatedNodes: list, elementsList: list, geometricNodes: list, integrationPoints: list, weights: list, poisson: float, G: float):
    print("  3 - Resolução do problema no contorno")
    start = timer()

    HMatrix, GMatrix = getHandGMatrices(sourcePoints, auxiliaryMesh, duplicatedNodes, elementsList, geometricNodes, integrationPoints, weights, poisson, G)
    FHMatrix, FGMatrix, FVector = applyBoundaryConditions(HMatrix, GMatrix, prescribedU, prescribedP, sourcePoints)
    resultsVector = np.linalg.solve(FHMatrix, np.dot(FGMatrix, FVector))
    boundaryDisplacements, boundaryForces = handleResultantVectors(resultsVector, FVector, prescribedU)

    end = timer()
    print("    -- Tempo total: ", "%.5f" % (end - start), " segundos.\n")

    return boundaryDisplacements, boundaryForces

#### Resolução do problema nos pontos internos
# D[] e S[]
def getDandSMatrices(sourcePoints, elementsList, auxiliaryMesh, duplicatedNodes, geometricNodes, integrationPoints, weights, poisson, G):
    start = timer()
    
    DMatrix = np.zeros((4 * len(sourcePoints), 2 * len(auxiliaryMesh)))
    SMatrix = np.zeros((4 * len(sourcePoints), 2 * len(auxiliaryMesh)))

    for sp in range(len(sourcePoints)):
        for el in range(len(elementsList)):
            elementNodes = elementsList[el].getElementNodesRealCoordinates(geometricNodes)
            dimensionlessPoints = elementsList[el].getDimensionlessPointsBasedOnGeometricCoordinates()
            auxiliaryDimensionlessPoints = elementsList[el].getDimensionlessPointsBasedOnElementContinuity(duplicatedNodes)

            Dcontribution = np.zeros((4,2 * len(elementNodes)))
            Scontribution = np.zeros((4,2 * len(elementNodes)))   

            for ip in range(len(integrationPoints)):           
                integrationPointsRealCoordinates = getIntegrationPointCoordinates(integrationPoints[ip], elementNodes, dimensionlessPoints)

                integrationPointRadius = getRadius(sourcePoints[sp], integrationPointsRealCoordinates)
                radius = integrationPointRadius[1]
                radiusDiff = [integrationPointRadius[0][0] / integrationPointRadius[1], integrationPointRadius[0][1] / integrationPointRadius[1]]

                _, normalVector, jacobian = getPointProperties(ip, elementNodes, dimensionlessPoints)
                
                DRDN = radiusDiff[0] * normalVector[0] +  radiusDiff[1] * normalVector[1]
                
                D = np.zeros((4,2))
                S = np.zeros((4,2))

                # preenchimento das matrizes D e S
                for k in range(2):
                    D[0, k] += (1 / (4 * math.pi * (1 - poisson) * radius)) * ((1 - 2 * poisson) * (dirac[k][0] * radiusDiff[0] + dirac[k][0] * radiusDiff[0] - dirac[0][0] * radiusDiff[k]) + 2 * radiusDiff[0] * radiusDiff[0] * radiusDiff[k]) #i=0 e j=0
                    D[1, k] += (1 / (4 * math.pi * (1 - poisson) * radius)) * ((1 - 2 * poisson) * (dirac[k][0] * radiusDiff[1] + dirac[k][1] * radiusDiff[0] - dirac[0][1] * radiusDiff[k]) + 2 * radiusDiff[0] * radiusDiff[1] * radiusDiff[k]) #i=0 e j=1
                    D[2, k] += (1 / (4 * math.pi * (1 - poisson) * radius)) * ((1 - 2 * poisson) * (dirac[k][1] * radiusDiff[0] + dirac[k][0] * radiusDiff[1] - dirac[1][0] * radiusDiff[k]) + 2 * radiusDiff[1] * radiusDiff[0] * radiusDiff[k]) #i=1 e j=0
                    D[3, k] += (1 / (4 * math.pi * (1 - poisson) * radius)) * ((1 - 2 * poisson) * (dirac[k][1] * radiusDiff[1] + dirac[k][1] * radiusDiff[1] - dirac[1][1] * radiusDiff[k]) + 2 * radiusDiff[1] * radiusDiff[1] * radiusDiff[k]) #i=1 e j=1

                    # k0, i0, j0
                    partial1S0 = 2 * DRDN * ((1 - 2 * poisson) * dirac[0][0] * radiusDiff[k] + poisson * (dirac[0][k] * radiusDiff[0] + dirac[0][k] * radiusDiff[0]) - 4 * radiusDiff[0] * radiusDiff[0] * radiusDiff[k])
                    partial2S0 = 2 * poisson * (normalVector[0] * radiusDiff[0] * radiusDiff[k] + normalVector[0] * radiusDiff[0] * radiusDiff[k])
                    partial3S0 = (1 - 2 * poisson) * (2 * normalVector[k] * radiusDiff[0] * radiusDiff[0] + normalVector[0] * dirac[0][k] + normalVector[0] * dirac[0][k])
                    partial4S0 = (- 1 + 4 * poisson) * normalVector[k] * dirac[0][0] 
                    S[0, k] += (G / (2 * math.pi * (1 - poisson) * radius * radius)) * (partial1S0 + partial2S0 + partial3S0 + partial4S0)

                    # k0, i0, j1
                    partial1S1 = 2 * DRDN * ((1 - 2 * poisson) * dirac[0][1] * radiusDiff[k] + poisson * (dirac[0][k] * radiusDiff[1] + dirac[1][k] * radiusDiff[0]) - 4 * radiusDiff[0] * radiusDiff[1] * radiusDiff[k])
                    partial2S1 = 2 * poisson * (normalVector[0] * radiusDiff[1] * radiusDiff[k] + normalVector[1] * radiusDiff[0] * radiusDiff[k])
                    partial3S1 = (1 - 2 * poisson) * (2 * normalVector[k] * radiusDiff[0] * radiusDiff[1] + normalVector[1] * dirac[0][k] + normalVector[0] * dirac[1][k])
                    partial4S1 = (- 1 + 4 * poisson) * normalVector[k] * dirac[0][1] 
                    S[1, k] += (G / (2 * math.pi * (1 - poisson) * radius * radius)) * (partial1S1 + partial2S1 + partial3S1 + partial4S1)

                    # k0, i1, j0
                    partial1S2 = 2 * DRDN * ((1 - 2 * poisson) * dirac[1][0] * radiusDiff[k] + poisson * (dirac[1][k] * radiusDiff[0] + dirac[0][k] * radiusDiff[1]) - 4 * radiusDiff[1] * radiusDiff[0] * radiusDiff[k])
                    partial2S2 = 2 * poisson * (normalVector[1] * radiusDiff[0] * radiusDiff[k] + normalVector[0] * radiusDiff[1] * radiusDiff[k])
                    partial3S2 = (1 - 2 * poisson) * (2 * normalVector[k] * radiusDiff[1] * radiusDiff[0] + normalVector[0] * dirac[1][k] + normalVector[1] * dirac[0][k])
                    partial4S2 = (- 1 + 4 * poisson) * normalVector[k] * dirac[1][0] 
                    S[2, k] += (G / (2 * math.pi * (1 - poisson) * radius * radius)) * (partial1S2 + partial2S2 + partial3S2 + partial4S2)

                    # k0, i1, j1
                    partial1S3 = 2 * DRDN * ((1 - 2 * poisson) * dirac[1][0] * radiusDiff[k] + poisson * (dirac[1][k] * radiusDiff[1] + dirac[1][k] * radiusDiff[1]) - 4 * radiusDiff[1] * radiusDiff[1] * radiusDiff[k])
                    partial2S3 = 2 * poisson * (normalVector[1] * radiusDiff[1] * radiusDiff[k] + normalVector[1] * radiusDiff[1] * radiusDiff[k])
                    partial3S3 = (1 - 2 * poisson) * (2 * normalVector[k] * radiusDiff[1] * radiusDiff[1] + normalVector[1] * dirac[1][k] + normalVector[1] * dirac[1][k])
                    partial4S3 = (- 1 + 4 * poisson) * normalVector[k] * dirac[1][1] 
                    S[3, k] += (G / (2 * math.pi * (1 - poisson) * radius * radius)) * (partial1S3 + partial2S3 + partial3S3 + partial4S3)
                
                fi = np.zeros((2, 2 * len(elementNodes))) 

                for en in range(len(elementNodes)):
                    shapeFunctionValueOnIP = getShapeFunctionValueOnNode(integrationPoints[ip], en, auxiliaryDimensionlessPoints)
                    
                    fi[0, 2 * en] = shapeFunctionValueOnIP
                    fi[1, 2 * en + 1] = shapeFunctionValueOnIP
                                
                Dcontribution = Dcontribution + np.dot(D, fi) * jacobian * weights[ip]
                Scontribution = Scontribution + np.dot(S, fi) * jacobian * weights[ip]
            
            # colocar DH e DG na matriz global
            for en in range(len(elementNodes)):
                for i in range(2):
                    DMatrix[4*sp][2 * elementsList[el].nodeList[en] + i] += Dcontribution[0][2 * en + i]
                    DMatrix[4*sp + 1][2 * elementsList[el].nodeList[en] + i] += Dcontribution[1][2 * en + i]
                    DMatrix[4*sp + 2][2 * elementsList[el].nodeList[en] + i] += Dcontribution[2][2 * en + i]
                    DMatrix[4*sp + 3][2 * elementsList[el].nodeList[en] + i] += Dcontribution[3][2 * en + i]

                    SMatrix[4*sp][2 * elementsList[el].nodeList[en] + i] += Scontribution[0][2 * en + i]
                    SMatrix[4*sp + 1][2 * elementsList[el].nodeList[en] + i] += Scontribution[1][2 * en + i]
                    SMatrix[4*sp + 2][2 * elementsList[el].nodeList[en] + i] += Scontribution[2][2 * en + i]
                    SMatrix[4*sp + 3][2 * elementsList[el].nodeList[en] + i] += Scontribution[3][2 * en + i]

    end = timer()
    print("    -- Geração das matrizes D e S: ", "%.5f" % (end - start), " segundos.")

    return DMatrix, SMatrix

# Deslocamentos dos pontos internos
def handleInternalDisplacements(internalPoints, boundaryForces, boundaryDisplacements, auxiliaryMesh, duplicatedNodes, elementsList, geometricNodes, integrationPoints, weights, poisson, G):
    start = timer()

    IntHMatrix, IntGMatrix = getHandGMatrices(internalPoints, auxiliaryMesh, duplicatedNodes, elementsList, geometricNodes, integrationPoints, weights, poisson, G)
    internalDisplacements = np.dot(IntGMatrix, boundaryForces) - np.dot(IntHMatrix, boundaryDisplacements)
    
    end = timer()
    print("    -- Deslocamentos dos pontos internos: ", "%.5f" % (end - start), " segundos.")

    return internalDisplacements

#Tensões nos pontos internos
def handleInternalStress(internalPoints, boundaryForces, boundaryDisplacements, elementsList, auxiliaryMesh, duplicatedNodes, geometricNodes, integrationPoints, weights, poisson, G):
    start = timer()

    DMatrix, SMatrix = getDandSMatrices(internalPoints, elementsList, auxiliaryMesh, duplicatedNodes, geometricNodes, integrationPoints, weights, poisson, G)
    internalStress = np.dot(DMatrix, boundaryForces) - np.dot(SMatrix, boundaryDisplacements)
    
    end = timer()
    print("    -- Tensões dos pontos internos: ", "%.5f" % (end - start), " segundos.")

    return internalStress

# Resolução do problema para os pontos internos
def solveInternalPoints(internalPoints, boundaryForces, boundaryDisplacements, auxiliaryMesh, duplicatedNodes, elementsList, geometricNodes, integrationPoints, weights, poisson, G):
    print("  4 - Resolução do problema no contorno")
    start = timer()

    internalDisplacements = handleInternalDisplacements(internalPoints, boundaryForces, boundaryDisplacements, auxiliaryMesh, duplicatedNodes, elementsList, geometricNodes, integrationPoints, weights, poisson, G)
    internalStress = handleInternalStress(internalPoints, boundaryForces, boundaryDisplacements, elementsList, auxiliaryMesh, duplicatedNodes, geometricNodes, integrationPoints, weights, poisson, G)

    end = timer()
    print("    -- Tempo total: ", "%.5f" % (end - start), " segundos.\n")

    return internalDisplacements, internalStress
