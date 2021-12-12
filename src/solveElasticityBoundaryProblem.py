import numpy as np
from src.handleAuxiliaryMesh import handleAuxiliaryMesh
from src.getHandGMatrices import solveBoundaryProblem, solveInternalPoints

def solveElasticityBoundaryProblem(
    prescribedDisplacements: list, 
    prescribedForces: list, 
    material: list, 
    geometricNodes: list, 
    internalPoints: list, 
    elements: list,
    numberOfIntegrationPoints: int
):
    poisson = material[1]
    G = material[0] / (2 * (1 + poisson))
    integrationPoints, weights = np.polynomial.legendre.leggauss(numberOfIntegrationPoints)

    # Cria malha de colocação
    duplicatedNodes, auxiliaryMesh = handleAuxiliaryMesh(elements, geometricNodes)
    sourcePoints = auxiliaryMesh

    # Resolução do problema    
    boundaryDisplacements, boundaryForces = solveBoundaryProblem(sourcePoints, prescribedDisplacements, prescribedForces, auxiliaryMesh, duplicatedNodes, elements, geometricNodes, integrationPoints, weights, poisson, G)
    internalDisplacements, internalStress = solveInternalPoints(internalPoints, boundaryForces, boundaryDisplacements, auxiliaryMesh, duplicatedNodes, elements, geometricNodes, integrationPoints, weights, poisson, G)

    return boundaryDisplacements, boundaryForces, internalDisplacements, internalStress