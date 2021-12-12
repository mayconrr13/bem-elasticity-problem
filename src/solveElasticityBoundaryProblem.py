import numpy as np
from src.getSourcePointsOutOfDomain import getSourcePoints
from src.handleAuxiliaryMesh import handleAuxiliaryMesh
from src.processSteps import solveBoundaryProblem, solveInternalPoints

def solveElasticityBoundaryProblem(
    prescribedDisplacements: list, 
    prescribedForces: list, 
    material: list, 
    geometricNodes: list, 
    internalPoints: list, 
    elements: list,
     duplicatedNodes: list, 
     auxiliaryMesh: list,
    numberOfIntegrationPoints: int
):
    poisson = material[1]
    G = material[0] / (2 * (1 + poisson))
    integrationPoints, weights = np.polynomial.legendre.leggauss(numberOfIntegrationPoints)

    sourcePoints = getSourcePoints(duplicatedNodes, geometricNodes, elements)

    # Resolução do problema    
    boundaryDisplacements, boundaryForces = solveBoundaryProblem(sourcePoints, prescribedDisplacements, prescribedForces, auxiliaryMesh, duplicatedNodes, elements, geometricNodes, integrationPoints, weights, poisson, G)
    internalDisplacements, internalStress = solveInternalPoints(internalPoints, boundaryForces, boundaryDisplacements, auxiliaryMesh, duplicatedNodes, elements, geometricNodes, integrationPoints, weights, poisson, G)
    
    return boundaryDisplacements, boundaryForces, internalDisplacements, internalStress
    