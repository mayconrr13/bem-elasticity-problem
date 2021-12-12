from timeit import default_timer as timer
from src.readInputFile import readInputFile
from src.handleAuxiliaryMesh import handleAuxiliaryMesh
from src.solveElasticityBoundaryProblem import solveElasticityBoundaryProblem
from src.createParaviewFile import createParaviewFile
import matplotlib.pyplot as plt

def elasticityProblemBEM(file: str):
    print("Início do processo")
    start = timer()
    # Leitura do arquivo de entrada
    (prescribedDisplacements, 
    prescribedForces, 
    material, 
    geometricNodes, 
    internalPoints, 
    elements) = readInputFile(file)  

    # Cria malha de colocação
    duplicatedNodes, auxiliaryMesh = handleAuxiliaryMesh(elements, geometricNodes)

    # Resolução do problema de elasticidade por MEC   
    (
        boundaryDisplacements, 
        _, 
        internalDisplacements, 
        internalStress
    ) = solveElasticityBoundaryProblem(
        prescribedDisplacements, 
        prescribedForces, 
        material, 
        geometricNodes, 
        internalPoints, 
        elements, 
        duplicatedNodes, 
        auxiliaryMesh,
        12
    )
    print(boundaryDisplacements)

    # Criação do arquivo de saíde em Paraview
    createParaviewFile(boundaryDisplacements, internalDisplacements, internalStress,auxiliaryMesh, internalPoints, elements)

    end = timer()
    print("Fim do processo. Tempo total: ", "%.5f" % (end - start), " segundos.")

    return

elasticityProblemBEM("src/generateInputFile/squarePlate1EO2.txt")