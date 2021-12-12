from timeit import default_timer as timer
from src.readInputFile import readInputFile
from src.solveElasticityBoundaryProblem import solveElasticityBoundaryProblem
from src.createParaviewFile import createParaviewFile

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

    # Resolução do problema de elasticidade por MEC   
    (
        boundaryDisplacements, 
        boundaryForces, 
        internalDisplacements, 
        internalStress
    ) = solveElasticityBoundaryProblem(
        prescribedDisplacements, 
        prescribedForces, 
        material, 
        geometricNodes, 
        internalPoints, 
        elements,
        12
    )

    # Criação do arquivo de saíde em Paraview
    createParaviewFile(boundaryDisplacements, boundaryForces, internalDisplacements, internalStress)

    end = timer()
    print("Fim do processo. Tempo total: ", "%.5f" % (end - start), " segundos.")

    return

elasticityProblemBEM("src/ex1_inputFileEP.txt")