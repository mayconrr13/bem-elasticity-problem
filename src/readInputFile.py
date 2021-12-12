from timeit import default_timer as timer

from src.Element import Element

def lineToIgnore(inputFile: str):
    inputFile.readline()
    return

def readInputFile(file: str):
    start = timer()

    inputFile = open(file, "r")

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    numberOfExternalNodes = int(inputFile.readline().split(":")[1])
    numberOfInternalNodes = int(inputFile.readline().split(":")[1])   
    nodesPerElement = int(inputFile.readline().split(":")[1]) 
    numberOfElements = int(inputFile.readline().split(":")[1])
    displacements = int(inputFile.readline().split(":")[1])
    forces = int(inputFile.readline().split(":")[1])

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    geometricNodes = []    
    for i in range(numberOfExternalNodes):
        nodesCoordinates = inputFile.readline().split(",")
        geometricNodes.append([float(nodesCoordinates[1]), float(nodesCoordinates[2])])

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    internalPoints = []
    for j in range(numberOfInternalNodes):
        points = inputFile.readline().split(",")
        internalPoints.append([float(points[1]), float(points[2])])

    
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    elements = []
    for k in range(numberOfElements):
        item = inputFile.readline().split(",")
        
        element = []
        for kk in range(nodesPerElement + 1):
            if kk != 0:
                element.append(int(item[kk]))
        
        elements.append(Element(element))

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    prescribedDisplacements = []
    for l in range(displacements):
        displacement = inputFile.readline().split(",")
        prescribedDisplacements.append([int(displacement[0]), [float(displacement[1]), float(displacement[2])]])

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    prescribedForces = []
    for l in range(forces):
        force = inputFile.readline().split(",")
        prescribedForces.append([int(force[0]), [float(force[1]), float(force[2])]])

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    material = inputFile.readline().split(",")
    material = [float(material[0]), float(material[1])]

    inputFile.close()

    end = timer()
    print("  1 - Leitura de dados: ", "%.5f" % (end - start), " segundos.\n")

    return prescribedDisplacements, prescribedForces, material, geometricNodes, internalPoints, elements

# prescribedDisplacements, prescribedForces, material, geometricNodes, internalPoints, elements = readInputFile("src/ex1_inputFileEP.txt")