def lineToIgnore(inputFile: str):
    inputFile.readline()
    return

def readInputFile(file: str):
    inputFile = open(file, "r")

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    numberOfExternalNodes = int(inputFile.readline().split(":")[1])
    numberOfInternalNodes = int(inputFile.readline().split(":")[1])   
    nodesPerElement = int(inputFile.readline().split(":")[1]) 
    numberOfElements = int(inputFile.readline().split(":")[1])
    prescribedValues = int(inputFile.readline().split(":")[1])

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
        
        elements.append(element)

    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)
    lineToIgnore(inputFile)

    u= []
    q= []
    for l in range(prescribedValues):
        potentialFlow = inputFile.readline().split(",")

        if int(potentialFlow[1]) == 1:
            q.append([int(potentialFlow[0]), [float(potentialFlow[2]), float(potentialFlow[3])]])

        elif int(potentialFlow[1]) == 0:
            u.append([int(potentialFlow[0]), [float(potentialFlow[2]), float(potentialFlow[3])]])

    inputFile.close()

    return u, q, geometricNodes, internalPoints, elements

# u, q, geometricNodes, internalPoints, elements = readInputFile("ex1_inputFile.txt")