import math
import matplotlib.pyplot as plt
import numpy as np

def generateHollowCylinder():
    elementsOrder = 1
    elementsOnLine = 3

    elementSize = 10

    l1Size = 100
    l2Size = 400
    l3Size = 400
    l4Size = 100
    l5Size = 500
    l6Size = 500

    # NODES
    nodes = []

    
    # l1
    for i in range(elementSize + 1):
        node = [0, 0]

        node[0] = 
        node[1] = height

        nodes.append(node)

    # VerticalLine 
    VLLelementsSize =  height / elementsOnLine
    for i in range((elementsOnLine ) * elementsOrder + 1):
        node = [0, 0]

        node[1] = height - (i * VLLelementsSize) / elementsOrder

        nodes.append(node)
        
    plotX = []
    plotY = []

    for i in range(len(nodes)):
        plotX.append(nodes[i][0])
        plotY.append(nodes[i][1])
    

    plt.scatter(plotX, plotY, label= "stars", color= "green", marker= "*", s=30)
    plt.savefig("teste")


    generatedFile = open('squareplate.txt', "a")

    # nodes
    for i in range(len(nodes)):
        generatedFile.write('%.0f' % i + ', ' + '%.3f' % nodes[i][0] + ', ' + '%.3f' % nodes[i][1] + '\n')

    # connectividade
    for i in range(elementsOnLine * (elementsOrder + 1)):
        dummy = 0
        if i > 2:
            dummy = 1
        if i > 5:
            dummy = 2
        if i > 8:
            dummy = 3
        generatedFile.write('%.0f' % i + ', ' + '%.0f' % (3*i + dummy + 0) + ', ' + '%.0f' % (3*i + dummy + 1) + ', ' + '%.0f' % (3*i + dummy + 2) + ', ' + '%.0f' % (3*i + dummy + 3) + '\n')

    # deslocamento
    for i in range(elementsOnLine * elementsOrder + 1):
        generatedFile.write('%.0f' % (30 + i) + ', 0, 0\n')

    # força
    for i in range(elementsOnLine * elementsOrder + 1):
        generatedFile.write('%.0f' % (10 + i) + ', 1, 0\n')

    generatedFile.close()



generateHollowCylinder()

