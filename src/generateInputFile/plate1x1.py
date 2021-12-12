import math
import matplotlib.pyplot as plt
import numpy as np

def generateHollowCylinder():
    elementsOrder = 3
    elementsOnLine = 3
    elementsOnCurve = 10

    width = 1 #mm
    height = 1 #mm

    Fy = 2
    E = 1 #N/mm2
    poisson = 0

    # NODES
    nodes = []

    # horizontalLine bottom
    HLBelementsSize =  width / elementsOnLine
    for i in range((elementsOnLine ) * elementsOrder + 1):
        node = [0, 0]

        node[0] = (i * HLBelementsSize) / elementsOrder

        nodes.append(node)

    # vertical line right
    VLRelementsSize =  height / elementsOnLine
    for i in range((elementsOnLine ) * elementsOrder + 1):
        node = [0, 0]

        node[0] = width
        node[1] = (i * VLRelementsSize) / elementsOrder

        nodes.append(node)

    # horizontalLine top
    HLTelementsSize =  width / elementsOnLine
    for i in range((elementsOnLine) * elementsOrder + 1):
        node = [0, 0]

        node[0] = width - (i * HLTelementsSize) / elementsOrder
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

