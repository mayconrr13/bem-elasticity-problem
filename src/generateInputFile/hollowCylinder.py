import math
import matplotlib.pyplot as plt
import numpy as np

def generateHollowCylinder():
    elementsOrder = 3
    elementsOnLine = 2
    elementsOnCurve = 2

    internalRadius = 20 #mm
    width = 50 #mm
    height = 50 #mm

    Fx = 100
    Fy = 100
    E = 200000 #N/mm2
    poisson = 0.3

    # NODES
    nodes = []

    # horizontalLine bottom
    HLBelementsSize =  (width - internalRadius) / elementsOnLine
    for i in range((elementsOnLine ) * elementsOrder + 1):
        node = [0, 0]

        node[0] = internalRadius + (i * HLBelementsSize) / elementsOrder

        nodes.append(node)

    # vertical line right
    VLRelementsSize =  height / elementsOnLine
    for i in range((elementsOnLine ) * elementsOrder + 1):
        node = [0, 0]

        node[0] = width
        node[1] = (i * VLRelementsSize) / elementsOrder

        nodes.append(node)

    # horizontalLine bottom
    HLTelementsSize =  width / elementsOnLine
    for i in range((elementsOnLine) * elementsOrder + 1):
        node = [0, 0]

        node[0] = width - (i * HLTelementsSize) / elementsOrder
        node[1] = height

        nodes.append(node)

    # VerticalLine 
    VLLelementsSize =  (height - internalRadius) / elementsOnLine
    for i in range((elementsOnLine ) * elementsOrder + 1):
        node = [0, 0]

        node[1] = height - (i * VLLelementsSize) / elementsOrder

        nodes.append(node)

    # raio interno
    # elementArc = (math.pi * 0.5) / (elementsOrder * elementsOnCurve)
    # for i in range((elementsOnCurve) * elementsOrder + 1):
    #     #elemento
    #     if i == 0 or i == elementsOnCurve:
    #         nodei = [0, 0]
    #         nodef = [0, 0]

    #         nodei[0] = internalRadius * math.sin((elementArc * i))
    #         nodei[1] = internalRadius * math.cos((elementArc * i))

    #         nodef[0] = internalRadius * math.sin((elementArc * i))
    #         nodef[1] = internalRadius * math.cos((elementArc * i))

    #         nodes.append(nodei)
    #         nodes.append(nodef)
    
    #     else: 
    #         node = [0, 0]

    #         node[0] = internalRadius * math.sin(elementArc * i)
    #         node[1] = internalRadius * math.cos(elementArc * i)

    #         nodes.append(node)
        
    # plotX = []
    # plotY = []

    # for i in range(len(nodes)):
    #     plotX.append(nodes[i][0])
    #     plotY.append(nodes[i][1])
    

    # plt.scatter(plotX, plotY, label= "stars", color= "green", marker= "*", s=30)
    # plt.savefig("teste")

    print(nodes[0][0])

    generatedFile = open('hollowPlate.txt', "a")
    for i in range(len(nodes)):
        generatedFile.write('%.0f' % i + ', ' + '%.3f' % nodes[i][0] + ', ' + '%.3f' % nodes[i][1] + '\n')
    generatedFile.close()



generateHollowCylinder()

