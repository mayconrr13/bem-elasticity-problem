import math
import matplotlib.pyplot as plt
import numpy as np

def generateHollowCylinder():
    elementsOrder = 1
    

    elementSize = 100

    l1Size = 150
    l2Size = 250
    l3Size = 250
    l4Size = 150
    l5Size = 100 * 2 ** 0.5

    print(l5Size)

    # NODES
    nodes = []
    connec = []
    acc = 0
    
    # l1
    for i in range(100 + 1):
        node = [0, 0]

        node[0] = 100 + (l1Size / 100) * i 
        node[1] = 0

        nodes.append(node)

        if i % 2 != 0:
            connec.append(str(i) + '\n')
            connec.append(str(i) + ',')
        else:
            connec.append(str(i) + ',')

    connect = []
    for j in range(100):
        connect.append([j, j + 1])
    for j in range(100):
        connect.append([100 + j + 1,100 + j + 2])
    for j in range(100):
        connect.append([200 + j + 2,200 + j + 3])
    for j in range(100):
        connect.append([300 + j + 3,300 + j + 4])
    for j in range(100):
        connect.append([400 + j + 4,400 + j + 5])
    # print(connect)

    
  
    acc = (len(connec))
    # l2
    for i in range(100 + 1):
        node = [0, 0]

        node[0] = 250 
        node[1] = (l2Size / 100) * i

        nodes.append(node)
        
        if i % 2 != 0:
            connec.append(str(acc + i) + '\n')
        else:
            connec.append(str(acc + i) + ',')
    
    print(connec)

    acc = (len(connec))
    # l3
    for i in range(100 + 1):
        node = [0, 0]

        node[0] = 250 - (l3Size / 100) * i 
        node[1] = 250

        nodes.append(node)

        if i % 2 != 0:
            connec.append(str(acc + i) + '\n')
        else:
            connec.append(str(acc + i) + ',')

    acc = (len(connec))
    # l4
    for i in range(100 + 1):
        node = [0, 0]

        node[0] = 0
        node[1] = 250 - (l4Size / 100) * i

        nodes.append(node)

        if i % 2 != 0:
            connec.append(str(acc + i) + '\n')
        else:
            connec.append(str(acc + i) + ',')

    acc = (len(connec))
    # l5
    for i in range(100 + 1):
        node = [0, 0]

        node[0] = (100 / 100) * i
        node[1] = 100 - (100 / 100) * i

        nodes.append(node)

        if i % 2 != 0:
            connec.append(str(acc + i) + '\n')
        else:
            connec.append(str(acc + i) + ',')
    list = " ".join(connec)
    
    for i in range(101): 
        print(str(300 + i + 3) + ', 0, 0')
    # for i in range(101): 
    #     print(str(101 + i) + ', 0, 100')
    # print(len(nodes))
    plotX = []
    plotY = []

    for i in range(len(nodes)):
        plotX.append(nodes[i][0])
        plotY.append(nodes[i][1])
    

    plt.scatter(plotX, plotY, label= "stars", color= "green", marker= "*", s=30)
    plt.savefig("teste")


    generatedFile = open('hollowPlate.txt', "a")
    for i in range(len(connect)):
        generatedFile.write('%.0f' % i + ', ' + '%.0f' % connect[i][0] + ', ' + '%.0f' % connect[i][1] + '\n')
        
    # nodes
    # for i in range(len(nodes)):
    #     generatedFile.write('%.0f' % i + ', ' + '%.9f' % nodes[i][0] + ', ' + '%.9f' % nodes[i][1] + '\n')

    # connectividade
    # for i in range(101):
    #     print([])
        
    # # deslocamento
    # for i in range(elementsOnLine * elementsOrder + 1):
    #     generatedFile.write('%.0f' % (30 + i) + ', 0, 0\n')

    # # força
    # for i in range(elementsOnLine * elementsOrder + 1):
    #     generatedFile.write('%.0f' % (10 + i) + ', 1, 0\n')

    generatedFile.close()



generateHollowCylinder()

