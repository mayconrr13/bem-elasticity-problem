import numpy as np
from src.shapeFunctions import *

def getIntegrationPointCoordinates(points: float, elementNodes: list, dimensionlessPoints: list):
        elementPointsCoordinates = []
        xComponent = 0
        yComponent = 0

        for j in range(len(elementNodes)):
            xComponent += getShapeFunctionValueOnNode(points, j, dimensionlessPoints) * elementNodes[j][0]
            yComponent += getShapeFunctionValueOnNode(points, j, dimensionlessPoints) * elementNodes[j][1]

        elementPointsCoordinates = [xComponent, yComponent]

        elementPointsCoordinates = np.array(elementPointsCoordinates)

        return elementPointsCoordinates