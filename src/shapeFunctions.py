def getShapeFunctionValueOnNode(ksi: float, node: int, adimensionalPoints: list):
    shapeFunctionValueOnPoint = 1

    for i in range(len(adimensionalPoints)):
        if i != node:
            shapeFunctionValueOnPoint *= (ksi - adimensionalPoints[i]) / (adimensionalPoints[node] - adimensionalPoints[i])

    return shapeFunctionValueOnPoint
