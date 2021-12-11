def getDiffShapeFunction(ksi: float, node: int, adimensionalPoints: list):
    diffShapeFunctionValueOnPoint = 0

    for i in range(len(adimensionalPoints)):
        if i != node:
            partial = 1

            for j in range(len(adimensionalPoints)):
                if i == j:
                    partial *= (1) / (adimensionalPoints[node] - adimensionalPoints[j])

                elif i != j and j != node:
                    partial *= (ksi - adimensionalPoints[j]) / (adimensionalPoints[node] - adimensionalPoints[j])

            diffShapeFunctionValueOnPoint += partial

    return diffShapeFunctionValueOnPoint
    