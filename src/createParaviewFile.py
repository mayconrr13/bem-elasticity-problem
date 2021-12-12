from timeit import default_timer as timer

def createParaviewFile(boundaryDisplacements, boundaryForces, internalDisplacements, internalStress):
    start = timer()

    # ...

    end = timer()
    print("  5 - Arquivo de saída para Paraview: ", "%.5f" % (end - start), " segundos.\n")
    return