from sympy.integrals.quadrature import gauss_legendre

from src.utils.shapeFunctions import *
from src.utils.diffShapeFunction import *

def elasticityProblemBEM():
    # ler arquivo de entrada
    # GMESH???

    # criar malha de colocação

    # RESOLVER PROBLEMA
    # Definição do número de pontos de integraçãoque serão utilizados
    integrationPoints, weights = gauss_legendre(12, 10)

    # criar arquivo de saída
    # PARAVIEW

    return

elasticityProblemBEM()