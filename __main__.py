from sympy.integrals.quadrature import gauss_legendre

from src.utils.shapeFunctions import *
from src.utils.diffShapeFunction import *

def elasticityProblemBEM():
    # Definição do número de pontos de integraçãoque serão utilizados
    integrationPoints, weights = gauss_legendre(12, 10)

    return

elasticityProblemBEM()