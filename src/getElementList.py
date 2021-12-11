import numpy as np
from Element import Element

def getElementsList(elements: list): 
    elementsList = np.zeros(len(elements), dtype=list)
    
    for el in range(len(elements)):
        elementsList[el] = Element(elements[el])

    return elementsList