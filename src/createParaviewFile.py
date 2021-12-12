from timeit import default_timer as timer

def createParaviewFile(boundaryDisplacements, boundaryForces, internalDisplacements, internalStress,auxiliaryMesh, internalPoints, elements):
    start = timer()

    outputFile = open('results.vtu', 'a')

    outputFile.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid">\n  <UnstructuredGrid>\n  <Piece NumberOfPoints="' + str(len(auxiliaryMesh) + len(internalPoints)) + '" NumberOfCells="' + str(len(elements) + len(internalPoints)) + '">\n')
    
    # nós
    outputFile.write('    <Points>\n      <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
    for i in range(len(auxiliaryMesh)):
        x = auxiliaryMesh[i][0] + boundaryDisplacements[2 * i]
        y = auxiliaryMesh[i][1] + boundaryDisplacements[2 * i + 1]
        outputFile.write(str("%.5f" % x) + ' ' + str("%.5f" % y) + ' 0\n')

    for j in range(len(internalPoints)):
        x = internalPoints[j][0] + internalDisplacements[2 * j]
        y = internalPoints[j][1] + internalDisplacements[2 * j + 1]
        outputFile.write(str("%.5f" % x) + ' ' + str("%.5f" % y) + ' 0\n')
    outputFile.write('      </DataArray>\n    </Points>\n')

    # elementos
    outputFile.write('    <Cells>\n      <DataArray type="Int32" Name="connectivity" format="ascii">\n')
    for i in range(len(elements)):
        nodeList = elements[i].nodeList                
        outputFile.write(str(nodeList[0]) + ' ' + str(nodeList[1]) + '\n')
    
    for j in range(len(internalPoints)):              
        outputFile.write(str(len(auxiliaryMesh) + j) + '\n')
    outputFile.write('      </DataArray>\n')
    
    outputFile.write('      <DataArray type="Int32" Name="offsets" format="ascii">\n')
    for i in range(len(elements)):                
        outputFile.write(str((i + 1) * 2) + '\n')

    for j in range(len(internalPoints)):                
        outputFile.write(str(len(auxiliaryMesh) + 1 + j) + '\n')
    outputFile.write('      </DataArray>\n')
    
    outputFile.write('      <DataArray type="UInt8" Name="types" format="ascii">\n')
    for i in range(len(elements)):                
        outputFile.write('68\n')

    for i in range(len(internalPoints)):                
        outputFile.write('1\n')
    outputFile.write('      </DataArray>\n')

    outputFile.write('    </Cells>\n')
    outputFile.write('    <PointData>\n')

    # deslocamentos
    outputFile.write('      <DataArray type="Float64" NumberOfComponents="2" Name="Displacements" format="ascii">\n')    
    for i in range(len(auxiliaryMesh)):  
        x = boundaryDisplacements[2 * i]
        y = boundaryDisplacements[2 * i + 1]              
        outputFile.write(str("%.5f" % x) + ' ' + str("%.5f" % y) + '\n')

    for j in range(len(internalPoints)):  
        x = internalDisplacements[2 * j]
        y = internalDisplacements[2 * j + 1]              
        outputFile.write(str("%.5f" % x) + ' ' + str("%.5f" % y) + '\n')
    outputFile.write('      </DataArray>\n')

    # tensoes
    outputFile.write('      <DataArray type="Float64" NumberOfComponents="4" Name="Stress" format="ascii">\n')    
    for i in range(len(auxiliaryMesh)):          
        outputFile.write('0 0 0 0\n')

    for j in range(len(internalPoints)):  
        sigxx = internalStress[4 * j + 0]             
        sigxy = internalStress[4 * j + 1]             
        sigyx = internalStress[4 * j + 2]             
        sigyy = internalStress[4 * j + 3]             
        outputFile.write(str("%.5f" % sigxx) + ' ' + str("%.5f" % sigxy) + ' '  + str("%.5f" % sigyx) + ' ' + str("%.5f" % sigyy) +'\n')
    outputFile.write('      </DataArray>\n')

    outputFile.write('    </PointData>\n')
    outputFile.write('  </Piece>\n  </UnstructuredGrid>\n</VTKFile>')
    
    outputFile.close()

    end = timer()
    print("  5 - Arquivo de saída para Paraview: ", "%.5f" % (end - start), " segundos.\n")
    return

#     write(11,'(a)') '      <DataArray type="Float64" NumberOfComponents="1" Name="SigmaX" format="ascii">'
#     do i = 1,nnos
#         write(11,'(es12.5)') tensao(i,1)
#     end do
#     write(11,'(a)') '      </DataArray>'
 