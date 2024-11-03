from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy

poblacion = []
path = r"C:\Users\Fer_m\OneDrive\Documentos\UAdeC\8voSemestre\Administracion de Proyectos de Software\BFOA-main\multiFasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)
tempBacteria = bacteria(path)
original = bacteria(path)
globalNFE = 0

# Ajustar parámetros de atracción y repulsión
dAttr = 0.15  # incremento de atracción
wAttr = 0.25
hRep = dAttr
wRep = 10    # reducción de repulsión para mejorar convergencia

def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

def validaSecuencias(path, veryBest):
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

for i in range(numeroDeBacterias):
    poblacion.append(bacteria(path))

for _ in range(iteraciones):
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        bacteria.autoEvalua()

    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    
    # Solo incrementar NFE si hay mejora en fitness
    
    globalNFE += chemio.parcialNFE

    best = max(poblacion, key=lambda x: x.fitness)
    if veryBest is None or best.fitness > veryBest.fitness:
        clonaBest(veryBest, best)
        
    print("Interacción:", veryBest.interaction, "Fitness:", veryBest.fitness, "NFE:", globalNFE)

    if _ % 2 == 0:  # Inserta bacterias aleatorias cada dos iteraciones
        chemio.eliminarClonar(path, poblacion)
        chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)
    
    print("Población:", len(poblacion))

veryBest.showGenome()
validaSecuencias(path, veryBest)
