import math
import random

from bacteria import bacteria


class chemiotaxis():
    def __init__(self):
        self.parcialNFE = 0
        self.best_interaction = float('-inf')
    
    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        total = 0.0
        # Optimización: considera más bacterias en el cálculo
        sorted_poblacion = sorted(poblacion, 
                                key=lambda x: abs(x.blosumScore - bacteria.blosumScore))[:10]  # considerar las 10 más cercanas

        for other in sorted_poblacion:
            if other is bacteria:
                continue
            diff = (bacteria.blosumScore - other.blosumScore) ** 2.0
            total += d * math.exp(w * diff)
        
        # Podrías promediar las interacciones
        return total / len(sorted_poblacion) if sorted_poblacion else 0

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)
        repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)
        return attract + repel
    
    
    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction
        self.best_interaction = max(self.best_interaction, bacteria.interaction)
      

    
    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
        self.parcialNFE = 0
        for bacteria in poblacion:
            self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
            self.parcialNFE += bacteria.NFE
            bacteria.NFE = 0
        

    def eliminarClonar(self, path, poblacion):
        """Versión mejorada de eliminación y clonación"""
        # Ordenamos por fitness
        poblacion.sort(key=lambda x: x.fitness)
        
        # Mantener un mínimo de diversidad
        min_diversidad = 3  # número de bacterias a mantener para diversidad
        unique_bacteria = list({bacteria: None for bacteria in poblacion[:min_diversidad]}.keys())
        
        # Eliminar bacterias inferiores, manteniendo al menos min_diversidad
        del poblacion[min_diversidad:]
        
        # Clonación de las mejores bacterias
        clones = self.clonacion(path, poblacion)
        poblacion.extend(clones)
        
        # Asegúrate de no superar el tamaño de la población original
        if len(poblacion) > len(unique_bacteria) + len(clones):
            poblacion = unique_bacteria + poblacion[:len(poblacion) - len(unique_bacteria) - len(clones)]
        
    def clonacion(self, path, poblacion):
        poblacionClones = []
        best = max(poblacion, key=lambda x: x.fitness)
        
        for bacteria in poblacion:
            newBacteria = bacteria.clonar(path)
            # Mutación adaptativa basada en la diferencia de fitness
            fitness_diff = max(1, int((best.fitness - bacteria.fitness)/5))
            mutacion = min(5, max(1, fitness_diff))  # Entre 1 y 5
            newBacteria.tumboNado(mutacion)
            newBacteria.autoEvalua()
            poblacionClones.append(newBacteria)
        
        return poblacionClones

    
         
    def randomBacteria(self, path):
        bact = bacteria(path)
        # Aplicamos un número aleatorio de tumbles para mayor diversidad
        bact.tumboNado(random.randint(1, 3))
        bact.autoEvalua()
        return bact
   
    def insertRamdomBacterias(self, path, num, poblacion):
        """Versión mejorada de inserción de bacterias aleatorias"""
        for _ in range(num):
            new_bact = self.randomBacteria(path)
            # Aseguramos que la nueva bacteria sea diferente a las existentes
            while any(self.is_similar(new_bact, b) for b in poblacion):
                new_bact = self.randomBacteria(path)
            
            poblacion.append(new_bact)
            # Eliminamos la peor bacteria
            poblacion.sort(key=lambda x: x.fitness)
            del poblacion[0]
    
    def is_similar(self, bact1, bact2, threshold=0.9):
        """Comprueba si dos bacterias son muy similares"""
        if len(bact1.matrix.seqs) != len(bact2.matrix.seqs):
            return False
        
        similarity = 0
        total = 0
        for seq1, seq2 in zip(bact1.matrix.seqs, bact2.matrix.seqs):
            for c1, c2 in zip(seq1, seq2):
                if c1 == c2:
                    similarity += 1
                total += 1
        
        return (similarity / total) > threshold
         
    

# Path: BFOA_MSAv2/evaluadorBlosum.py
