import random
from bitarray import bitarray  # /!\ Need to install bitarray /!\
import binascii
import pickle 
import os
import time as t

import zinc_grammar 
import score_util
import cfg_util
import sascorer
import optimize_J as opt

from rdkit.Chem import Draw
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import AllChem as Chem


# chemGE version
def generate_bit_vector(P):  # Generate a new bit vector from the probability vector
    
    bit_vector = []
    for j in range(0, len(P)):
        if random.random() < P[j]:
            bit_vector.append(1)

        else:
            bit_vector.append(0)
    
    return bitarray(bit_vector)


# chemGE version
def BITtoGene(bit_vector):  # Takes a bitarray object
    # The bit_vector is composed of multiple bytes aligned one after the other
    # 8 bits = 1 byte 
    # 1 byte is translated into 1 integer

    #  /!\ ONLY WORKS WITH INTEGER /!\

    byte_list = [bit_vector[x:x+8] for x in range(0,len(bit_vector),8)]  # Splits the bitvector in a list of byte

    gene = []
    
    for byte in byte_list: # Translates the byte list in list of int
        out = 0
        for bit in byte:
            out = (out << 1) | bit  # Magic : byte to int
        gene.append(out)

    return gene


# chemGE version
def evaluate(bit_vector):  # Takes a bitarray object
    # gene is a list of int, bit_vector a bitarray
    gene = BITtoGene(bit_vector) 
    smile = opt.canonicalize(cfg_util.decode(opt.GenetoCFG(gene)))  # Transform the gene into a smile
    score = score_util.calc_score(smile)  # Calculate the J score
    return score


def convergence(P):
    # Take the vector of probability and check if all the elements are
    # inferior to 0.05 or superior to 0.95
    boolean = all((i >0.05 and not i<0.95) or (not i>0.05 and i <0.95) for i in P)
    # If TRUE then it stop the while loop
    if boolean == True:
        print('CONVERGENCE')
    return boolean 
 

def save_log(population):
    save = input("Save logs and image of final population ? (press 'y' or 'n') : ")
    if save == 'n':
        pass

    else:
        directory = input("Please input log file name (or directory) : ")
        
        # Creating a folder for this log

        os.system('mkdir '+ directory)
        file_name = directory

        # Stocking the final population in a pickle object

        f = open(directory + '/' + file_name + ".p",'wb')
        
        # Remove double and non valid smile from the list befor stocking it

        ms = []
        smile_list = []
        for bit_vector in population:
            gene = BITtoGene(bit_vector) 
            smile = opt.canonicalize(cfg_util.decode(opt.GenetoCFG(gene)))
            
            if smile != '' and smile != None and smile not in smile_list:            
                if MolFromSmiles(smile) != None:
                    smile_list.append(smile)
                    ms.append(MolFromSmiles(smile))
        pickle.dump(smile_list,f)
        f.close()

        # Stocking the random.seed of this experiement in a text file 

        f = open(directory + '/' + 'seed.txt','w')
        f.write(str(time)+'\n')
        #Stocking the final population and their score in the same file
        f.write('smile' +'\t'+ 'score'+'\n')
        for smile in smile_list :
            score = score_util.calc_score(smile)
            f.write(smile + '\t' + str(score)+'\n')
        f.close()
    
    # Saving population Image

    if save == 'n':
        pass
    
    else:
        for i in range(len(ms)):
            Draw.MolToFile(ms[i], directory+'/'+ str(i) + '.png' , size=(120,120))
        os.system('montage ' +directory+'/*.png ' +directory+'/final.png')  # Execute this command in the shell. Put all images of the molecules in a unique image
        # Need image magick to work...


def main():

    global time
    time = t.time()
    print(time)
    random.seed(time)

    max_generation = 1000000
    max_time = 8*3600  # 8 hours
    population_size = 100
    bit_vector_size = 2400  # Maximum length of vectors in the population (should be a multiple of 8)

    P = [0.5 for _ in range(0,bit_vector_size)]  # Probability vector
    LR = 0.1  # Learning Rate (typically 0.1–0.4)
    MS = 0.05  # Degree of mutation (typical value is 0.05)
    Pr_mutation = 0.08  # Probability of mutation (typically 0.02)

    mu = 2  # Number of vector used to make P evolve

    k = 0
    duration = 0
    converge = False
    best_fitness = -1e11
    best_bit_vector = None
    while converge is not True and duration < max_time:  # k < max_generation or duration < max_time depending on what you want
        
        population = []
        score_smile = []
        best_bit_vector = None
        best_fitness = -1e10
        for i in range(0, population_size):

            bit_vector = generate_bit_vector(P)  # Create a new vector which represents an individual
            population.append(bit_vector)
            fitness = evaluate(population[i])  # Evaluate the fitness of the new vectorsr
            if fitness > -1e10:
                score_smile.append([fitness,bit_vector])
            #print(fitness)
            #gene = BITtoGene(bit_vector) 
            #smile = opt.canonicalize(cfg_util.decode(opt.GenetoCFG(gene)))
            #print(smile)
            if fitness > best_fitness:  # /!\ '<' and '>'

                best_fitness = fitness  # Update the best individual (i.e. max fitness)
                best_bit_vector = bit_vector
                
                gene = BITtoGene(best_bit_vector) 
                smile = opt.canonicalize(cfg_util.decode(opt.GenetoCFG(gene)))
                print(best_fitness,'=',smile)
        try :
            # Evolution
            for j in range(0, len(P)):

                #P[j] = P[j]*(1 - LR) + int(best_bit_vector[j])*LR  # Update the probability vector with the best indiv

                score_smile = sorted(score_smile, key = lambda x: x[0], reverse = False)  # The best smile is a the end of the list 
                print(score_smile)
                if len(score_smile) < mu:
                    N = len(score_smile)
                else:
                    N = mu
                X = 0
                for i in range(N):
                    X += (i+1)*(score_smile[i][1][j]-P[j])
                
                P[j] = P[j] + LR/(P[j]*(1-P[j]))*X  # Information Geometric implementation
                
            # Mutation
            for j in range(0,len(P)):

                if random.random() < Pr_mutation:

                    P[j] = P[j]*(1 - MS) + random.randint(0, 1)*MS
        except :
            print('No valid SMILE generated : pass')
        converge = convergence(P)        
        k += 1
        duration = t.time()-time
        print(k,' time : ',duration, ' s')
        
    gene = BITtoGene(best_bit_vector) 
    smile = opt.canonicalize(cfg_util.decode(opt.GenetoCFG(gene)))
    print(smile)

    
    save_log(population)

    return best_bit_vector


if __name__ == '__main__':
    main()  


# O=C(Nc1ccc(CBr)cc1)OSc1ccc(-c2ccccc2)c(-c2ccccc2)c1