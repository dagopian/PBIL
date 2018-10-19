import random
from bitarray import bitarray  # /!\ Need to install bitarray /!\
import binascii

import zinc_grammar 
import score_util
import cfg_util
import sascorer
import optimize_J as opt


def generate_bit_vector(P):  # Generate a new bit vector from the probability vector

    bit_vector = []
    for j in range(0, len(P)):
        if random.random() < P[j]:
            bit_vector.append(1)
    
        else:
            bit_vector.append(0)
    
    return bitarray(bit_vector)




def evaluate_test(bit_vector):
    x = int(bit_vector[0:5].to01(), 2) #Convert the 5 first digit of the bitarray in int
    y = int(bit_vector[5:10].to01(), 2) #Convert the 5 last digit of the bitarray in int
    #print(x,y)
    log_x.append(x)
    log_y.append(y)
    #return x**2 + y**2  
    return x**4 - x**3 - 20*x**2 + x + 1 + y**4 - y**3 - 20*y**2 + y + 1

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


## Taken from stackoverflow forum ##
def text_to_bits(text, encoding='utf-8', errors='surrogatepass'):
    bits = bin(int(binascii.hexlify(text.encode(encoding, errors)), 16))[2:]
    return bits.zfill(8 * ((len(bits) + 7) // 8))

def bits_to_text(bits, encoding='utf-8', errors='surrogatepass'):
    n = int(bits, 2)
    return int2bytes(n).decode(encoding, errors)

def int2bytes(i):
    hex_string = '%x' % i
    n = len(hex_string)
    return binascii.unhexlify(hex_string.zfill(n + (n & 1)))
##                               ##

def BITtoGene2(bit_vector):

    gene = bits_to_text(bit_vector)

    return gene



def evaluate(bit_vector):  # Takes a bitarray object
    # gene is a list of int, bit_vector a bitarray
    gene = BITtoGene(bit_vector) 
    smile = opt.canonicalize(cfg_util.decode(opt.GenetoCFG(gene)))  # Transform the gene into a smile
    score = score_util.calc_score(smile)  # Calculate the J score
    return score

def main():

    global log_x     # Optional, to keep a trace of the evolution
    global log_y

    log_x = []
    log_y = []

    max_generation = 100000
    population_size = 100
    bit_vector_size = 160  # Maximum length of vectors in the population (should be a multiple of 8)

    P = [0.5 for _ in range(0,bit_vector_size)]  # Probability vector
    LR = 0.2  # Learning Rate (typically 0.1–0.4)
    MS = 0.05  # Degree of mutation (typical value is 0.05)
    Pr_mutation = 0.02  # Probability of mutation (typically 0.02)

    k = 0
    converge = False  # TO DO
    best_fitness = -10**11
    while converge == False and k < max_generation:
        
        population = []
        
        for i in range(0, population_size):

            bit_vector = generate_bit_vector(P)  # Create a new vector which represent an individual
            population.append(bit_vector)
            fitness = evaluate(population[i])  # Evaluate the fitness of the new vector
            print(fitness)
            if fitness > best_fitness:  # /!\ '<' and '>'

                best_fitness = fitness  # Update the best individual (i.e. max fitness)
                best_bit_vector = bit_vector

        for j in range(0, len(P)):

            P[j] = P[j]*(1 - LR) + best_bit_vector[j]*LR  # Update the probability vector with the best indiv

        for j in range(0,len(P)):

            if random.random() < Pr_mutation:

                P[j] = P[j]*(1 - MS) + random.randint(0, 1)*MS

        k += 1
        print(k)
        #print(best_bit_vector, int(best_bit_vector[0:5].to01(), 2), int(best_bit_vector[5:10].to01(), 2)) #Debugging to delete
    #   print('\n',log_x,'\n',log_y)    
    return best_bit_vector


if __name__ == '__main__':
    main()  


# O=C(Nc1ccc(CBr)cc1)OSc1ccc(-c2ccccc2)c(-c2ccccc2)c1