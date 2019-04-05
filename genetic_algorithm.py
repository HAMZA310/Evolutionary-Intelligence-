
from bitstring import *
import random
from collections import Counter
from itertools import chain

# global vars
TERMINATE = False
SHOW_FITNESS = False
COUNT = 0
    
def main():
    Holy_BitString = '01010101011001101101010111011001100101010101100101010101'
    Holy_BitString = ConstBitStream(bin = Holy_BitString)
    POPULATION_SIZE = 64
    CHROMOSOME_LENGTH = 56
    _show_fitness = input('Do you want to see Fitness of each candidate solution? yes/no ')
    global SHOW_FITNESS
    if _show_fitness in ['yes', 'y', '1', 'Y']:
        SHOW_FITNESS = True
    genetic_alg(Holy_BitString, POPULATION_SIZE, CHROMOSOME_LENGTH)
    print(COUNT)

def genetic_alg(Holy_BitString, POPULATION_SIZE, CHROMOSOME_LENGTH):
    solutions = generate_candidate_sols(POPULATION_SIZE)
    parents = evaluate_candidates(solutions, Holy_BitString, POPULATION_SIZE)

    
    while(not TERMINATE):
        pairs_of_parents = select_parents(parents)
        recombinded_parents = list(chain(*map(lambda pair: recombine_pairs_of_parents(pair[0], pair[1], \
                                                                                         CHROMOSOME_LENGTH),                                                                                                    pairs_of_parents))) # chain does: convert into list of offspring not of tuples
        mutated_offspring = list(map(lambda offspring: mutate_offspring(offspring, POPULATION_SIZE, \
                                                                            CHROMOSOME_LENGTH), recombinded_parents))
        parents = evaluate_candidates(mutated_offspring, Holy_BitString, POPULATION_SIZE) # new parents (offspring)

def random_num():
    random.seed()
    return random.randrange(2**14) ## for fitting in 14 bits. 
    
def generate_candidate_sols(n): 
    '''
        retrun a list of n candidate solutions (each a 56 bit string)
        functionality: generate 4 bitstrings of length 14 bits and then merge them to get 56 bit Strings
        
    '''
    if (n > 0):
        bslist = [ConstBitStream(uint=random_num(), length=14) for n in range(4)] 
        candidate_sol = ConstBitStream('0').join(bslist)
        return [candidate_sol] + generate_candidate_sols(n - 1)
                
    return []

 # Used in finding fiteness of a solution
def get_matching_bit_pairs(_Holy_BitString, _candidate_sol):
    _Holy_BitString.pos = 0
    _candidate_sol.pos = 0
    
    matching_bit_pairs = 0
    global TERMINATE, COUNT
    try:
        if not TERMINATE:
            while (_Holy_BitString.read(2).bin == _candidate_sol.read(2).bin):
                matching_bit_pairs = matching_bit_pairs + 1
            if (SHOW_FITNESS):
                print('Fitness: ', round((matching_bit_pairs)/28*100, 2), '%')
            COUNT = COUNT + 1
    except:
        print('Optimum Path Discovered: \n', _candidate_sol.bin)
        print('Fitness: 100 %')
        TERMINATE = True
    return matching_bit_pairs


def expected_counts(cur_sol_matching_bit_pairs, prev_gen_count_of_cur_sol, \
                            all_sols_matching_bit_pairs_count, POPULATION_SIZE):
    """
    1. all_sols_matching_bit_pairs :: int -> is sum of matching bit pairs of all solutions.
    2. Add-one smoothing is assumed. 
    3. n(i, t+1) = n(i, t) * f(i, t) / f'avg(t)
    4. returns: a float, take ceiling or floor later. 
    
    """
    total_possible_matching_pairs = 28
    f_sol = cur_sol_matching_bit_pairs / total_possible_matching_pairs # 
    f_avg = (all_sols_matching_bit_pairs_count / total_possible_matching_pairs) / POPULATION_SIZE
    _exp_counts = (f_sol / f_avg) * prev_gen_count_of_cur_sol
    
    return _exp_counts

def append_copies(fittest_parents, POPULATION_SIZE, candidates):
    parents = []
    for elm in fittest_parents:
        for i in range(elm[1]):
            if len(parents) >= POPULATION_SIZE:
                break
            parents.append(candidates[elm[0]])
    return parents

def evaluate_candidates(candidates, Holy_BitString, POPULATION_SIZE): 
    '''
        Additional Functions Used: 'round' converts to nearest int. 
        In 'cands_and_copies_count' Index 0 corresponds to 1st sol, 1 to 2nd, so on and so forth. 
        Optional functionality Not implemented: if more than half the pop is 0 
                                                in all_sols_matching_bit_pairs, add-one smooth. 
        
    '''
    all_cands_matching_bit_pairs = list(map(lambda cand: get_matching_bit_pairs(Holy_BitString, cand), candidates))
    if TERMINATE:
        return
    all_sols_matching_bit_pairs_count = sum(all_cands_matching_bit_pairs)
    all_expected_counts = list(map(lambda bit_pairs: round(expected_counts(bit_pairs, 1, \
                                all_sols_matching_bit_pairs_count, POPULATION_SIZE)),all_cands_matching_bit_pairs)) 
    cands_and_copies_count = Counter(dict(zip([x for x in range(0,POPULATION_SIZE)], all_expected_counts))) 
    fittest_parents = cands_and_copies_count.most_common(POPULATION_SIZE)
    parents = append_copies(fittest_parents, POPULATION_SIZE, candidates)
    return parents

def select_parents(parents):
    return tuple(zip(parents[0::2], parents[1::2])) # pairs of consecutive parents 

def recombine_pairs_of_parents(p1, p2, CHROMOSOME_LENGTH):
    """
        args: p1, and p2 are 56 Bit strings 
        split at .6-.9 of 56 bits (CHROMOSOME_LENGTH). i.e. between 31-50 bits
        
    """
    split_point = random.randrange((int) (CHROMOSOME_LENGTH * .6), (int)(CHROMOSOME_LENGTH*.9))

    p1 = BitStream(p1)
    p2 = BitStream(p2)
    p1[:split_point], p2[:split_point] = p2[:split_point], p1[:split_point]
    
    return p1, p2

def mutate_offspring(p, POPULATION_SIZE, CHROMOSOME_LENGTH):
   p = BitStream(p)
   if len(p) > 0:
       r = random.random()
       if (r >= 1/CHROMOSOME_LENGTH and r <= 1/POPULATION_SIZE) or (r <= 1/CHROMOSOME_LENGTH and r >= 1/POPULATION_SIZE):  # either number could be bigger. 
           p[:1] = BitStream(bin=str((int)(p[:1].bin) ^ 1))
       return p[:1] + mutate_offspring(p[1:],  POPULATION_SIZE, CHROMOSOME_LENGTH)
   return []


if __name__ == '__main__':
    main()