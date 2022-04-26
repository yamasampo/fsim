import math
from multiprocessing.sharedctypes import Value
from pandas import value_counts
from scipy.integrate import quad
from typing import Union

def time_to_fixation_haploid(
        init_freq:float, n:int, s:Union[int, float]):
    """ Returns expected time to fixation for mutations that ultimately go to 
    fixation. This function uses equation (17) of Kimura and Ohta (1968).

    Parameters
    ----------
    init_freq: float
        A mutant allele frequency in a population at the initial generation. '
    n: int
        Population size. Number of individuals of haploid organisms. 
    s: int or float
        Selection coefficient. Mutant allele has a selective advantage s 
        over wildtype allele. 

    Return
    ------
    expected time to fixation: float

    Reference
    ---------
    Kimura M., and T. Ohta, 1969 The Average Number of Generations until 
    Fixation of a Mutant Gene in a Finite Population. Genetics 61: 763-771. 
    https://doi.org/10.1093/genetics/61.3.763
    """
    assert s != 0
    
    # Adjust parameters to fit the formula.
    # Kimura and Ohta (1969) considered the mutant gene A1 has the selective 
    # advantage s/2 over its allele A1. But the input argument s itself is 
    # defined as the selective advantage. 
    s = 2 * s

    # In Kimura and Ohta (1969) model, n is the number of individuals of 
    # diploid organisms. 
    # NOTE: Here, we assume that unbiased gamete production and no dominance 
    # for the diploid organism.
    n = n / 2

    # Calculate population size scaled selection coefficient (selection intensity)
    ns = n * s

    # Calculate the first term of equation (17)
    first_term = J1(init_freq, s, ns)

    # Calculate the second term of equation (17)
    up = u(init_freq, ns)
    second_term = (1-up)/up * J2(init_freq, s, ns)

    return first_term + second_term

def u(init_freq, ns):
    denom = 1 - math.exp(-2 * ns)
    nume = 1 - math.exp(-2 * ns * init_freq)
    
    return nume / denom
    
def J1(init_freq, s, ns):
    def func1(x, ns):
        return (math.exp(2*ns*x)-1)*(math.exp(-2*ns*x)-math.exp(-2*ns)) / x / (1-x)
    
    term1 = 2 / s / (1-math.exp(-2*ns))
    term2 = quad(func1, init_freq, 1, args=(ns))
    
    return term1 * term2[0]
        
def J2(init_freq, s, ns):
    def func2(x, ns):
        return (math.exp(2*ns*x)-1)*(1-math.exp(-2*ns*x)) / x / (1-x)
    
    term1 = 2 / s / (1-math.exp(-2*ns))
    term2 = quad(func2, 0, init_freq, args=(ns))
    
    return term1 * term2[0]

def time_to_fixation_diploid(init_freq, ns, s):
    """ From equation (17) of Kimura and Ohta (1968) 
    """
    assert ns != 0
    up = u(init_freq, ns)
    return J1(init_freq, s, ns) + (1-up)/up * J2(init_freq, s, ns)

def fixation_probability_haploid(init_freq, ns):
    if ns == 0:
        return init_freq
    
    denom = 1 - math.exp(-2 * ns)
    nume = 1 - math.exp(-2 * ns * init_freq)
    
    return nume / denom

def prob_of_obtaining_given_number_of_neutral_mutants_next_gen(
        pop_size:int,
        parental_mut_num:int,
        offspring_mut_num:int
        ):
    """ Returns a probability of obtaining a given number of mutant alleles at 
    the next generation, given a number of mutant alleles in a parental 
    population. This function considers haploid organism (i.e., population size 
    is equal to a total number of alleles in the population).
    
    Assumption
    ----------
    1. No mutation: 
        A new mutation does not occur.
    2. Random mating: 
        Each individual has an equal probability of becoming a parent. 
    3. Perfect reproduction: 
        All the parents successfully produce offsprings.
    4. Constant population size: 
        The population size does not change over time. 
        parental pop size = offspring pop size = pop_size
    5. Mutually exclusive configurations
        Observing different types of configurations (order of mutant and 
        wildtype alleles) do not occur simultaneously. 
    
    Parameters
    ----------
    pop_size: int
        A number of individuals in a population.
    parental_mut_num: int
        A number of mutant alleles in the parent population.
    offspring_mut_num: int
        A number of mutant alleles in the next (offspring) population.

    """
    # Probability of "mating" (haveing an opportunity to become a parent) 
    # for a mutant individual at t generation. 
    #   Assumption 1: no mutation
    #   Assumption 2: random mating (1 / pop_size x parental_mut_num)
    pt_mating_mut = parental_mut_num / pop_size

    # Probability of reproduction (making an offspring) for the parent
    #   Assumption 3: perfect reproduction
    pt_mating_success_mut = pt_mating_mut
    pt_mating_success_wild = 1 - pt_mating_success_mut
    
    # Probability that mutant and wildtype individuals made offsprings in a 
    # given order (configuration)
    #   Assumption 2: Random mating
    #   Assumption 3: Perfect reproduction
    #   Assumption 4: Constant population size
    pt_observing_one_config = \
        pt_mating_success_mut ** offspring_mut_num \
        * pt_mating_success_wild ** (pop_size - offspring_mut_num)
    
    # A total number of configurations
    #   Assumption 4: Constant population size
    #   Assumption 5: Mutually exclusive configurations
    total_configs = number_of_combinations(pop_size, offspring_mut_num)
    return pt_observing_one_config * total_configs

def number_of_combinations(n, r):
    if r == 0:
        return 1
    elif r == n:
        return 1

    return math.factorial(n) / (math.factorial(n - r) * math.factorial(r))

