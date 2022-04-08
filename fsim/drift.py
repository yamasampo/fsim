import numpy as np

import os
from typing import Union
from fsim.utils import read_control_file, write_settings, write_info_to_file

def main(control_file):
    kw_args = read_control_file(control_file)
    fsim_genetic_drift(**kw_args)

def fsim_genetic_drift(
        pop_size:int, 
        ns:Union[float, int], 
        init_mut_num:int, 
        generation_num:int, 
        output_path_prefix:str,
        output_only_fixation:bool=False,
        total_site_num=0, 
        var_site_num=0, 
        poly_site_num=0, 
        fix_site_num=0
    ):
    """ Returns a list of allele frequency trajectories. `output_only_fixation`
    will control whether the function outputs trajectories for mutations that
    went to fixation. 

    Parameters
    ----------
    pop_size: int
        Population size (number of individuals in a population).
    ns: float or int
        Population size-scaled selection coefficient. 
    init_mut_num: int
        Number of mutant alleles in a population at the initial generation. 
    generation_num: int
        Number of generations. If 0 is given, a simulation will be continued 
        until a mutation goes to fixation or is lost from a population. 
    output_path_prefix: str
        A string text that specifies the output file paths. Two files, 
        "<output_path_prefix>.traj.txt" and "<output_path_prefix>.summ.txt" will 
        be newly generated. 
    output_only_fixation: bool
        Control whether or not write trajectories for mutations that went to
        fixations into an output file. 

    Returns
    -------
    site_count: int
        Number of simulation runs. 
    mutant_freq_trajectories: list
        A list of allele frequency trajectories for all the simulated mutations.    

    """
    # Initialize a list of allele frequency trajectories
    # Trajectory of allele frequency for every simulation replicate will be 
    # added to this list
    mutant_freq_trajectories = []

    # Initialize number of simulated sites (replicates)
    site_count = 0

    # Calculate selection coefficient
    selection_coeff = ns / pop_size

    # Get output file paths
    out_traj_path, out_summ_path = get_output_file_paths(output_path_prefix)

    # Check data type of output_only_fixation
    assert type(output_only_fixation) == bool, \
        'Unknwon value for output_only_fixation argument was found. '\
        'Only "True" or "False" is supported.'

    # Open trajectory output file
    # Simulated frequency trajectory will be consequently added to the file.
    output_fh = open(out_traj_path, 'a')

    # Output setting
    write_settings(
        output_fh, 
        pop_size=pop_size, ns=ns, init_mut_num=init_mut_num, 
        generation_num=generation_num, 
    )
    
    # "[Result]" defines a block for simulation results
    print('[Result]', file=output_fh)
    if total_site_num > 0:
        for _ in range(total_site_num):
            mutant_freq_list = single_rep(
                pop_size, selection_coeff, init_mut_num, generation_num)
            
            mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1

            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue
            
            # Write to file
            res_str = [
                '0' if r == 0 else str(round_num(r, 5)) 
                for r in mutant_freq_list]
            write_info_to_file(output_fh, '\t', *res_str)

    elif var_site_num > 0:
        var_site_count = 0
        
        while var_site_count < var_site_num:
            mutant_freq_list = single_rep(
                pop_size, selection_coeff, init_mut_num, generation_num)
            
            mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1
            if mutant_freq_list[-1] > 0:
                var_site_count += 1
                
            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            # Write to file
            res_str = [
                '0' if r == 0 else str(round_num(r, 5)) 
                for r in mutant_freq_list]
            write_info_to_file(output_fh, '\t', *res_str)
            
    elif poly_site_num > 0:
        poly_site_count = 0
        
        while poly_site_count < poly_site_num:
            mutant_freq_list = single_rep(
                pop_size, selection_coeff, init_mut_num, generation_num)
            
            mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1
            if mutant_freq_list[-1] > 0:
                if mutant_freq_list[-1] < 1:
                    poly_site_count += 1
                
            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            # Write to file
            res_str = [
                '0' if r == 0 else str(round_num(r, 5)) 
                for r in mutant_freq_list]
            write_info_to_file(output_fh, '\t', *res_str)

    elif fix_site_num > 0:
        fix_site_count = 0
        
        while fix_site_count < fix_site_num:
            mutant_freq_list = single_rep(
                pop_size, selection_coeff, init_mut_num, generation_num)
            
            mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1
            if mutant_freq_list[-1] == 1:
                fix_site_count += 1
                
            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            # Write to file
            res_str = [
                '0' if r == 0 else str(round_num(r, 5)) 
                for r in mutant_freq_list]
            write_info_to_file(output_fh, '\t', *res_str)
        
    else:
        raise Exception('Please input integers to total_site_num, var_site_num or poly_site_num.')
    
    output_fh.close()

    assert len(mutant_freq_trajectories) == site_count
    fix_num = len([mut for mut in mutant_freq_trajectories if mut[-1] == 1])
    with open(out_summ_path, 'w') as f:
        print(f'total_rep_num: {site_count}', file=f)
        print(f'fixation_num: {fix_num}', file=f)

    return site_count, mutant_freq_trajectories

def get_output_file_paths(output_path_prefix):
     # Check if a file with the same name as a given output file path does not exist.
    out_traj_path = output_path_prefix+'.traj.txt'
    out_summ_path = output_path_prefix+'.summary.txt'
    if os.path.isfile(out_traj_path):
        raise FileExistsError(out_traj_path)
    if os.path.isfile(out_summ_path):
        raise FileExistsError(out_summ_path)

    return out_traj_path, out_summ_path


def single_rep(
        pop_size, selection_coeff, init_mut_num, 
        generation_num=0, implementation='binomial'
    ):
    """ Returns trajectory of mutant allele frequency. 

    Parameters
    ----------
    pop_size: int
        Population size (number of individuals in a population).
    selection_coeff: float or int
        Selection coefficient. Individuals with mutant and wildtype alleles have
        fitness difference by selection coefficient. 
    init_mut_num: int
        Number of mutant individuals in a population at the initial generation. 
    generation_num: int (default: 0)
        Stops simulation if generation_num greater than 0 is given. 

    Return
    ------
    mutant_freq_list: list
        A trajectory of mutant allele frequency over time. 

    """
    relative_fitness = {
        0: 1 / (2+selection_coeff), # wild type
        1: (1+selection_coeff) / (2+selection_coeff) # mutant
    }
    mutant_freq = init_mut_num / pop_size
    wildtype_freq = (pop_size - init_mut_num) / pop_size

    mutant_freq_list = [mutant_freq]

    gen_num = 0

    while mutant_freq > 0 and mutant_freq < 1:
        assert mutant_freq + wildtype_freq == 1
        exp_mutant_freq = calculate_expected_next_gen_mutant_freq(
            mutant_freq, wildtype_freq, relative_fitness)

        mutant_count = random_sampling(
            pop_size, exp_mutant_freq, implementation)
        
        # Calculate allele frequency
        mutant_freq = mutant_count / pop_size
        wildtype_freq = (pop_size - mutant_count) / pop_size
        
        mutant_freq_list.append(mutant_freq)
        gen_num += 1

        if generation_num > 0:
            # If a given number of generations passed, return tajectory
            if gen_num == generation_num:
                return mutant_freq_list
        
    return mutant_freq_list

def random_sampling(pop_size, exp_mutant_freq, implementation):
    if implementation == 'binomial':
        return np.random.binomial(pop_size, exp_mutant_freq)

    elif implementation == 'exhaustive':
        offsprings = []

        # Construct a population on the next generation
        while len(offsprings) < pop_size:
            if np.random.rand() < exp_mutant_freq:
                offspring = 1 # mutant allele is sampled for next gen
            else:
                offspring = 0 # wildtype allele is sampled for next gen

            offsprings.append(offspring)
                
        # Check if the number of offsprings is the same value as population size
        assert len(offsprings) == pop_size
        return offsprings.count(1)

    raise ValueError(f'Unknown value for `implementation`: {implementation}. '\
        'Only "binomial" or "exhaustive" are supported.')

def calculate_expected_next_gen_mutant_freq(
        mutant_freq, wildtype_freq, relative_fitness):
    """ Returns an expected frequency of mutant alleles in a population at the 
    next generation.
    """

    wildtype_freq_weighted = wildtype_freq * relative_fitness[0]
    mutant_freq_weighted = mutant_freq * relative_fitness[1]

    mutant_freq_weighted_corr = round_num(mutant_freq_weighted / \
        (wildtype_freq_weighted + mutant_freq_weighted), 5)

    wildtype_freq_weighted_corr = round_num(wildtype_freq_weighted / \
        (wildtype_freq_weighted + mutant_freq_weighted), 5)

    # Make sure to sum to 1
    assert mutant_freq_weighted_corr + wildtype_freq_weighted_corr == 1, \
        f"{mutant_freq_weighted_corr} and {wildtype_freq_weighted_corr}"

    return mutant_freq_weighted_corr

def round_num(a, ndigits):
    n = 10 ** ndigits
    return (a * n * 2 + 1) // 2 / n

if __name__ == '__main__':
    import sys
    control_file = sys.argv[1]
    main(control_file)

