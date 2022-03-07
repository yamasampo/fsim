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
        output_path:str,
        get_only_fixation:bool=False,
        total_site_num=0, 
        var_site_num=0, 
        poly_site_num=0, 
        fix_site_num=0
    ):
    mutant_freq_trajectories = []
    site_count = 0

    selection_coeff = ns / pop_size

    if os.path.isfile(output_path):
        raise FileExistsError(output_path)

    output_fh = open(output_path, 'a')

    write_settings(
        output_fh, 
        pop_size=pop_size, ns=ns, init_mut_num=init_mut_num, 
        generation_num=generation_num, 
    )

    print('[Result]', file=output_fh)
    if total_site_num > 0:
        for _ in range(total_site_num):
            mutant_freq_list = single_rep(
                pop_size, selection_coeff, init_mut_num, generation_num)
            
            site_count += 1

            if get_only_fixation:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue
            mutant_freq_trajectories.append(mutant_freq_list)
            
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
            
            if mutant_freq_list[-1] > 0:
                var_site_count += 1
                
            if get_only_fixation:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            mutant_freq_trajectories.append(mutant_freq_list)

            # Write to file
            res_str = [
                '0' if r == 0 else str(round_num(r, 5)) 
                for r in mutant_freq_list]
            write_info_to_file(output_fh, '\t', *res_str)

            site_count += 1
            
    elif poly_site_num > 0:
        poly_site_count = 0
        
        while poly_site_count < poly_site_num:
            mutant_freq_list = single_rep(
                pop_size, selection_coeff, init_mut_num, generation_num)
            
            if mutant_freq_list[-1] > 0:
                if mutant_freq_list[-1] < 1:
                    poly_site_count += 1
                
            if get_only_fixation:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            mutant_freq_trajectories.append(mutant_freq_list)

            # Write to file
            res_str = [
                '0' if r == 0 else str(round_num(r, 5)) 
                for r in mutant_freq_list]
            write_info_to_file(output_fh, '\t', *res_str)
            site_count += 1
            
    elif fix_site_num > 0:
        fix_site_count = 0
        
        while fix_site_count < fix_site_num:
            mutant_freq_list = single_rep(
                pop_size, selection_coeff, init_mut_num, generation_num)
            
            if mutant_freq_list[-1] == 1:
                fix_site_count += 1
                
            if get_only_fixation:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue
            mutant_freq_trajectories.append(mutant_freq_list)

            # Write to file
            res_str = [
                '0' if r == 0 else str(round_num(r, 5)) 
                for r in mutant_freq_list]
            write_info_to_file(output_fh, '\t', *res_str)
            site_count += 1
        
    else:
        raise Exception('Please input integers to total_site_num, var_site_num or poly_site_num.')
    
    output_fh.close()

    return site_count, mutant_freq_trajectories

def single_rep(pop_size, selection_coeff, init_mut_num, generation_num=0):
    """ Returns trajectory of mutant allele frequency. 

    Parameters
    ----------
    generation_num: int (default: 0)
        Stops simulation if generation_num greater than 0 is given. 

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
        offsprings = [] # 0 is wildtype and 1 is mutant

        assert mutant_freq + wildtype_freq == 1
        _, wildtype_sample_prob = weight_sampling_probability(
            mutant_freq, wildtype_freq, relative_fitness)

        # Construct a population on the next generation
        while len(offsprings) < pop_size:
            if np.random.rand() < wildtype_sample_prob:
                offspring = 0
            else:
                offspring = 1

            offsprings.append(offspring)
            
            # NOTE: wildtype_sample_prob is already weighted by relative fitness.
            # if np.random.rand() < relative_fitness[offspring]:
            #     offsprings.append(offspring)
                
        # Check if the number of offsprings is the same value as population size
        assert len(offsprings) == pop_size
        
        # Calculate allele frequency
        mutant_freq = offsprings.count(1) / pop_size
        wildtype_freq = offsprings.count(0) / pop_size
        
        mutant_freq_list.append(mutant_freq)
        gen_num += 1

        if generation_num > 0:
            # If a given number of generations passed, return tajectory
            if gen_num == generation_num:
                return mutant_freq_list
        
    return mutant_freq_list

def weight_sampling_probability(mutant_freq, wildtype_freq, relative_fitness):
    wildtype_freq_weighted = wildtype_freq * relative_fitness[0]
    mutant_freq_weighted = mutant_freq * relative_fitness[1]

    mutant_freq_weighted_corr = round_num(mutant_freq_weighted / \
        (wildtype_freq_weighted + mutant_freq_weighted), 5)

    wildtype_freq_weighted_corr = round_num(wildtype_freq_weighted / \
        (wildtype_freq_weighted + mutant_freq_weighted), 5)

    assert mutant_freq_weighted_corr + wildtype_freq_weighted_corr == 1, \
        f"{mutant_freq_weighted_corr} and {wildtype_freq_weighted_corr}"

    return mutant_freq_weighted_corr, wildtype_freq_weighted_corr

def round_num(a, ndigits):
    n = 10 ** ndigits
    return (a * n * 2 + 1) // 2 / n

if __name__ == '__main__':
    import sys
    control_file = sys.argv[1]
    main(control_file)

