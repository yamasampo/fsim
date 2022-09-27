import numpy as np

import os
from typing import Union
from fsim import __version__
from fsim.utils import read_control_file, write_settings, write_info_to_file

def main(control_file):
    kw_args = read_control_file(control_file)
    fsim_genetic_drift(**kw_args)

def fsim_genetic_drift(
        pop_size: int, 
        ns: Union[float, int], 
        init_mut_num: int, 
        generation_num: int, 
        output_path_prefix: str,
        output_only_fixation: bool = False,
        implementation: str = 'binomial',
        total_site_num = 0, 
        var_site_num = 0, 
        poly_site_num = 0, 
        fix_site_num = 0
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
    output_path_prefix: str
        A string text that specifies the output file paths. Two files, 
        "<output_path_prefix>.traj.txt" and "<output_path_prefix>.summ.txt" will 
        be newly generated. 
    generation_num: int
        Number of generations. If 0 is given, a simulation will be continued 
        until a mutation goes to fixation or is lost from a population. 
    output_scale: str
        "number" or "frequency" is supported. "number" will output the number of 
        mutant alleles at each generation, whereas "frequency" will output the
        number divided by population size. 
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
    # mutant_freq_trajectories = []

    # Initialize number of simulated sites (replicates)
    site_count = 0

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
        output_only_fixation=output_only_fixation
    )

    # Calculate relative fitness
    selection_coeff = ns / pop_size
    relative_fitness = {
        0: 1 / (2+selection_coeff), # wild type
        1: (1+selection_coeff) / (2+selection_coeff) # mutant
    }
    exp_freq_d = get_expected_next_gen_freqs_dict(pop_size, relative_fitness)
    
    # "[Result]" defines a block for simulation results
    print('[Result]', file=output_fh)
    if total_site_num > 0:
        fix_site_count = 0

        for _ in range(total_site_num):
            mutant_freq_list = single_rep(
                pop_size, init_mut_num, exp_freq_d, generation_num, 
                implementation)
            
            # mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1
            if mutant_freq_list[-1] == 1:
                fix_site_count += 1

            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue
            
            # Write to file
            res_str = format_result_string(mutant_freq_list, 'frequency')
            write_info_to_file(output_fh, '\t', *res_str)

    elif var_site_num > 0:
        var_site_count = 0
        fix_site_count = 0
        
        while var_site_count < var_site_num:
            mutant_freq_list = single_rep(
                pop_size, init_mut_num, exp_freq_d, generation_num, 
                implementation)
            
            # mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1
            if mutant_freq_list[-1] > 0:
                var_site_count += 1
                if mutant_freq_list[-1] == 1:
                    fix_site_count += 1
                
            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            # Write to file
            res_str = format_result_string(mutant_freq_list, 'frequency')
            write_info_to_file(output_fh, '\t', *res_str)
            
    elif poly_site_num > 0:
        poly_site_count = 0
        fix_site_count = 0
        
        while poly_site_count < poly_site_num:
            mutant_freq_list = single_rep(
                pop_size, init_mut_num, exp_freq_d, generation_num, 
                implementation)
            
            # mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1
            if mutant_freq_list[-1] > 0:
                if mutant_freq_list[-1] == 1:
                    fix_site_count += 1
                elif mutant_freq_list[-1] < 1:
                    poly_site_count += 1
                
                
            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            # Write to file
            res_str = format_result_string(mutant_freq_list, 'frequency')
            write_info_to_file(output_fh, '\t', *res_str)

    elif fix_site_num > 0:
        fix_site_count = 0
        
        while fix_site_count < fix_site_num:
            mutant_freq_list = single_rep(
                pop_size, init_mut_num, exp_freq_d, generation_num, 
                implementation)
            
            # mutant_freq_trajectories.append(mutant_freq_list)
            site_count += 1
            if mutant_freq_list[-1] == 1:
                fix_site_count += 1
                
            if output_only_fixation == True:
                # If a mutation did not go to fixation, do not store.
                if mutant_freq_list[-1] != 1:
                    continue

            # Write to file
            res_str = format_result_string(mutant_freq_list, 'frequency')
            write_info_to_file(output_fh, '\t', *res_str)
        
    else:
        raise Exception('Please input integers to total_site_num, '\
                        'var_site_num or poly_site_num.')
    
    output_fh.close()

    # assert len(mutant_freq_trajectories) == site_count
    # fix_num = len([mut for mut in mutant_freq_trajectories if mut[-1] == 1])
    with open(out_summ_path, 'w') as f:
        print(f'Created by fsim package (version {__version__})', file=f)
        print(f'total_rep_num: {site_count}', file=f)
        print(f'fixation_num: {fix_site_count}', file=f)

    # return site_count, mutant_freq_trajectories
    return site_count

def format_result_string(mutant_freq_list, output_scale):
    res_str = [
        '0' if r == 0 else (
            str(int(r)) if output_scale == 'number' else str(round_num(r, 5)))
        for r in mutant_freq_list
    ]
    return res_str

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
        pop_size: int, init_mut_num: int, exp_mutant_freqs: dict,
        generation_num: int = 0, implementation: str = 'binomial'
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
    output_scale: str
        "number" or "frequency" is supported. "number" will output the number of 
        mutant alleles at each generation, whereas "frequency" will output the
        number divided by population size. 

    Return
    ------
    mutant_freq_list: list
        A trajectory of mutant allele frequency over time. 

    """
    # Store mutant allele frequency in a population at the first 
    # generation.
    mutant_freq = init_mut_num / pop_size
    wildtype_freq = (pop_size - init_mut_num) / pop_size

    mutant_freq_list = [mutant_freq]

    # Initialize the number of generations that are passed. 
    # This will be incremented as simulation proceed. 
    gen_num = 0

    # Continue simulation unless mutant allele is lost or fixed. 
    while mutant_freq > 0 and mutant_freq < 1:
        assert mutant_freq + wildtype_freq == 1
        # Calculate expected mutant allele frequency for the next generation.
        exp_mutant_freq = exp_mutant_freqs[int(mutant_freq*pop_size)]

        # Conduct random binomial sampling of number of mutant alleles for the
        # next generation. 
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

def check_output_scale_arg(output_scale):
    assert output_scale in {'frequency', 'number'}

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

def get_expected_next_gen_freqs_dict(pop_size: int, relative_fitness: float):
    exp_freq_d = {}
    for mutant_count in range(0, pop_size+1):
        wildtype_count = pop_size - mutant_count
        exp_mut_freq = calculate_expected_next_gen_mutant_freq(
            mutant_count/pop_size, wildtype_count/pop_size, relative_fitness)

        exp_freq_d[mutant_count] = exp_mut_freq

    return exp_freq_d

def calculate_expected_next_gen_mutant_freq(
        mutant_freq: float, wildtype_freq: float, relative_fitness: float):
    """ Returns an expected frequency of mutant alleles in a population at the 
    next generation.
    """
    # mutant_freq = mutant_num / pop_size
    # wildtype_freq = (pop_size - mutant_num) / pop_size

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

def round_num(a: Union[int, float], ndigits: int):
    n = 10 ** ndigits
    return (a * n * 2 + 1) // 2 / n

if __name__ == '__main__':
    import sys
    control_file = sys.argv[1]
    main(control_file)

