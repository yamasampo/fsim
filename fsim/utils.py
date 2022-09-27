
import configparser
from warnings import warn
from collections import namedtuple

def read_control_file(control_file):
    # Initialize ConfigParser object
    config = configparser.ConfigParser(
        strict=True,
        comment_prefixes=('/*', ';', '#'),
        inline_comment_prefixes=('/*', ';', '#')
    )

    # Parse control file
    paths = config.read(control_file)

    # Check number of read control files.
    if len(paths) == 0:
        raise FileNotFoundError(
        f'Specified control file, {control_file}, is not found.')
    elif len(paths) > 1:
        raise TypeError(f'Iterable {type(control_file)} is given as a control '\
            'file. Only one control file is supported.')

    # Check sections. Only 'REQUIRED' and 'OPTIONAL' sections will be used.
    assert 'REQUIRED' in config.sections(), \
        f'REQUIRED section is not found in {control_file}.'

    expected_sections = ['REQUIRED', 'OPTIONAL']
    not_expected_sections = [
        s for s in config.sections() if s not in expected_sections]
    if len(not_expected_sections) >= 1:
        msg = f'Unexpected sections, {", ".join(not_expected_sections)}, '\
              'were found. These are not used in '\
              'the analysis. If you wish to include in the analysis, please '\
              'specify in "REQUIRED" or "OPTIONAL" sections.'
        warn(msg)

    converters_d = {
        'pop_size': int,
        'ns': float,
        'init_mut_num': int,
        'generation_num': int,
        'total_site_num': int,
        'var_site_num': int,
        'poly_site_num': int,
        'fix_site_num': int,
        'output_only_fixation': lambda s: True if s == 'True' else (False if s == 'False' else -9)
    }

    flattened  = [
        (opt, converters_d[opt](v)) 
        if opt in converters_d.keys() else (opt, v) 
            for s in expected_sections
                for opt, v in config[s].items()
    ]

    return dict(flattened)

SimulationResult = namedtuple('SimulationResult', ['setting', 'result', 'total_run_num'])

def read_data(file_path, fix_only=True):
    setting_d = {}
    result = []
    run_count = 0
    
    with open(file_path, 'r') as f:
        for l in f:
            if l.startswith('['):
                section = cut_square_parenthesis(l[:-1])
            else:
                if section == 'Setting':
                    key, value = get_argument(l[:-1])
                    setting_d[key] = value
                
                elif section == 'Result':
                    run_count += 1
                    values = get_values(l[:-1])
                    
                    if fix_only:
                        if values[-1] == 1:
                            result.append(values)
                    else:
                        result.append(values)
                    
    return SimulationResult(setting_d, result, run_count)

def cut_square_parenthesis(s):
    return s.split('[')[1].split(']')[0]

def get_argument(s):
    parts = s.split('=')
    assert len(parts) == 2
    return [part.strip() for part in parts]

def get_values(s):
    return [float(v) for v in s.split()]

def write_info_to_file(file_handle, separator, *args, **kw_args):
    """ Write arguments or keyword arguments to a file. Values will be 
    separated by a given separator. 
    """
    output_lines = []
    if len(args) > 0:
        output_lines.append(separator.join(args))

    if len(kw_args) > 0:
        for k, v in kw_args.items():
            output_lines.append(f'{k}{separator}{v}')

    print('\n'.join(output_lines), file=file_handle)

def write_settings(file_handle, **kw_args):
    print('[Setting]', file=file_handle)
    write_info_to_file(file_handle, separator=' = ', **kw_args)

