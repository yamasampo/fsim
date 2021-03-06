
import os
import configparser
from warnings import warn

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

