from typing import Dict, Tuple, List, Callable


def flatten(dictionary, sep='', parent_key=''):
    items = []
    for curr_key, curr_val in dictionary.items():
        new_key = parent_key + sep + curr_key if parent_key else curr_key
        if isinstance(curr_val, dict):
            items.extend(flatten(curr_val, sep=sep, parent_key=new_key).items())
        else:
            items.append((new_key, curr_val))
    return dict(items)


def ask_attribute(field: str, magma_attr: str = '') -> str:
    answer_options_pre = {
        'y': 'yes',
        'Or type the correct attribute name': '  ',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answer_options = format_options(answer_options_pre)
    question = f'\nIs {magma_attr} correct for {field}?\n' \
               f'ANSWER OPTIONS:\n{answer_options}'
    answer = input(question)
    while not verify_map_format(answer) and answer not in answer_options_pre.keys():
        verification = f'\nThe answer {answer} is not in format model:attribute\n' \
                       f'Provide a correctly formatted answer'
        answer = input(verification)
    return answer


def ask_characteristics() -> Tuple:
    answer_options_pre = {
        '1': 'tissue',
        '2': 'cell type',
        '3': 'treatment',
        '4': 'genotype',
        '5': 'disease state'
    }
    answer_options = format_options(answer_options_pre)
    answer_key = input(answer_options)
    while answer_key not in answer_options_pre.keys():
        verification = f'\nThe answer {answer_key} is not in allowed answers\n' \
                       f'Provide an allowed answer'
        answer_key = input(verification)
    answer_value = input('Enter attribute name for this characteristic')
    while not verify_map_format(answer_value):
        verification = f'\nThe answer {answer_value} is not in format model:attribute\n' \
                       f'Provide a correctly formatted answer'
        answer_value = input(verification)
    return answer_options_pre[answer_key], answer_value


def add_another() -> str:
    answer_options_pre = {
        'y': 'yes',
        'n': 'no - same as STOP',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answer_options = format_options(answer_options_pre)
    question = f'\nDo you want to add a sample characteristic?\n' \
               f'ANSWER OPTIONS:\n{answer_options}'
    answer = input(question)
    while answer not in answer_options_pre.keys():
        verification = f'\nThe answer {answer} is not in allowed answers\n' \
                       f'Provide an allowed answer'
        answer = input(verification)
    return answer


def format_options(options: Dict) -> str:
    formatted = [':\t'.join([x, options[x]]) for x in options]
    answer_options = '\n'.join(formatted)
    return answer_options


def characteristics(add_another_func: Callable, d: Dict) -> Dict:
    aw = add_another_func()
    if aw == 'y':
        update_info = ask_characteristics()
        dc = d.copy()
        dc.update({update_info[0]: update_info[1]})
        return characteristics(add_another_func, dc)
    else:
        return d


def verify_map_format(attr_map: str) -> bool:
    if attr_map == '':
        return True
    return len(attr_map.split(':')) == 2


def parse_subject_name(subject_name: str, conventions: List, sep: str = '') -> Dict:
    name_split = subject_name.split(sep)
    name_sections = {x: name_split[n] for n, x in enumerate(conventions)}
    return name_sections


# TODO


def walk_meetis_dir(metis_dir: str) -> List:
    pass


def get_checksums() -> Dict:
    pass
