from typing import Dict, Tuple, List, Callable


def flatten(dictionary, sep='', parent_key=''):
    items = []
    for currKey, currVal in dictionary.items():
        new_key = parent_key + sep + currKey if parent_key else currKey
        if isinstance(currVal, dict):
            items.extend(flatten(currVal, sep=sep, parent_key=new_key).items())
        else:
            items.append((new_key, currVal))
    return dict(items)


def ask_attribute(field: str, magma_attr: str = '') -> str:
    answerOptionsPre = {
        'y': 'yes',
        'Or type the correct attribute name': '  ',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answer_options = format_options(answerOptionsPre)
    question = f'\nIs {magma_attr} correct for {field}?\n' \
               f'ANSWER OPTIONS:\n{answer_options}'
    answer = input(question)
    while not verify_map_format(answer) and answer not in answerOptionsPre.keys():
        verification = f'\nThe answer {answer} is not in format model:attribute\n' \
                       f'Provide a correctly formatted answer'
        answer = input(verification)
    return answer


def ask_characteristics() -> Tuple:
    answerOptionsPre = {
        '1': 'tissue',
        '2': 'cell type',
        '3': 'treatment',
        '4': 'genotype',
        '5': 'disease state'
    }
    answerOptions = format_options(answerOptionsPre)
    answerKey = input(answerOptions)
    while answerKey not in answerOptionsPre.keys():
        verification = f'\nThe answer {answerKey} is not in allowed answers\n' \
                       f'Provide an allowed answer'
        answerKey = input(verification)
    answerValue = input('Enter attribute name for this characteristic')
    while not verify_map_format(answerValue):
        verification = f'\nThe answer {answerValue} is not in format model:attribute\n' \
                       f'Provide a correctly formatted answer'
        answerValue = input(verification)
    return answerOptionsPre[answerKey], answerValue


def add_another() -> str:
    answerOptionsPre = {
        'y': 'yes',
        'n': 'no - same as STOP',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answerOptions = format_options(answerOptionsPre)
    question = f'\nDo you want to add a sample characteristic?\n' \
               f'ANSWER OPTIONS:\n{answerOptions}'
    answer = input(question)
    while answer not in answerOptionsPre.keys():
        verification = f'\nThe answer {answer} is not in allowed answers\n' \
                       f'Provide an allowed answer'
        answer = input(verification)
    return answer


def format_options(options: Dict) -> str:
    formatted = [':\t'.join([x, options[x]]) for x in options]
    answerOptions = '\n'.join(formatted)
    return answerOptions


def characteristics(add_another_func: Callable, d: Dict) -> Dict:
    aw = add_another_func()
    if aw == 'y':
        updateInfo = ask_characteristics()
        dc = d.copy()
        dc.update({updateInfo[0]: updateInfo[1]})
        return characteristics(add_another_func, dc)
    else:
        return d


def verify_map_format(attr_map: str) -> bool:
    if attr_map == '':
        return True
    return len(attr_map.split(':')) == 2


def parse_subject_name(subject_name: str, conventions: List, sep: str = '') -> Dict:
    nameSplit = subject_name.split(sep)
    nameSections = {x: nameSplit[n] for n, x in enumerate(conventions)}
    return nameSections


# TODO


def walk_meetis_dir(metis_dir: str) -> List:
    pass


def get_checksums() -> Dict:
    pass
