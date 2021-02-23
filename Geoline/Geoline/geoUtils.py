from typing import Dict, Tuple


def flatten(dictionary, sep='', parent_key=''):
    items = []
    for currKey, currVal in dictionary.items():
        new_key = parent_key + sep + currKey if parent_key else currKey
        if isinstance(currVal, dict):
            items.extend(flatten(currVal, sep=sep, parent_key=new_key).items())
        else:
            items.append((new_key, currVal))
    return dict(items)


def askAttribute(field: str, magmaAttr: str='') -> str:
    answerOptionsPre = {
        'Y/y': 'yes',
        'Or type the correct attribute name': '_',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answerOptions = formatOptions(answerOptionsPre)
    question = f'\nIs {magmaAttr} correct for {field}?\n' \
               f'ANSWER OPTIONS:\n{answerOptions}'
    answer = input(question)
    return answer


def askCharacteristics() -> Tuple:
    answerOptionsPre = {
        '1': 'tissue',
        '2': 'cell type',
        '3': 'treatment',
        '4': 'genotype',
        '5': 'disease state',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answerOptions = formatOptions(answerOptionsPre)
    answerKey = input(answerOptions)
    if answerKey not in answerOptionsPre.keys():
        raise KeyError(f'askCharacteristics(): selected answer {answerKey} is not an allowed option'
                       f'Choose one from {list(answerOptionsPre.keys())}')
    answerValue = input('Enter attribute name for this characteristic')
    return answerOptionsPre[answerKey], answerValue


def addAnother() -> str:
    answerOptionsPre = {
        'y': 'yes',
        'n': 'no - same as STOP',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answerOptions = formatOptions(answerOptionsPre)
    answer = input(answerOptions)
    if answer not in answerOptionsPre.keys():
        raise KeyError(f'addAnother(): selected answer {answer} is not an allowed option'
                       f'Choose one from {list(answerOptionsPre.keys())}')
    return answer


def formatOptions(options: Dict) -> str:
    formatted = [':\t'.join([x, options[x]]) for x in options]
    answerOptions = '\n'.join(formatted)
    return answerOptions
