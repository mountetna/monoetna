from typing import Dict


def flatten(dictionary, parent_key='', sep=' '):
    items = []
    for currKey, currVal in dictionary.items():
        new_key = parent_key + sep + currKey if parent_key else currKey
        if isinstance(currVal, dict):
            items.extend(flatten(currVal, new_key, sep=sep).items())
        else:
            items.append((new_key, currVal))
    return dict(items)


def askAttribute(field: str, magmaAttr: str) -> str:
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


def askCharacteristics() -> str:
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
    answer = input(answerOptions)
    return answerOptionsPre[answer]


def addAnother() -> str:
    answerOptionsPre = {
        'Y/y': 'yes',
        'n': 'no - same as STOP',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answerOptions = formatOptions(answerOptionsPre)
    answer = input(answerOptions)
    return answer


def formatOptions(options: Dict):
    formatted = [':\t'.join([x, options[x]]) for x in options]
    answerOptions = '\n'.join(formatted)
    return answerOptions
