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


def askAttribute(field: str, magmaAttr: str='') -> str:
    answerOptionsPre = {
        'y': 'yes',
        'Or type the correct attribute name': '  ',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answerOptions = formatOptions(answerOptionsPre)
    question = f'\nIs {magmaAttr} correct for {field}?\n' \
               f'ANSWER OPTIONS:\n{answerOptions}'
    answer = input(question)
    while not verifyMapFormat(answer) and answer not in answerOptionsPre.keys():
        verification = f'\nThe answer {answer} is not in format model:attribute\n' \
                       f'Provide a correctly formatted answer'
        answer = input(verification)
    return answer

def askCharacteristics() -> Tuple:
    answerOptionsPre = {
        '1': 'tissue',
        '2': 'cell type',
        '3': 'treatment',
        '4': 'genotype',
        '5': 'disease state'
    }
    answerOptions = formatOptions(answerOptionsPre)
    answerKey = input(answerOptions)
    while answerKey not in answerOptionsPre.keys():
        verification = f'\nThe answer {answerKey} is not in allowed answers\n' \
                       f'Provide an allowed answer'
        answerKey = input(verification)
    answerValue = input('Enter attribute name for this characteristic')
    while not verifyMapFormat(answerValue):
        verification = f'\nThe answer {answerValue} is not in format model:attribute\n' \
                       f'Provide a correctly formatted answer'
        answerValue = input(verification)
    return answerOptionsPre[answerKey], answerValue


def addAnother() -> str:
    answerOptionsPre = {
        'y': 'yes',
        'n': 'no - same as STOP',
        '0': 'back',
        'STOP': 'cancel pipeline'
    }
    answerOptions = formatOptions(answerOptionsPre)
    question = f'\nDo you want to add a sample characteristic?\n' \
               f'ANSWER OPTIONS:\n{answerOptions}'
    answer = input(question)
    while answer not in answerOptionsPre.keys():
        verification = f'\nThe answer {answer} is not in allowed answers\n' \
                       f'Provide an allowed answer'
        answer = input(verification)
    return answer


def formatOptions(options: Dict) -> str:
    formatted = [':\t'.join([x, options[x]]) for x in options]
    answerOptions = '\n'.join(formatted)
    return answerOptions

def characteristics(addAnother: Callable, d: Dict) -> Dict:
    aw = addAnother()
    if aw == 'y':
        updateInfo = askCharacteristics()
        dc = d.copy()
        dc.update({updateInfo[0]: updateInfo[1]})
        return characteristics(addAnother, dc)
    else:
        return d


def verifyMapFormat(attrMap: str) -> bool:
    if attrMap == '':
        return True
    return len(attrMap.split(':')) == 2


def parseSubjectName(subjectName: str, conventions: List, sep: str='') -> Dict:
    nameSplit = subjectName.split(sep)
    nameSections = {x: nameSplit[n] for n,x in enumerate(conventions)}
    return nameSections


### TODO


def walkMetisDir(metisDir: str) -> List:
    pass


def getChecksums() -> Dict:
    pass

