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
    answerOptions = '\n'.join(
        ['Y/y:\tyes',
         'Or type the correct attribute name',
         '0:\tback',
         'STOP:\tcancel pipeline']
    )
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

    formatted = [':\t'.join([x, answerOptionsPre[x]]) for x in answerOptionsPre]
    answerOptions = '\n'.join(formatted)
    answer = input(answerOptions)
    return answerOptionsPre[answer]
    



