/* This isn't used in Rule creation—only reading
from Magma for Name creation.
This should be split into reducers per key,
and eventually used for creation/editing. */

import { Rule, Synonym, Token } from '../models';
import { ACTION_TYPE, ADD_RULES } from '../actions/rules';



interface RulesState {
    rules: Rule[]
    tokens: Token[]
    synonyms: Synonym[]
}

const initialState: RulesState = {
    rules: [],
    tokens: [],
    synonyms: [],
}


export function rulesReducer(state: RulesState = initialState, action: ACTION_TYPE): RulesState {
    switch (action.type) {
        case ADD_RULES:
            return {
                rules: [...state.rules, ...action.rules],
                tokens: [...state.tokens, ...action.tokens],
                synonyms: [...state.synonyms, ...action.synonyms],
            }
        default: {
            return state;
        }
    }
}