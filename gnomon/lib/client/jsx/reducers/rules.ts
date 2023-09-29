/* This isn't used in Rule creation—only reading
from Magma for Name creation.
This should be split into reducers per key,
and eventually used for creation/editing. */

import { Rule, RuleParent, RuleToken, Synonym, Token, TokenValue } from '../models';
import { ACTION_TYPE, ADD_RULES_FROM_MAGMA } from '../actions/rules';
import { listToIdObject, listToIdListObject } from './utils';


// all these child states have duplicate values
// but they're all refs
interface RuleParentsState {
    "byRuleName": Record<string, RuleParent[]>
    "byParentRuleName": Record<string, RuleParent[]>
}


interface RuleTokensState {
    "byRuleName": Record<string, RuleToken[]>
    "byTokenName": Record<string, RuleToken[]>
}


interface TokenValuesState {
    "byName": Record<string, TokenValue>
    "byTokenName": Record<string, TokenValue[]>
}


interface SynonymsState {
    "byValue": Record<string, Synonym>
    "byTokenName": Record<string, Synonym[]>
}


export interface RulesState {
    rules: Record<string, Rule>
    ruleParents: RuleParentsState
    tokens: Record<string, Token>
    ruleTokens: RuleTokensState
    tokenValues: TokenValuesState
    synonyms: SynonymsState
}

const initialState: RulesState = {
    rules: {},
    ruleParents: { byRuleName: {}, byParentRuleName: {} },
    tokens: {},
    ruleTokens: { byTokenName: {}, byRuleName: {} },
    tokenValues: { byName: {}, byTokenName: {} },
    synonyms: { byValue: {}, byTokenName: {} },
}


export function rulesReducer(state: RulesState = initialState, action: ACTION_TYPE): RulesState {
    switch (action.type) {
        case ADD_RULES_FROM_MAGMA:
            return {
                ...state,
                rules: {
                    ...state.rules,
                    ...listToIdObject(action.rules, "name")
                },
                ruleParents: {
                    byRuleName: {
                        ...state.ruleParents.byRuleName,
                        ...listToIdListObject(action.ruleParents, "ruleName"),
                    },
                    byParentRuleName: {
                        ...state.ruleParents.byParentRuleName,
                        ...listToIdListObject(action.ruleParents, "parentRuleName"),
                    },
                },
                tokens: {
                    ...state.tokens,
                    ...listToIdObject(action.tokens, "name")
                },
                ruleTokens: {
                    byTokenName: {
                        ...state.ruleTokens.byTokenName,
                        ...listToIdListObject(action.ruleTokens, "tokenName"),
                    },
                    byRuleName: {
                        ...state.ruleTokens.byRuleName,
                        ...listToIdListObject(action.ruleTokens, "ruleName"),
                    },
                },
                tokenValues: {
                    byName: {
                        ...state.tokenValues.byName,
                        ...listToIdObject(action.tokenValues, "name"),
                    },
                    byTokenName: {
                        ...state.tokenValues.byTokenName,
                        ...listToIdListObject(action.tokenValues, "tokenName"),
                    },
                },
                synonyms: {
                    byValue: {
                        ...state.synonyms.byValue,
                        ...listToIdObject(action.synonyms, "value")
                    },
                    byTokenName: {
                        ...state.synonyms.byTokenName,
                        ...listToIdListObject(action.synonyms, "tokenName")
                    },
                },
            }
        default: {
            return state;
        }
    }
}


