/* This isn't used in Rule creation—only reading
from Magma for Name creation.
This should be split into reducers per key,
and eventually used for creation/editing. */

import { Rule, RuleParent, RuleToken, Synonym, Token, TokenValue } from '../models';
import { ACTION_TYPE, ADD_RULES_FROM_MAGMA } from '../actions/rules';
import { listToIdObject, listToIdGroupObject } from '../utils/object';


interface RuleParentsState {
    "byLocalId": Record<string, RuleParent>
    "byRuleName": Record<string, string[]>
    "byParentRuleName": Record<string, string[]>
}


interface RuleTokensState {
    "byLocalId": Record<string, RuleToken>
    "byRuleName": Record<string, string[]>
    "byTokenName": Record<string, string[]>
}


interface TokenValuesState {
    "byName": Record<string, TokenValue>
    "byTokenName": Record<string, string[]>  // { TokenValue.tokenName: TokenValue.name }
}


interface SynonymsState {
    "byValue": Record<string, Synonym>
    "byTokenName": Record<string, string[]>  // { Synonym.tokenName: Synonym.value }
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
    ruleParents: { byLocalId: {}, byRuleName: {}, byParentRuleName: {} },
    tokens: {},
    ruleTokens: { byLocalId: {}, byTokenName: {}, byRuleName: {} },
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
                    byLocalId: {
                        ...state.ruleParents.byLocalId,
                        ...listToIdObject(action.ruleParents, "localId"),
                    },
                    byRuleName: {
                        ...state.ruleParents.byRuleName,
                        ...listToIdGroupObject(action.ruleParents, "ruleName", "localId"),
                    },
                    byParentRuleName: {
                        ...state.ruleParents.byParentRuleName,
                        ...listToIdGroupObject(action.ruleParents, "parentRuleName", "localId"),
                    },
                },
                tokens: {
                    ...state.tokens,
                    ...listToIdObject(action.tokens, "name")
                },
                ruleTokens: {
                    byLocalId: {
                        ...state.ruleTokens.byLocalId,
                        ...listToIdObject(action.ruleTokens, "localId")
                    },
                    byTokenName: {
                        ...state.ruleTokens.byTokenName,
                        ...listToIdGroupObject(action.ruleTokens, "tokenName", "localId"),
                    },
                    byRuleName: {
                        ...state.ruleTokens.byRuleName,
                        ...listToIdGroupObject(action.ruleTokens, "ruleName", "localId"),
                    },
                },
                tokenValues: {
                    byName: {
                        ...state.tokenValues.byName,
                        ...listToIdObject(action.tokenValues, "name"),
                    },
                    byTokenName: {
                        ...state.tokenValues.byTokenName,
                        ...listToIdGroupObject(action.tokenValues, "tokenName", "name"),
                    },
                },
                synonyms: {
                    byValue: {
                        ...state.synonyms.byValue,
                        ...listToIdObject(action.synonyms, "value")
                    },
                    byTokenName: {
                        ...state.synonyms.byTokenName,
                        ...listToIdGroupObject(action.synonyms, "tokenName", "value")
                    },
                },
            }
        default: {
            return state;
        }
    }
}


