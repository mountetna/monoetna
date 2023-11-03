import { Rule, RuleParent, RuleToken, Synonym, Token, TokenValue } from '../models';
import { makeActionObject } from './utils';



export const ADD_RULES_FROM_MAGMA = "ADD_RULES_FROM_MAGMA"


export function addRulesFromMagma(
    rules: Rule[],
    ruleParents: RuleParent[],
    tokens: Token[],
    ruleTokens: RuleToken[],
    tokenValues: TokenValue[],
    synonyms: Synonym[],
) {

    return makeActionObject(ADD_RULES_FROM_MAGMA, {
        rules,
        ruleParents,
        tokens,
        ruleTokens,
        tokenValues,
        synonyms,
    })
}


export type ACTION_TYPE =
    | ReturnType<typeof addRulesFromMagma>
