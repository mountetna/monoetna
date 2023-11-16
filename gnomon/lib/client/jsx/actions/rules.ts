import { Rule, RuleParent, RuleToken, Synonym, Token, TokenValue } from '../models';
import { makeActionObject } from './utils';
import { fetchRulesFromMagma as _fetchRulesFromMagma } from '../utils/rules';
import { useDispatch } from '../utils/redux';


export const ADD_RULES_FROM_MAGMA = 'ADD_RULES_FROM_MAGMA';


export async function fetchRulesFromMagma(projectName: string) {
    const rules = await _fetchRulesFromMagma(projectName);

    return (dispatch: ReturnType<typeof useDispatch>) => {
        dispatch(addRulesFromMagma(
            rules.rules,
            rules.ruleParents,
            rules.tokens,
            rules.ruleTokens,
            rules.tokenValues,
            rules.synonyms,
        ));
    };
}


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
    });
}


export type ACTION_TYPE =
    | ReturnType<typeof addRulesFromMagma>
