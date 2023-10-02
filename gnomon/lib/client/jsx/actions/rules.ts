import { v4 as uuidv4 } from 'uuid';

import { json_get } from "etna-js/utils/fetch";
import { magmaPath } from "etna-js/api/magma_api";

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


interface MagmaToken {
    name: string
    label: string
    values: Record<string, string>
}


interface MagmaRulesResponse {
    rules: Record<string, string>
    tokens: Record<string, MagmaToken>
    synonyms: string[][]
}


interface ParsedRules {
    rules: Rule[]
    ruleParents: RuleParent[]
    tokens: Token[]
    ruleTokens: RuleToken[]
    tokenValues: TokenValue[]
    synonyms: Synonym[]
}

const parseMagmaRulesResponse = (res: MagmaRulesResponse): ParsedRules => {
    const tokenValues: TokenValue[] = []

    const tokens = Object.values(res.tokens).map((res_token) => {

        Object.entries(res_token.values).forEach(([name, label]) => {
            tokenValues.push({ name, label, tokenName: res_token.name });
        })

        return { name: res_token.name, label: res_token.label } as Token;
    });

    const synonyms: Synonym[] = []

    res.synonyms.forEach((magmaSynonym) => {
        const tokenName = magmaSynonym[0]

        magmaSynonym.slice(1).forEach((synonym) => {
            synonyms.push({ value: synonym, tokenName })
        })
    })

    const rules: Rule[] = []
    const allRuleParents: RuleParent[] = []
    const allRuleTokens: RuleToken[] = []

    Object.entries(res.rules).forEach(([name, element_str]) => {
        let hasCounter = false;
        const ruleParents: RuleParent[] = []
        const ruleTokens: RuleToken[] = []

        try {
            element_str.split(" ").forEach((el) => {
                if (el.startsWith(".")) {
                    // it's a counter
                    if (el[1] == "n") { hasCounter = true; return }

                    // it's a rule reference
                    // TODO: create RuleParents after all rules to validate parent exists
                    ruleParents.push({
                        localId: uuidv4(),
                        ruleName: name,
                        parentRuleName: el.slice(1),
                        ord: ruleParents.length
                    })
                    return
                }

                // check if synonym, then replace with true token name
                for (const synonym of synonyms) {
                    if (synonym.value == el) {
                        el = synonym.tokenName
                        break
                    }
                }

                // it's a token
                const token = tokens.find(token => token.name == el);
                if (!token) {
                    throw new Error(`Token ${el} does not exist for Rule ${name}. Cannot parse Rule.`)
                }
                ruleTokens.push({
                    localId: uuidv4(),
                    ruleName: name,
                    tokenName: token.name,
                    ord: ruleTokens.length
                })
            });

            rules.push({ name, hasCounter });
            allRuleParents.push(...ruleParents)
            allRuleTokens.push(...ruleTokens)
        } catch (e) {
            console.error(e);
        }
    });

    return {
        rules,
        ruleParents: allRuleParents,
        tokens,
        ruleTokens: allRuleTokens,
        tokenValues,
        synonyms,
    }
}


export function fetchRulesFromMagma(project_name) {
    return (dispatch) => {
        json_get(magmaPath(`gnomon/${project_name}`)).then(
            ({ config }: { config: MagmaRulesResponse }) => {
                const parsedRules = parseMagmaRulesResponse(config);
                dispatch(addRulesFromMagma(
                    parsedRules.rules,
                    parsedRules.ruleParents,
                    parsedRules.tokens,
                    parsedRules.ruleTokens,
                    parsedRules.tokenValues,
                    parsedRules.synonyms,
                ))
            }
        )
    }
}


export type ACTION_TYPE =
    | ReturnType<typeof addRulesFromMagma>
