import { json_get } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';

import { Rule, RuleParent, RuleToken, Synonym, Token, TokenValue } from '../models';
import { createLocalId } from './models';



export interface MagmaToken {
    name: string
    label: string
    values: Record<string, string>
}


export interface MagmaRulesResponse {
    config: MagmaRules
}


export interface MagmaRules {
    rules: Record<string, string>
    tokens: Record<string, MagmaToken>
    synonyms: string[][]
}


export interface ParsedRules {
    rules: Rule[]
    ruleParents: RuleParent[]
    tokens: Token[]
    ruleTokens: RuleToken[]
    tokenValues: TokenValue[]
    synonyms: Synonym[]
}

const parseMagmaRulesResponse = (res: MagmaRules): ParsedRules => {
    const tokenValues: TokenValue[] = [];

    const tokens = Object.values(res.tokens).map((res_token) => {

        Object.entries(res_token.values).forEach(([name, label]) => {
            tokenValues.push({
                localId: createLocalId(),
                name,
                label,
                tokenName: res_token.name
            });
        });

        return { name: res_token.name, label: res_token.label } as Token;
    });

    const synonyms: Synonym[] = [];

    res.synonyms.forEach((magmaSynonym) => {
        const tokenName = magmaSynonym[0];

        magmaSynonym.slice(1).forEach((synonym) => {
            synonyms.push({ value: synonym, tokenName });
        });
    });

    const rules: Rule[] = [];
    const allRuleParents: RuleParent[] = [];
    const allRuleTokens: RuleToken[] = [];

    Object.entries(res.rules).forEach(([name, element_str]) => {
        let hasCounter = false;
        const ruleParents: RuleParent[] = [];
        const ruleTokens: RuleToken[] = [];

        try {
            element_str.split(' ').forEach((el) => {
                if (el.startsWith('.')) {
                    // it's a counter
                    if (el[1] == 'n') { hasCounter = true; return; }

                    // it's a rule reference
                    // TODO: create RuleParents after all rules to validate parent exists
                    ruleParents.push({
                        localId: createLocalId(),
                        ruleName: name,
                        parentRuleName: el.slice(1),
                        ord: ruleParents.length
                    });
                    return;
                }

                // check if synonym, then replace with true token name
                for (const synonym of synonyms) {
                    if (synonym.value == el) {
                        el = synonym.tokenName;
                        break;
                    }
                }

                // it's a token
                const token = tokens.find(token => token.name == el);
                if (!token) {
                    throw new Error(`Token ${el} does not exist for Rule ${name}. Cannot parse Rule.`);
                }
                ruleTokens.push({
                    localId: createLocalId(),
                    ruleName: name,
                    tokenName: token.name,
                    ord: ruleTokens.length
                });
            });

            rules.push({ name, hasCounter });
            allRuleParents.push(...ruleParents);
            allRuleTokens.push(...ruleTokens);
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
    };
};


export async function fetchRulesFromMagma(projectName: string): Promise<ParsedRules> {
    const res: MagmaRulesResponse = await json_get(magmaPath(`gnomon/${projectName}`));

    return parseMagmaRulesResponse(res.config);
}