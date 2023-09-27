import { json_get } from "etna-js/utils/fetch";
import { magmaPath } from "etna-js/api/magma_api";

import { Rule, Synonym, Token, TokenValue } from '../models';
import { makeActionObject } from './utils';



export const ADD_RULES = "ADD_RULES"


export function addRules(rules: Rule[], tokens: Token[], synonyms: Synonym[]) {
    return makeActionObject(ADD_RULES, { rules, tokens, synonyms })
}


interface MagmaRulesResponse {
    rules: object
    tokens: object
    synonyms: Synonym[]
}


interface ParsedRules {
    rules: Rule[]
    tokens: Token[]
    synonyms: Synonym[]
}

const parseMagmaRulesResponse = (res: MagmaRulesResponse): ParsedRules => {
    const tokens = Object.values(res.tokens).map((res_token) => {
        const values = Object.entries(res_token.values).map(([name, label]) => {
            return { name, label } as TokenValue;
        })
        return { name: res_token.name, label: res_token.label, values } as Token;
    });


    const rules: Rule[] = [];
    Object.entries(res.rules).forEach(([name, element_str]) => {
        let hasCounter = false;
        const tokenNames: string[] = []
        const parentRuleNames: string[] = []

        try {
            element_str.split(" ").forEach((el: string) => {
                if (el[0] == ".") {
                    // it's a counter
                    if (el[1] == "n") { hasCounter = true; return }

                    // it's a rule reference
                    parentRuleNames.push(el.slice(1)); return
                }

                // it's a token
                // check if it's a synonym, and replace el with true token name
                // TODO: handle more than 2 values per synonym
                const synonym = res.synonyms.find((synonym) => synonym[0] == el);
                if (synonym) { el = synonym[1] }

                const token = tokens.find(token => token.name == el);
                if (!token) {
                    throw new Error(`Token ${el} does not exist for Rule ${name}. Cannot parse Rule.`)
                }
                tokenNames.push(token.name)
            });

            rules.push({ name, tokenNames, hasCounter, parentRuleNames });
        } catch (e) {
            console.error(e);
        }
    });

    return { rules, tokens, synonyms: res.synonyms };
}


export function fetchRulesFromMagma(project_name) {
    return (dispatch) => {
        json_get(magmaPath(`gnomon/${project_name}`)).then(
            ({ config }: { config: MagmaRulesResponse }) => {
                const parsedRules = parseMagmaRulesResponse(config);
                dispatch(addRules(parsedRules.rules, parsedRules.tokens, parsedRules.synonyms))
            }
        )
    }
}


export type ACTION_TYPE =
    | ReturnType<typeof addRules>
