import { createSelector } from 'reselect'

import { CompleteCreateName, Rule, RuleParent, RuleToken, Token, TokenValue } from "../models";
import { RulesState } from '../reducers/rules';
import { State } from '../store';



export const selectRulesByName = (state: State): Record<string, Rule> => state.rules.rules

export const selectRuleWithName = (name: string) => (state: State): Rule => state.rules.rules[name]

export const selectRuleParentLocalIdsByRuleName = (state: State): Record<string, string[]> => state.rules.ruleParents.byRuleName

export const selectRuleParentsByLocalId = (state: State): Record<string, RuleParent> => state.rules.ruleParents.byLocalId

export const selectRuleNamesHierarchicalListByPrimaryRuleName: (state: State) => Record<string, string[]> = createSelector(
    [(state: State) => state.rules],
    (rules: RulesState): Record<string, string[]> => {
        const orderedRuleNamesByPrimaryRuleName: Record<string, string[]> = {}

        Object.values(rules.rules).forEach(primaryRule => {

            const orderedRuleNames: string[] = []
            const ruleNamesToScan: string[] = [primaryRule.name]

            while (ruleNamesToScan.length) {
                const ruleName = ruleNamesToScan.pop()
                if (!ruleName) {
                    break
                }
                const rule = rules.rules[ruleName]

                orderedRuleNames.unshift(rule.name)

                // check for rule parents
                if (rules.ruleParents.byRuleName[ruleName]) {
                    const parentRuleNames = rules.ruleParents.byRuleName[ruleName].map(rpLocalId => {

                        return rules.ruleParents.byLocalId[rpLocalId].parentRuleName
                    })
                    ruleNamesToScan.push(...parentRuleNames)
                }
            }

            orderedRuleNamesByPrimaryRuleName[primaryRule.name] = orderedRuleNames
        })

        return orderedRuleNamesByPrimaryRuleName
    }
)

export const selectTokens = (state: State): Record<string, Token> => state.rules.tokens

export const selectTokenWithName = (name: string) => (state: State): Token => state.rules.tokens[name]

export const selectTokenValuesByLocalId = (state: State): Record<string, TokenValue> => state.rules.tokenValues.byLocalId

export const selectTokenValueLocalIdsByTokenName = (state: State): Record<string, string[]> => state.rules.tokenValues.byTokenName

export const selectTokenValueLocalIdsWithTokenName = (tokenName: string) => (state: State): string[] => state.rules.tokenValues.byTokenName[tokenName]

export const selectRuleTokensByLocalId = (state: State): Record<string, RuleToken> => state.rules.ruleTokens.byLocalId

export const selectRuleTokenLocalIdsByRuleName = (state: State): Record<string, string[]> => state.rules.ruleTokens.byRuleName

export const selectRuleTokenLocalIdsWithRuleName = (ruleName: string) => (state: State): string[] => state.rules.ruleTokens.byRuleName[ruleName]