import { createSelector } from 'reselect'

import { Rule, RuleParent, RuleToken, Token, TokenValue } from "../models";
import { RulesState } from '../reducers/rules';


interface State {
    rules: RulesState
}


export const selectRulesByName = (state: State): Record<string, Rule> => state.rules.rules

export const selectRuleWithName = (name: string) => (state: State): Rule => state.rules.rules[name]

export const selectRuleParentLocalIdsByRuleName = (state: State): Record<string, string[]> => state.rules.ruleParents.byRuleName

export const selectRuleParentsByLocalId = (state: State): Record<string, RuleParent> => state.rules.ruleParents.byLocalId

export const selectTokens = (state: State): Record<string, Token> => state.rules.tokens

export const selectTokenWithName = (name: string) => (state: State): Token => state.rules.tokens[name]

export const selectTokenValuesByLocalId = (state: State): Record<string, TokenValue> => state.rules.tokenValues.byLocalId

export const selectTokenValueLocalIdsByTokenName = (state: State): Record<string, string[]> => state.rules.tokenValues.byTokenName

export const selectTokenValueLocalIdsWithTokenName = (tokenName: string) => (state: State): string[] => state.rules.tokenValues.byTokenName[tokenName]

export const selectRuleTokensByLocalId = (state: State): Record<string, RuleToken> => state.rules.ruleTokens.byLocalId

export const selectRuleTokenLocalIdsByRuleName = (state: State): Record<string, string[]> => state.rules.ruleTokens.byRuleName

export const selectRuleTokenLocalIdsWithRuleName = (ruleName: string) => (state: State): string[] => state.rules.ruleTokens.byRuleName[ruleName]