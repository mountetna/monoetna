import { createSelector } from 'reselect'

import { Rule, Token } from "../models";
import { RulesState } from '../reducers/rules';


interface State {
    rules: RulesState
}


export const selectRules = (state: State): Record<string, Rule> => state.rules.rules

export const selectRuleByName = (name: string) => (state: State): Rule => state.rules.rules[name]

export const selectTokens = (state: State): Record<string, Token> => state.rules.tokens

export const selectTokenByName = (name: string) => (state: State): Token => state.rules.tokens[name]

// const selectShopItems = state => state.shop.items
// const selectTaxPercent = state => state.shop.taxPercent

// const selectSubtotal = createSelector(selectShopItems, items =>
//   items.reduce((subtotal, item) => subtotal + item.value, 0)
// )