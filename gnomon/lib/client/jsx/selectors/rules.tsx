import { createSelector } from 'reselect'

import { Rule } from "../models";

export const selectRules = (state): Rule[] => state.rules.rules

// const selectShopItems = state => state.shop.items
// const selectTaxPercent = state => state.shop.taxPercent

// const selectSubtotal = createSelector(selectShopItems, items =>
//   items.reduce((subtotal, item) => subtotal + item.value, 0)
// )