import { createSelector } from 'reselect'

import { RulesState } from '../reducers/rules'
import { NamesState } from '../reducers/names'
import { selectCreateNameGroupIdsWithPrimaryRule, selectSelectedCreateNameGroupIds } from './names'
import { CreateName, Rule } from '../models'
import { defaultDict } from '../utils/object'
import { selectRulesNamesHierarchicalListByPrimaryRuleName } from './rules'


interface State {
    rules: RulesState
    names: NamesState
}


// returns Rules from all selected CreateNames
// that are seen in all those CreateNames
export const selectCommonRulesFromSelection: (state: State) => Rule[] = createSelector(
    [(state: State) => state],
    (state: State): Rule[] => {
        const selectedCreateNameGroupLocalIds = selectSelectedCreateNameGroupIds(state)

        // get CreateNames from selected CreateNameGroups
        const selectedCreateNames: CreateName[] = selectedCreateNameGroupLocalIds
            .map(cngLocalId => state.names.createNames.byCreateNameGroupLocalId[cngLocalId]
                .map(cnLocalId => state.names.createNames.byLocalId[cnLocalId])
            )
            .reduce((prev, curr) => [...prev, ...curr], [])

        const selectedRulesByOccurrence = defaultDict<string, number>(_ => 0)

        selectedCreateNames.forEach(cn => {
            ++selectedRulesByOccurrence[cn.ruleName]
        })

        // get Rules from CreateNames that are seen across all CreateNameGroups
        // (relies on the fact that a Rule can occur only once per CreateNameGroup)
        const fullyRepresentedRules: Rule[] = []

        Object.entries(selectedRulesByOccurrence).forEach(([ruleName, occurrences]) => {
            const rule = state.rules.rules[ruleName]

            if (occurrences == selectedCreateNameGroupLocalIds.length) {
                fullyRepresentedRules.push(rule)
            }
        })

        return fullyRepresentedRules
    }
)

export const selectVisibleRules: (state: State) => Rule[] = createSelector(
    [(state: State) => state],
    (state: State): Rule[] => {
        const visiblePrimaryRuleNames = Object.keys(selectCreateNameGroupIdsWithPrimaryRule(state))
        const rulesNamesHierarchicalListByPrimaryRuleName = selectRulesNamesHierarchicalListByPrimaryRuleName(state)

        const visibleRules = new Set<Rule>()

        visiblePrimaryRuleNames.forEach(primaryRuleName => {
            rulesNamesHierarchicalListByPrimaryRuleName[primaryRuleName].forEach(ruleName => {

                visibleRules.add(state.rules.rules[ruleName])
            })
        })

        return [...visibleRules]
    }
)