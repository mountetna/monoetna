import { createSelector } from 'reselect'
import * as _ from "lodash"

import { RulesState } from '../reducers/rules'
import { NamesState } from '../reducers/names'
import { selectCreateNameGroupIdsByPrimaryRule } from './names'
import { CreateName, Rule } from '../models'
import { defaultDict } from '../utils/object'
import { selectRulesNamesHierarchicalListByPrimaryRuleName } from './rules'
import { State } from '../store'



// returns Rules from all selected CreateNames
// that are seen in all those CreateNames
export const selectCommonRulesFromSelection: (state: State) => Rule[] = createSelector(
    [(state: State) => state],
    (state: State): Rule[] => {
        const selectedCreateNameGroupLocalIds: string[] = [...state.names.createNameGroups.selectionLocalIds]

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
        const visiblePrimaryRuleNames = Object.keys(selectCreateNameGroupIdsByPrimaryRule(state))
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

export const selectReplaceRuleFromSelection: (state: State) => Rule | undefined = createSelector(
    [(state: State) => state],
    (state: State): Rule | undefined => {
        const commonRulesFromSelection = selectCommonRulesFromSelection(state)
        const rulesNamesHierarchicalListByPrimaryRuleName: Record<string, string[]> = selectRulesNamesHierarchicalListByPrimaryRuleName(state)
        const replaceRuleNamesToHierarchy: Record<string, string[]> = {}

        // only the selection's common rules are candidates
        for (const rule of commonRulesFromSelection) {
            replaceRuleNamesToHierarchy[rule.name] = rulesNamesHierarchicalListByPrimaryRuleName[rule.name]
        }

        // find the rule that's deepest in the hierarchy
        if (Object.keys(replaceRuleNamesToHierarchy).length) {
            const replaceRuleName = _.sortBy(replaceRuleNamesToHierarchy, [hierarchy => hierarchy.length]).at(-1).at(-1)
            return state.rules.rules[replaceRuleName]
        }

        return undefined
    }
)