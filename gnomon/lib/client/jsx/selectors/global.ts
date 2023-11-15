import { createSelector } from 'reselect'
import _ from "lodash"

import { selectCreateNameGroupIdsByPrimaryRule } from './names'
import { CreateName, Rule, RuleToken, TokenValue } from '../models'
import { defaultDict } from '../utils/object'
import { selectRuleNamesHierarchicalListByPrimaryRuleName, selectRuleTokenLocalIdsByRuleName, selectRuleTokensByLocalId, selectTokenValuesByLocalId } from './rules'
import { State } from '../store'



export const selectGlobalState = (state: State): State => state

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
        const rulesNamesHierarchicalListByPrimaryRuleName = selectRuleNamesHierarchicalListByPrimaryRuleName(state)

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
        const rulesNamesHierarchicalListByPrimaryRuleName: Record<string, string[]> = selectRuleNamesHierarchicalListByPrimaryRuleName(state)
        const replaceRuleNamesToHierarchy: Record<string, string[]> = {}

        // only the selection's common rules are candidates
        for (const rule of commonRulesFromSelection) {
            replaceRuleNamesToHierarchy[rule.name] = rulesNamesHierarchicalListByPrimaryRuleName[rule.name]
        }

        // find the rule that's deepest in the hierarchy
        const hierarchyLength = Object.keys(replaceRuleNamesToHierarchy).length
        if (Object.keys(replaceRuleNamesToHierarchy).length) {
            const replaceRuleNameHierarchy = _.sortBy(replaceRuleNamesToHierarchy, [hierarchy => hierarchy.length])[hierarchyLength - 1]
            const replaceRuleName = replaceRuleNameHierarchy[replaceRuleNameHierarchy.length - 1]
            return state.rules.rules[replaceRuleName]
        }

        return undefined
    }
)


export interface RulesStateSliceForCompleteCreateNames {
    ruleNamesHierarchicalListByPrimaryRuleName: Record<string, string[]>,
    ruleTokenLocalIdsByRuleName: Record<string, string[]>,
    ruleTokensByLocalId: Record<string, RuleToken>,
    tokenValuesByLocalId: Record<string, TokenValue>,
    counterRequiredByRuleName: Record<string, boolean>,
}


export const selectRulesStateSliceForCompleteCreateNames: (state: State) => RulesStateSliceForCompleteCreateNames = createSelector(
    [(state: State) => state],
    (state: State): RulesStateSliceForCompleteCreateNames => {
        const counterRequiredByRuleName = {}

        for (const [ruleName, rule] of Object.entries(state.rules.rules)) {
            counterRequiredByRuleName[ruleName] = rule.hasCounter
        }

        return {
            ruleNamesHierarchicalListByPrimaryRuleName: selectRuleNamesHierarchicalListByPrimaryRuleName(state),
            ruleTokenLocalIdsByRuleName: selectRuleTokenLocalIdsByRuleName(state),
            ruleTokensByLocalId: selectRuleTokensByLocalId(state),
            tokenValuesByLocalId: selectTokenValuesByLocalId(state),
            counterRequiredByRuleName,
        }
    }
)