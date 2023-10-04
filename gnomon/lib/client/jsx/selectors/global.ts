import { createSelector } from 'reselect'
import { RulesState } from '../reducers/rules'
import { NamesState } from '../reducers/names'
import { selectSelectedCreateNameGroupIds } from './names'
import { CreateName, Rule } from '../models'
import { defaultDict } from '../utils/object'


interface State {
    rules: RulesState
    names: NamesState
}


// returns Rules from all selected CreateNames
// that are seen in all those CreateNames
export const selectIterableRules = createSelector(
    [(state: State) => state],
    (state: State) => {
        const selectedCreateNameGroupLocalIds: string[] = selectSelectedCreateNameGroupIds(state)

        // get CreateNames from selected CreateNameGroups
        const selectedCreateNames: CreateName[] = selectedCreateNameGroupLocalIds
            .map(cngLocalId => state.names.createNames.byCreateNameGroupLocalId[cngLocalId]
                .map(cnLocalId => state.names.createNames.byLocalId[cnLocalId])
            )
            .reduce((prev, curr) => [...prev, ...curr], [])

        const selectedRulesByOccurrence: Record<string, number> = defaultDict<string, number>(_ => 0)

        selectedCreateNames.forEach(cn => {
            ++selectedRulesByOccurrence[cn.ruleName]
        })

        // get Rules from CreateNames that only have one counterValue
        // and seen across all CreateNameGroups
        // (relies on the fact that a Rule can occur only once per CreateNameGroup)
        const fullyRepresentedRules: Rule[] = []

        Object.entries(selectedRulesByOccurrence).forEach(([ruleName, occurrences]) => {
            const rule = state.rules.rules[ruleName]

            if (rule.hasCounter && occurrences == selectedCreateNameGroupLocalIds.length) {
                fullyRepresentedRules.push(rule)
            }
        })

        return fullyRepresentedRules
    }
)