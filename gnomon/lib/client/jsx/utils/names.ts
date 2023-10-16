import { NamesState } from "../reducers/names"



export interface RuleSearchReplaceCriteria {
    createNameTokenValueLocalIds: string[]
    ruleCounterValue: number | undefined
}


export interface SearchReplaceCriteria {
    byRuleName: Record<string, RuleSearchReplaceCriteria>
}


export function createSearchReplaceCriteriaFromGroups(state: NamesState, createNameGroupIds: Iterable<string>): SearchReplaceCriteria[] {
    const criteriaList: SearchReplaceCriteria[] = []

    for (const cngLocalId of createNameGroupIds) {

        const createNameLocalIds: string[] = state.createNames.byCreateNameGroupLocalId[cngLocalId]
        const searchCriteria: SearchReplaceCriteria = { byRuleName: {} }

        createNameLocalIds.forEach(cnLocalId => {
            const cn = state.createNames.byLocalId[cnLocalId]
            const cntvLocalIds = (state.createNameTokenValues.byCreateNameLocalId[cnLocalId] || [])

            searchCriteria.byRuleName[cn.ruleName] = {
                createNameTokenValueLocalIds: cntvLocalIds,
                ruleCounterValue: cn.ruleCounterValue,
            }
        })

        criteriaList.push(searchCriteria)
    }

    return criteriaList
}