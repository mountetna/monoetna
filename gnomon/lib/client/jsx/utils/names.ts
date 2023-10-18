import { json_post } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';

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

export function fetchNextCounterValueFromMagma(projectName: string, ruleName: string, tokenPrefix: string): Promise<number> {

    return json_post(magmaPath(`gnomon/${projectName}/increment/${ruleName}/${tokenPrefix}`))
        .then(value => value)
        // remove this once json_post is fixed
        .catch(err => Promise.resolve(err))
        .then(err => Promise.reject(err))
}