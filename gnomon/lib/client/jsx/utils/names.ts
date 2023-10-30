import _ from "lodash"

import { json_post } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';

import { NamesState } from "../reducers/names"
import { CreateName, CreateNameTokenValue, RuleToken, TokenValue } from '../models';



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
        .then(value => {
            if (Number.isInteger(Number(value))) {
                return Number.parseInt(value)
            }
            return Promise.reject(`value "${value}" is not an integer`)
        })
}

export function renderTokens(
    createNameTokenValues: CreateNameTokenValue[],
    ruleTokensByLocalId: Record<string, RuleToken>,
    tokenValuesByLocalId: Record<string, TokenValue>,
): string {

    const renderedTokens = [...createNameTokenValues]
        .sort((a, b) => {
            const ordA = ruleTokensByLocalId[a.ruleTokenLocalId].ord
            const ordB = ruleTokensByLocalId[b.ruleTokenLocalId].ord

            return ordA - ordB
        })
        .map(cntv => tokenValuesByLocalId[cntv.tokenValueLocalId].name)
        .reduce((prev, curr) => prev + curr, "")

    return renderedTokens
}

export function renderCounter(createName: CreateName): string {
    return createName.ruleCounterValue != undefined ? String(createName.ruleCounterValue) : ""
}

export function renderCreateName(
    createName: CreateName,
    createNameTokenValues: CreateNameTokenValue[],
    ruleTokensByLocalId: Record<string, RuleToken>,
    tokenValuesByLocalId: Record<string, TokenValue>,
): string {

    const renderedTokens = renderTokens(createNameTokenValues, ruleTokensByLocalId, tokenValuesByLocalId)
    const renderedCounter = createName.ruleCounterValue != undefined ? createName.ruleCounterValue : ""

    return renderedTokens + renderedCounter
}