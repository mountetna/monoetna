import _ from 'lodash';

import { json_post, json_get } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';

import { NamesState } from '../reducers/names';
import { CreateName, CreateNameTokenValue, RuleToken, TokenValue } from '../models';
import { createFnConcurrencyWrapper } from './async';
import { Status} from '../utils/models';



export interface RuleSearchReplaceCriteria {
    createNameTokenValueLocalIds: string[]
    ruleCounterValue: number | undefined
}


export interface SearchReplaceCriteria {
    byRuleName: Record<string, RuleSearchReplaceCriteria>
}


export function createSearchReplaceCriteriaFromGroups(state: NamesState, createNameGroupIds: Iterable<string>): SearchReplaceCriteria[] {
    const criteriaList: SearchReplaceCriteria[] = [];

    for (const cngLocalId of createNameGroupIds) {

        const createNameLocalIds: string[] = state.createNames.byCreateNameGroupLocalId[cngLocalId];
        const searchCriteria: SearchReplaceCriteria = { byRuleName: {} };

        createNameLocalIds.forEach(cnLocalId => {
            const cn = state.createNames.byLocalId[cnLocalId];
            const cntvLocalIds = (state.createNameTokenValues.byCreateNameLocalId[cnLocalId] || []);

            searchCriteria.byRuleName[cn.ruleName] = {
                createNameTokenValueLocalIds: cntvLocalIds,
                ruleCounterValue: cn.ruleCounterValue,
            };
        });

        criteriaList.push(searchCriteria);
    }

    return criteriaList;
}

export async function fetchNextCounterValueFromMagma(projectName: string, ruleName: string, tokenPrefix: string): Promise<number> {
    const response = await json_post(magmaPath(`gnomon/${projectName}/increment/${ruleName}/${tokenPrefix}`));

    if (Number.isInteger(Number(response))) {
        return Number.parseInt(response);
    }
    throw new Error(`value "${response}" is not an integer`);
}

export function renderTokens(
    createNameTokenValues: CreateNameTokenValue[],
    ruleTokensByLocalId: Record<string, RuleToken>,
    tokenValuesByLocalId: Record<string, TokenValue>,
): string {

    const renderedTokens = [...createNameTokenValues]
        .sort((a, b) => {
            const ordA = ruleTokensByLocalId[a.ruleTokenLocalId].ord;
            const ordB = ruleTokensByLocalId[b.ruleTokenLocalId].ord;

            return ordA - ordB;
        })
        .map(cntv => tokenValuesByLocalId[cntv.tokenValueLocalId].name)
        .reduce((prev, curr) => prev + curr, '');

    return renderedTokens;
}

export function renderCounter(createName: CreateName): string {
    return createName.ruleCounterValue != undefined ? String(createName.ruleCounterValue) : '';
}

export function renderCreateName(
    createName: CreateName,
    createNameTokenValues: CreateNameTokenValue[],
    ruleTokensByLocalId: Record<string, RuleToken>,
    tokenValuesByLocalId: Record<string, TokenValue>,
): string {

    const renderedTokens = renderTokens(createNameTokenValues, ruleTokensByLocalId, tokenValuesByLocalId);
    const renderedCounter = createName.ruleCounterValue != undefined ? createName.ruleCounterValue : '';

    return renderedTokens + renderedCounter;
}


export interface MagmaBulkGenerateName {
    rule_name: string
    name: string
}


export interface MagmaName {
    project_name: string
    rule: string
    author: string
    identifier: string
    grammar: string
}


export interface MagmaBulkGenerateResponse {
    created: MagmaName[]
    existing?: MagmaName[]
}


export async function postNameBatchToMagma(projectName: string, names: MagmaBulkGenerateName[]): Promise<MagmaBulkGenerateResponse> {
    return await json_post(magmaPath(`gnomon/${projectName}/generate`), { names });
}

export interface MagmaListName {
    identifier: string
    author: string
    name_created_at: string
    record_created_at: string
}

export async function fetchNamesWithRuleAndRegexFromMagma(projectName: string, ruleName: string, regex: string = '.*'): Promise<MagmaListName[]> {
    return await json_get(magmaPath(`gnomon/${projectName}/list/${ruleName}`) + `?regex=${encodeURIComponent(regex)}`);
}

async function _fetchWhetherNameExistsInMagma(projectName: string, ruleName: string, name: string): Promise<boolean> {
    const magmaNames = await fetchNamesWithRuleAndRegexFromMagma(projectName, ruleName, `^${name}$`);
    return magmaNames.length > 0;
}

// TODO: switch to global magma concurrency limit
export const fetchWhetherNameExistsInMagma = createFnConcurrencyWrapper(_fetchWhetherNameExistsInMagma, 4);

export interface MagmaRequestState<T> {
    status: Status
    statusMessage?: string
    response?: T
}