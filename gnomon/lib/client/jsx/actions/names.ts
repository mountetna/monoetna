import _ from 'lodash';

import { CreateName, CreateNameGroup, CreateNameTokenValue } from '../models';
import { makeActionObject } from './utils';
import { Status, createLocalId } from '../utils/models';
import { defaultDict } from '../utils/object';
import { NamesCreationRequestState, NamesListRequestState, NamesState } from '../reducers/names/reducer';
import { createSearchReplaceCriteriaFromGroups, postNameBatchToMagma } from '../utils/names';
import { State } from '../store';
import { RulesStateSliceForCompleteCreateNames, selectRulesStateSliceForCompleteCreateNames } from '../selectors/global';
import { CompleteCreateNameRequestPayload } from '../selectors/names';
import { useDispatch } from '../utils/redux';



export const ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES = 'ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES';
export const ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES = 'ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES';
export const DELETE_CREATE_NAME_TOKEN_VALUE = 'DELETE_CREATE_NAME_TOKEN_VALUE';
export const SET_COUNTER_VALUE_FOR_CREATE_NAMES = 'SET_COUNTER_VALUE_FOR_CREATE_NAMES';
export const ADD_CREATE_NAME_GROUPS_TO_SELECTION = 'ADD_CREATE_NAME_GROUPS_TO_SELECTION';
export const REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION = 'REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION';
export const CLEAR_CREATE_NAME_GROUPS_SELECTION = 'CLEAR_CREATE_NAME_GROUPS_SELECTION';
export const CLEAR_CREATE_NAME_GROUPS_FILTER = 'CLEAR_CREATE_NAME_GROUPS_FILTER';
export const DELETE_GROUPS_WITH_NAMES = 'DELETE_GROUPS_WITH_NAMES';
export const DELETE_SELECTED_GROUPS_WITH_NAMES = 'DELETE_SELECTED_GROUPS_WITH_NAMES';
export const ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA = 'ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA';
export const REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA = 'REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA';
export const SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA = 'SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA';
export const SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA = 'SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA';
export const ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA = 'ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA';
export const REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA = 'REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA';
export const SET_MAGMA_NAMES_CREATION_REQUEST = 'SET_MAGMA_NAMES_CREATION_REQUEST';
export const SET_MAGMA_NAMES_LIST_REQUEST = 'SET_MAGMA_NAMES_LIST_REQUEST';
export const SET_SEARCH_VISIBILITY = 'SET_SEARCH_VISIBILITY';
export const SET_REPLACE_VISIBILITY = 'SET_REPLACE_VISIBILITY';
export const SET_MAGMA_CHECK_DUPLICATE_NAME_REQUEST = 'SET_MAGMA_CHECK_DUPLICATE_NAME_REQUEST';
export const SET_MAGMA_INCREMENT_COUNTER_REQUEST = 'SET_MAGMA_INCREMENT_COUNTER_REQUEST';


export function addNamesWithGroupsWithTokenValues(
    createNames: CreateName[],
    createNameGroups: CreateNameGroup[],
    createNameTokenValues: CreateNameTokenValue[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
) {
    return makeActionObject(ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES, {
        createNames,
        createNameGroups,
        createNameTokenValues,
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames,
    });
}

export function createNamesWithGroupForRule(
    primaryRuleName: string,
    state: State,
    updateCompletionStatus: boolean,
) {

    const rulesState = state.rules;
    const createNames: CreateName[] = [];
    const ruleNamesToScan: string[] = [primaryRuleName];

    const createNameGroup: CreateNameGroup = {
        localId: createLocalId(),
        primaryCreateNameLocalId: 'TODO',
    };

    const createNameTokenValues: CreateNameTokenValue[] = [];

    // make all CreateNames
    while (ruleNamesToScan.length) {
        const ruleName = ruleNamesToScan.pop();
        // typescript needs this check even though
        // we're already checking ruleNamesToScan.length
        if (!ruleName) {
            break;
        }
        const createName: CreateName = {
            localId: createLocalId(),
            ruleName,
            createNameGroupLocalId: createNameGroup.localId
        };

        createNames.push(createName);
        if (ruleName == primaryRuleName) {
            createNameGroup.primaryCreateNameLocalId = createName.localId;
        }

        // make CreateNameTokenValue if only one value for a given token
        const ruleTokens = rulesState.ruleTokens.byRuleName[ruleName].map(rtLocalId => rulesState.ruleTokens.byLocalId[rtLocalId]);
        ruleTokens.forEach(ruleToken => {
            const tokenValueLocalIds = rulesState.tokenValues.byTokenName[ruleToken.tokenName];

            if (tokenValueLocalIds.length == 1) {
                createNameTokenValues.push({
                    localId: createLocalId(),
                    tokenValueLocalId: tokenValueLocalIds[0],
                    createNameLocalId: createName.localId,
                    ruleTokenLocalId: ruleToken.localId,
                });
            }
        });

        if (rulesState.ruleParents.byRuleName[ruleName]) {
            const parentRuleNames = rulesState.ruleParents.byRuleName[ruleName].map(rpLocalId => {
                return rulesState.ruleParents.byLocalId[rpLocalId].parentRuleName;
            });
            // TODO: properly handle multiple parents
            // (will be ok with one for now)
            ruleNamesToScan.push(...parentRuleNames);
        }
    }

    return addNamesWithGroupsWithTokenValues(
        createNames,
        [createNameGroup],
        createNameTokenValues,
        updateCompletionStatus,
        selectRulesStateSliceForCompleteCreateNames(state)
    );
}


interface DuplicateCreateNameGroupReturn {
    createNames: CreateName[]
    createNameGroup: CreateNameGroup
    createNameTokenValues: CreateNameTokenValue[]
}


function _duplicateCreateNameGroup(
    createNameGroup: CreateNameGroup,
    namesState: NamesState,
    ruleCounterValuesByRuleName: Record<string, number> = {},
): DuplicateCreateNameGroupReturn {

    const newCng: CreateNameGroup = {
        localId: createLocalId(),
        primaryCreateNameLocalId: 'TODO',
    };
    const newCntvs: CreateNameTokenValue[] = [];

    const newCns: CreateName[] = namesState.createNames.byCreateNameGroupLocalId[createNameGroup.localId].map(oldCnLocalId => {
        const oldCn = namesState.createNames.byLocalId[oldCnLocalId];

        const newCn: CreateName = {
            ...oldCn,
            localId: createLocalId(),
            createNameGroupLocalId: newCng.localId,
        };

        if (newCn.ruleName in ruleCounterValuesByRuleName) {
            newCn.ruleCounterValue = ruleCounterValuesByRuleName[newCn.ruleName];
        }

        if (createNameGroup.primaryCreateNameLocalId == oldCn.localId) {
            newCng.primaryCreateNameLocalId = newCn.localId;
        }

        const _newCnTvs: CreateNameTokenValue[] = (namesState.createNameTokenValues.byCreateNameLocalId[oldCn.localId] || []).map(oldCntvLocalId => {
            const oldCntv = namesState.createNameTokenValues.byLocalId[oldCntvLocalId];

            return {
                ...oldCntv,
                localId: createLocalId(),
                createNameLocalId: newCn.localId,
            };
        });
        newCntvs.push(..._newCnTvs);

        return newCn;
    });

    return { createNames: newCns, createNameGroup: newCng, createNameTokenValues: newCntvs };
}

// maybe just pass entire state?
export function duplicateCreateNameGroups(
    createNameGroups: CreateNameGroup[],
    state: State,
    updateCompletionStatus: boolean = true,
    quantity: number = 1,
    overrideRuleCounterValuesByRuleName: Record<string, number> = {},
) {

    const newCreateNames: CreateName[] = [];
    const newCreateNameGroups: CreateNameGroup[] = [];
    const newCreateNameTokenValues: CreateNameTokenValue[] = [];

    createNameGroups.forEach(cng => {

        for (let index = 1; index <= quantity; index++) {

            const {
                createNames: _newCreateNames,
                createNameGroup: _newCreateNameGroup,
                createNameTokenValues: _newCreateNameTokenValues
            } = _duplicateCreateNameGroup(
                cng,
                state.names,
                overrideRuleCounterValuesByRuleName,
            );

            newCreateNames.push(..._newCreateNames);
            newCreateNameGroups.push(_newCreateNameGroup);
            newCreateNameTokenValues.push(..._newCreateNameTokenValues);
        }
    });

    return addNamesWithGroupsWithTokenValues(
        newCreateNames,
        newCreateNameGroups,
        newCreateNameTokenValues,
        updateCompletionStatus,
        selectRulesStateSliceForCompleteCreateNames(state),
    );
}

export function iterateOnCreateNameGroupsByRule(
    createNameGroups: CreateNameGroup[],
    ruleName: string,
    startIdx: number,
    endIdx: number,
    state: State,
    updateCompletionStatus: boolean = true,
) {

    const createNameGroupsByPrimaryRuleName = defaultDict<string, CreateNameGroup[]>(_ => []);

    createNameGroups.forEach(cng => {
        const primaryRuleName = state.names.createNames.byLocalId[cng.primaryCreateNameLocalId].ruleName;

        createNameGroupsByPrimaryRuleName[primaryRuleName].push(cng);
    });

    const newCreateNames: CreateName[] = [];
    const newCreateNameGroups: CreateNameGroup[] = [];
    const newCreateNameTokenValues: CreateNameTokenValue[] = [];

    Object.entries(createNameGroupsByPrimaryRuleName).forEach(([_, createNameGroups]) => {

        for (let index = startIdx; index <= endIdx; index++) {

            createNameGroups.forEach(cng => {

                const overrideRuleCounterValuesByRuleName = { [ruleName]: index };

                const {
                    createNames: _newCreateNames,
                    createNameGroup: _newCreateNameGroup,
                    createNameTokenValues: _newCreateNameTokenValues
                } = _duplicateCreateNameGroup(
                    cng,
                    state.names,
                    overrideRuleCounterValuesByRuleName,
                );

                newCreateNames.push(..._newCreateNames);
                newCreateNameGroups.push(_newCreateNameGroup);
                newCreateNameTokenValues.push(..._newCreateNameTokenValues);
            });
        }
    });

    return addNamesWithGroupsWithTokenValues(
        newCreateNames,
        newCreateNameGroups,
        newCreateNameTokenValues,
        updateCompletionStatus,
        selectRulesStateSliceForCompleteCreateNames(state),
    );
}

export function addOrReplaceCreateNameTokenValues(
    newCreateNameTokenValues: CreateNameTokenValue[],
    oldCreateNameTokenValueLocalIds: string[],
    state: State,
    updateCompletionStatus: boolean = true,
) {
    return makeActionObject(ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES, {
        newCreateNameTokenValues,
        oldCreateNameTokenValueLocalIds,
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames: selectRulesStateSliceForCompleteCreateNames(state),
    });
}

export function deleteCreateNameTokenValue(
    createNameTokenValue: CreateNameTokenValue,
    state: State,
    updateCompletionStatus: boolean = true,
) {
    return makeActionObject(DELETE_CREATE_NAME_TOKEN_VALUE, {
        createNameTokenValue,
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames: selectRulesStateSliceForCompleteCreateNames(state),
    });
}

export function setCreateNameRuleCounterValues(
    ruleCounterValuesByCreateNameLocalId: Record<string, number | undefined>,
    state: State,
    updateCompletionStatus: boolean = true,
) {
    return makeActionObject(SET_COUNTER_VALUE_FOR_CREATE_NAMES, {
        ruleCounterValuesByCreateNameLocalId,
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames: selectRulesStateSliceForCompleteCreateNames(state),
    });
}

export function addCreateNameGroupsToSelection(createNameGroupIds: string[]) {
    return makeActionObject(ADD_CREATE_NAME_GROUPS_TO_SELECTION, { createNameGroupIds });
}

export function removeCreateNameGroupsFromSelection(createNameGroupIds: string[]) {
    return makeActionObject(REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION, { createNameGroupIds });
}

export function clearCreateNameGroupsSelection() {
    return makeActionObject(CLEAR_CREATE_NAME_GROUPS_SELECTION, {});
}

export function setCreateNameGroupsSelectionFromSearchCriteria() {
    return makeActionObject(SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA, {});
}

export function setCreateNameGroupsFilterFromSearchCriteria() {
    return makeActionObject(SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA, {});
}

export function clearCreateNameGroupsFilter() {
    return makeActionObject(CLEAR_CREATE_NAME_GROUPS_FILTER, {});
}

export function deleteGroupsWithNames(
    createNameGroupIds: string[],
    state: State,
    updateCompletionStatus: boolean = true,
) {
    return makeActionObject(DELETE_GROUPS_WITH_NAMES, {
        createNameGroupIds,
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames: selectRulesStateSliceForCompleteCreateNames(state),
    });
}

export function deleteSelectedGroupsWithNames(
    state: State,
    updateCompletionStatus: boolean = true,
) {
    return makeActionObject(DELETE_SELECTED_GROUPS_WITH_NAMES, {
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames: selectRulesStateSliceForCompleteCreateNames(state),
    });
}

export function addCreateNameGroupsToSearchCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA, { createNameGroupLocalIds });
}

export function removeCreateNameGroupsFromSearchCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA, { createNameGroupLocalIds });
}

export function addCreateNameGroupsToReplaceCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA, { createNameGroupLocalIds });
}

export function removeCreateNameGroupsFromReplaceCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA, { createNameGroupLocalIds });
}

export function addOrReplaceCreateNameTokenValuesAndRuleCounterValuesFromReplaceCriteria(
    state: State,
    updateCompletionStatus: boolean = true,
) {
    const namesState = state.names;

    // TODO: support multiple CreateNameGroups for criteria?
    const replaceCriteria = createSearchReplaceCriteriaFromGroups(
        namesState, namesState.createNameGroups.replaceLocalIds
    )[0].byRuleName;

    const newCntvs: CreateNameTokenValue[] = [];
    const oldCntvLocalIds: string[] = [];
    const newCounterValuesByCnLocalId: Record<string, number | undefined> = {};

    for (const selectionGroupId of namesState.createNameGroups.selectionLocalIds) {
        const selectionCnLocalIds = namesState.createNames.byCreateNameGroupLocalId[selectionGroupId];

        for (const selectionCnLocalId of selectionCnLocalIds) {
            const selectionCn = namesState.createNames.byLocalId[selectionCnLocalId];

            const selectionCntvs = (namesState.createNameTokenValues.byCreateNameLocalId[selectionCnLocalId] || [])
                .map(cntvLocalId => namesState.createNameTokenValues.byLocalId[cntvLocalId]);
            const replaceCntvs = (replaceCriteria[selectionCn.ruleName]?.createNameTokenValueLocalIds || [])
                .map(cntvLocalId => namesState.createNameTokenValues.byLocalId[cntvLocalId]);

            for (const replaceCntv of replaceCntvs) {

                // create or replace TokenValues
                const selectionCntv: CreateNameTokenValue | undefined = selectionCntvs.find(cntv => {
                    return cntv.ruleTokenLocalId == replaceCntv.ruleTokenLocalId;
                });

                if (selectionCntv?.tokenValueLocalId == replaceCntv.tokenValueLocalId) {
                    continue;
                }

                newCntvs.push({
                    localId: createLocalId(),
                    tokenValueLocalId: replaceCntv.tokenValueLocalId,
                    createNameLocalId: selectionCn.localId,
                    ruleTokenLocalId: replaceCntv.ruleTokenLocalId,
                });

                if (selectionCntv) {
                    oldCntvLocalIds.push(selectionCntv.localId);
                }
            }

            // update counterValue
            const newCounterValue = replaceCriteria[selectionCn.ruleName]?.ruleCounterValue;

            if (newCounterValue != undefined && selectionCn.ruleCounterValue != newCounterValue) {
                newCounterValuesByCnLocalId[selectionCn.localId] = newCounterValue;
            }
        }
    }

    return [
        addOrReplaceCreateNameTokenValues(newCntvs, oldCntvLocalIds, state, updateCompletionStatus),
        setCreateNameRuleCounterValues(newCounterValuesByCnLocalId, state, updateCompletionStatus),
    ];
}


function setMagmaNamesCreationRequest(requestState: NamesCreationRequestState) {
    return makeActionObject(SET_MAGMA_NAMES_CREATION_REQUEST, {
        status: requestState.status,
        statusMessage: requestState.statusMessage,
        response: requestState.response,
    });
}


export function makeCreateNamesCreationRequest(projectName: string, payloads: CompleteCreateNameRequestPayload[]) {
    const names = payloads.map(payload => ({
        rule_name: payload.ruleName,
        name: payload.renderedName,
    }));

    return async (dispatch: ReturnType<typeof useDispatch>) => {
        dispatch(setMagmaNamesCreationRequest({ status: 'inProgress' }));

        try {
            const response = await postNameBatchToMagma(projectName, names);

            dispatch(setMagmaNamesCreationRequest({
                status: 'success',
                response,
            }));
        } catch (error) {
            dispatch(setMagmaNamesCreationRequest({
                status: 'error',
                statusMessage: String(error)
            }));
        }
    };
}


export function setMagmaNamesListRequest(requestState: NamesListRequestState, ruleName: string) {
    return makeActionObject(SET_MAGMA_NAMES_LIST_REQUEST, {
        status: requestState.status,
        statusMessage: requestState.statusMessage,
        response: requestState.response,
        ruleName,
    });
}


export function setSearchVisibility(visible: boolean) {
    return makeActionObject(SET_SEARCH_VISIBILITY, { visible });
}


export function setReplaceVisibility(visible: boolean) {
    return makeActionObject(SET_REPLACE_VISIBILITY, { visible });
}


export function setMagmaCheckDuplicateNameRequest(createNameGroupLocalId: string, status: Status, hasDuplicate?: boolean) {
    return makeActionObject(SET_MAGMA_CHECK_DUPLICATE_NAME_REQUEST, { createNameGroupLocalId, status, hasDuplicate });
}


export function setMagmaIncrementCounterRequest(createNameGroupLocalId: string, status: Status) {
    return makeActionObject(SET_MAGMA_INCREMENT_COUNTER_REQUEST, { createNameGroupLocalId, status });
}


export type ACTION_TYPE =
    | ReturnType<typeof addNamesWithGroupsWithTokenValues>
    | ReturnType<typeof duplicateCreateNameGroups>
    | ReturnType<typeof addOrReplaceCreateNameTokenValues>
    | ReturnType<typeof deleteCreateNameTokenValue>
    | ReturnType<typeof setCreateNameRuleCounterValues>
    | ReturnType<typeof addCreateNameGroupsToSelection>
    | ReturnType<typeof removeCreateNameGroupsFromSelection>
    | ReturnType<typeof clearCreateNameGroupsSelection>
    | ReturnType<typeof setCreateNameGroupsSelectionFromSearchCriteria>
    | ReturnType<typeof setCreateNameGroupsFilterFromSearchCriteria>
    | ReturnType<typeof clearCreateNameGroupsFilter>
    | ReturnType<typeof deleteGroupsWithNames>
    | ReturnType<typeof deleteSelectedGroupsWithNames>
    | ReturnType<typeof addCreateNameGroupsToSearchCriteria>
    | ReturnType<typeof removeCreateNameGroupsFromSearchCriteria>
    | ReturnType<typeof addCreateNameGroupsToReplaceCriteria>
    | ReturnType<typeof removeCreateNameGroupsFromReplaceCriteria>
    | ReturnType<typeof setMagmaNamesCreationRequest>
    | ReturnType<typeof setMagmaNamesListRequest>
    | ReturnType<typeof setSearchVisibility>
    | ReturnType<typeof setReplaceVisibility>
    | ReturnType<typeof setMagmaCheckDuplicateNameRequest>
    | ReturnType<typeof setMagmaIncrementCounterRequest>