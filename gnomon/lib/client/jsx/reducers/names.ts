import _ from 'lodash'

import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameCompleteCreateName, CreateNameGroup, CreateNameTokenValue, UNSET_VALUE } from '../models';
import {
    ACTION_TYPE,
    ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES,
    ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES,
    DELETE_CREATE_NAME_TOKEN_VALUE,
    SET_COUNTER_VALUE_FOR_CREATE_NAMES,
    DELETE_GROUPS_WITH_NAMES,
    DELETE_SELECTED_GROUPS_WITH_NAMES,
    SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA,
    SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA,
    CLEAR_CREATE_NAME_GROUPS_SELECTION,
    ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA,
    REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA,
    ADD_CREATE_NAME_GROUPS_TO_SELECTION,
    REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION,
    CLEAR_CREATE_NAME_GROUPS_FILTER,
    ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA,
    REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA,
    SET_MAGMA_NAMES_CREATION_REQUEST,
    SET_COMPOSE_ERROR_FOR_CREATE_NAME_GROUP,
    SET_MAGMA_NAMES_LIST_REQUEST,
} from '../actions/names';
import { listToIdObject, listToIdGroupObject, defaultDict } from '../utils/object';
import { difference, intersection } from '../utils/set'
import { MagmaBulkGenerateResponse, MagmaListName, SearchReplaceCriteria, createSearchReplaceCriteriaFromGroups, renderCounter, renderTokens } from '../utils/names'
import { RulesStateSliceForCompleteCreateNames } from '../selectors/global';
import { createLocalId } from '../utils/models';
import { MagmaRequestState } from '../utils/names';



interface CreateNamesState {
    byLocalId: Record<string, CreateName>
    byCreateNameGroupLocalId: Record<string, string[]>
}


interface CompleteCreateNames {
    byLocalId: Record<string, CompleteCreateName>
}


interface CreateNameCompleteCreateNames {
    byLocalId: Record<string, CreateNameCompleteCreateName>
    byCreateNameLocalId: Record<string, string>
    byCompleteCreateNameLocalId: Record<string, string[]>
}


/*
    Example CompleteCreateNameParentsByRenderedValues
    - parent will be UNSET_VALUE if no parent
    - counterValue will be UNSET_VALUE if no counterValue on Rule
    - TODO: handle multiple parents
    { 
        [UNSET_VALUE]: {
            PROJECT: {
                [UNSET_VALUE]: ccnp-localId-0,
            }
        },
        ccn-local-id-0: {
            HS: {
                1: ccnp-localId-2,
                2: ccnp-localId-3,
            },
        },
        ccn-localId-1: {
            DNA: {
                0: ccnp-localId-4,
                1: ccnp-localId-5,
            }
        },
    }
*/
export type CompleteCreateNameParentsByRenderedValues = Record<
    string | typeof UNSET_VALUE, Record<
        string, Record<
            number | typeof UNSET_VALUE, string
        >
    >
>


interface CompleteCreateNameParents {
    byLocalId: Record<string, CompleteCreateNameParent>
    byParentLocalId: Record<string, string[]>
    byChildLocalId: Record<string, string[]>
    byParentLocalIdByChildRenderedTokensByChildCounterValue: CompleteCreateNameParentsByRenderedValues
}


interface CreateNameTokenValuesState {
    byLocalId: Record<string, CreateNameTokenValue>
    byCreateNameLocalId: Record<string, string[]>
    byTokenValueLocalId: Record<string, string[]>
    byRuleTokenLocalId: Record<string, string[]>
}


interface CreateNameGroupsState {
    byLocalId: Record<string, CreateNameGroup>
    searchLocalIds: Set<string>
    replaceLocalIds: Set<string>
    selectionLocalIds: Set<string>
    filterLocalIds: Set<string>
    filterEnabled: boolean
    composeErrorsByLocalId: Record<string, boolean>
}


export type NamesCreationRequestState = MagmaRequestState<MagmaBulkGenerateResponse>


export type NamesListRequestState = MagmaRequestState<MagmaListName[]>


export interface NamesState {
    createNames: CreateNamesState
    completeCreateNames: CompleteCreateNames
    createNameCompleteCreateNames: CreateNameCompleteCreateNames
    completeCreateNameParents: CompleteCreateNameParents
    createNameTokenValues: CreateNameTokenValuesState
    createNameGroups: CreateNameGroupsState
    creationRequest: NamesCreationRequestState
    magmaNamesListRequestsByRuleName: Record<string, NamesListRequestState>
}

const initialState: NamesState = {
    createNames: {
        byLocalId: {},
        byCreateNameGroupLocalId: {},
    },
    completeCreateNames: {
        byLocalId: {},
        // byRenderedTokensByChildCounterValue: {},
    },
    createNameCompleteCreateNames: {
        byLocalId: {},
        byCreateNameLocalId: {},
        byCompleteCreateNameLocalId: {},
    },
    completeCreateNameParents: {
        byLocalId: {},
        byParentLocalId: {},
        byChildLocalId: {},
        byParentLocalIdByChildRenderedTokensByChildCounterValue: {},
    },
    createNameTokenValues: {
        byLocalId: {},
        byCreateNameLocalId: {},
        byTokenValueLocalId: {},
        byRuleTokenLocalId: {},
    },
    createNameGroups: {
        byLocalId: {},
        searchLocalIds: new Set(),
        replaceLocalIds: new Set(),
        selectionLocalIds: new Set(),
        filterLocalIds: new Set(),
        filterEnabled: false,
        composeErrorsByLocalId: {},
    },
    creationRequest: { status: "idle" },
    magmaNamesListRequestsByRuleName: {},
}


export function namesReducer(state: NamesState = initialState, action: ACTION_TYPE): NamesState {
    switch (action.type) {
        case ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES:
            return addNamesWithGroupsAndTokensValues(
                action.createNames,
                action.createNameGroups,
                action.createNameTokenValues,
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state,
            )
        case DELETE_GROUPS_WITH_NAMES:
            return deleteGroupsWithNames(
                action.createNameGroupIds,
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state,
            )
        case DELETE_SELECTED_GROUPS_WITH_NAMES:
            return deleteSelectedGroupsWithNames(
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state,
            )
        case ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES:
            return addCreateNameTokenValues(
                action.newCreateNameTokenValues,
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                deleteCreateNameTokenValues(
                    action.oldCreateNameTokenValueLocalIds,
                    false,  // don't need to update completion status multiple times
                    action.rulesStateSliceForCompleteCreateNames,
                    state
                ),
            )
        case DELETE_CREATE_NAME_TOKEN_VALUE:
            return deleteCreateNameTokenValues(
                [action.createNameTokenValue.localId],
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state
            )
        case SET_COUNTER_VALUE_FOR_CREATE_NAMES:
            return setRuleCounterValueForCreateNames(
                action.ruleCounterValuesByCreateNameLocalId,
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state
            )
        case ADD_CREATE_NAME_GROUPS_TO_SELECTION:
            return addGroupsToSelection(action.createNameGroupIds, state)
        case REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION:
            return removeGroupsFromSelection(action.createNameGroupIds, state)
        case CLEAR_CREATE_NAME_GROUPS_SELECTION:
            return deselectAllGroups(state)
        case CLEAR_CREATE_NAME_GROUPS_FILTER:
            return disableGroupFilter(state)
        case ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA:
            return addGroupsToSearch(action.createNameGroupLocalIds, state)
        case REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA:
            return removeGroupsFromSearch(action.createNameGroupLocalIds, state)
        case SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA:
            return setGroupsSelectionFromSearchCriteria(state)
        case SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA:
            return setGroupsFilterFromSearchCriteria(state)
        case ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA:
            return addGroupsToReplace(action.createNameGroupLocalIds, state)
        case REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA:
            return removeGroupsFromReplace(action.createNameGroupLocalIds, state)
        case SET_MAGMA_NAMES_CREATION_REQUEST:
            return {
                ...state,
                creationRequest: {
                    ...state.creationRequest,
                    ..._.omit(action, ["type"])
                }
            }
        case SET_COMPOSE_ERROR_FOR_CREATE_NAME_GROUP:
            return setComposeErrorForCreateNameGroup(
                action.createNameGroupLocalId,
                action.hasError,
                state,
            )
        case SET_MAGMA_NAMES_LIST_REQUEST:
            return {
                ...state,
                magmaNamesListRequestsByRuleName: {
                    ...state.magmaNamesListRequestsByRuleName,
                    [action.ruleName]: {
                        ...state.magmaNamesListRequestsByRuleName[action.ruleName] || {},
                        ..._.omit(action, ["type", "ruleName"])
                    },
                }
            }
        default: {
            return state;
        }
    }
}


function addNamesWithGroupsAndTokensValues(
    createNames: CreateName[],
    createNameGroups: CreateNameGroup[],
    createNameTokenValues: CreateNameTokenValue[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {

    let newState = {
        ...state,
        createNames: {
            byLocalId: {
                ...state.createNames.byLocalId,
                ...listToIdObject(createNames, "localId"),
            },
            byCreateNameGroupLocalId: {
                ...state.createNames.byCreateNameGroupLocalId,
                ...listToIdGroupObject(createNames, "createNameGroupLocalId", "localId"),
            },
        },
        createNameGroups: {
            ...state.createNameGroups,
            byLocalId: {
                ...state.createNameGroups.byLocalId,
                ...listToIdObject(createNameGroups, "localId"),
            }
        },
    }

    newState.createNameTokenValues = addCreateNameTokenValues(
        createNameTokenValues,
        false,
        rulesStateSliceForCompleteCreateNames,
        newState,
    ).createNameTokenValues

    const cngLocalIds = createNameGroups.map(cng => cng.localId)

    if (updateCompletionStatus) {
        newState = addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
            cngLocalIds,
            rulesStateSliceForCompleteCreateNames,
            newState,
        )
    }

    return newState
}

function deleteGroupsWithNames(
    createNameGroupIds: string[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newGroupsById = { ...state.createNameGroups.byLocalId }
    const newNamesById = { ...state.createNames.byLocalId }
    const newNamesByGroupId = { ...state.createNames.byCreateNameGroupLocalId }

    const cntvLocalIdsToDelete: string[] = []

    createNameGroupIds.forEach((cngId) => {
        delete newGroupsById[cngId]

        state.createNames.byCreateNameGroupLocalId[cngId].forEach((cnId) => {
            delete newNamesById[cnId]

            // keep track of CreateNameTokenValue.localIds to cleanup
            cntvLocalIdsToDelete.push(...(state.createNameTokenValues.byCreateNameLocalId[cnId] || []))
        })
        delete newNamesByGroupId[cngId]
    })

    let newState = { ...state }

    if (updateCompletionStatus) {
        newState = removeCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(createNameGroupIds, newState)
    }

    return {
        ...newState,
        createNames: {
            byLocalId: newNamesById,
            byCreateNameGroupLocalId: newNamesByGroupId,
        },
        createNameTokenValues: deleteCreateNameTokenValues(
            cntvLocalIdsToDelete,
            false,
            rulesStateSliceForCompleteCreateNames,
            newState,
        ).createNameTokenValues,
        createNameGroups: {
            ...newState.createNameGroups,
            byLocalId: newGroupsById,
            searchLocalIds: removeGroupsFromSearch(createNameGroupIds, newState).createNameGroups.searchLocalIds,
            filterLocalIds: removeGroupsFromFilter(createNameGroupIds, newState).createNameGroups.filterLocalIds,
            replaceLocalIds: removeGroupsFromReplace(createNameGroupIds, newState).createNameGroups.replaceLocalIds,
            selectionLocalIds: removeGroupsFromSelection(createNameGroupIds, newState).createNameGroups.selectionLocalIds,
        },
    }
}

function deleteSelectedGroupsWithNames(
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    return deleteGroupsWithNames(
        [...state.createNameGroups.selectionLocalIds],
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames,
        state
    )
}

function addGroupsToSelection(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            selectionLocalIds: new Set([...state.createNameGroups.selectionLocalIds, ...createNameGroupIds]),
        }
    }
}

function removeGroupsFromSelection(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            selectionLocalIds: difference(state.createNameGroups.selectionLocalIds, new Set(createNameGroupIds)),
        }
    }
}

function deselectAllGroups(state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            selectionLocalIds: new Set(),
        }
    }
}

// TODO add tests
// TODO add procedure breakdown
function getMatchedGroupIdsFromSearchCriteria(
    searchCriteria: SearchReplaceCriteria,
    state: NamesState,
    respectFilter: boolean = true
): string[] {

    const searchCriteriaCreateNameTokenValues = Object.values(searchCriteria.byRuleName)
        .map(searchCriteriaByRule => {
            return searchCriteriaByRule.createNameTokenValueLocalIds
                .map(cntvLocalId => state.createNameTokenValues.byLocalId[cntvLocalId])
        })
        .reduce((cntvList, cntvsToAdd) => [...cntvList, ...cntvsToAdd], [])

    const groupIdsToMatchingValuesCount = defaultDict<string, number>(_ => 0)
    const createNamesToCheckRuleCounter = new Set<CreateName>()

    searchCriteriaCreateNameTokenValues.forEach(searchCriteriaCntv => {
        const cntvIdsToCheck = state.createNameTokenValues.byTokenValueLocalId[searchCriteriaCntv.tokenValueLocalId]

        cntvIdsToCheck.forEach(cntvId => {
            const cntvToCheck = state.createNameTokenValues.byLocalId[cntvId]
            const createName = state.createNames.byLocalId[cntvToCheck.createNameLocalId]
            const ruleSearchCriteria = searchCriteria.byRuleName[createName.ruleName]

            if (
                !ruleSearchCriteria
                || (
                    respectFilter
                    && state.createNameGroups.filterEnabled
                    && !state.createNameGroups.filterLocalIds.has(createName.createNameGroupLocalId)
                )
                || state.createNameGroups.searchLocalIds.has(createName.createNameGroupLocalId)
                || state.createNameGroups.replaceLocalIds.has(createName.createNameGroupLocalId)
            ) {
                return
            }

            // match if same token value and position in rule
            if (cntvToCheck.ruleTokenLocalId == searchCriteriaCntv.ruleTokenLocalId) {
                ++groupIdsToMatchingValuesCount[createName.createNameGroupLocalId]
            }

            createNamesToCheckRuleCounter.add(createName)
        })
    })

    // match if same ruleCounterValue
    // (relies on the fact that a Rule will show up only once per CreateNameGroup)
    createNamesToCheckRuleCounter.forEach(cn => {
        const ruleSearchCriteria = searchCriteria.byRuleName[cn.ruleName]

        if (
            ruleSearchCriteria.ruleCounterValue != undefined
            && cn.ruleCounterValue == ruleSearchCriteria.ruleCounterValue
        ) {
            ++groupIdsToMatchingValuesCount[cn.createNameGroupLocalId]
        }
    })

    const targetMatchCount = Object.values(searchCriteria.byRuleName)
        .reduce((totalCount, ruleSearchCriteria): number => {

            const ruleCount = ruleSearchCriteria.createNameTokenValueLocalIds.length
                + (ruleSearchCriteria.ruleCounterValue != undefined ? 1 : 0)

            return totalCount + ruleCount
        }, 0)



    const matchingGroupIds: string[] = []

    Object.entries(groupIdsToMatchingValuesCount).forEach(([groupId, matchCount]) => {
        if (matchCount == targetMatchCount) {
            matchingGroupIds.push(groupId)
        }
    })

    return matchingGroupIds
}

function setGroupsSelectionFromSearchCriteria(state: NamesState): NamesState {
    const searchCriteriaList = createSearchReplaceCriteriaFromGroups(state, state.createNameGroups.searchLocalIds)
    state = deselectAllGroups(state)

    searchCriteriaList.forEach(searchCriteria => {
        state = addGroupsToSelection(getMatchedGroupIdsFromSearchCriteria(searchCriteria, state), state)
    })

    return state
}

function setGroupsFilterFromSearchCriteria(state: NamesState): NamesState {
    const searchCriteriaList = createSearchReplaceCriteriaFromGroups(state, state.createNameGroups.searchLocalIds)
    state = disableGroupFilter(state)

    searchCriteriaList.forEach(searchCriteria => {
        state = addGroupsToFilter(getMatchedGroupIdsFromSearchCriteria(searchCriteria, state), state)
    })

    return state
}

function addGroupsToSearch(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            searchLocalIds: new Set([...state.createNameGroups.searchLocalIds, ...createNameGroupIds]),
        }
    }
}

function removeGroupsFromSearch(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            searchLocalIds: difference(state.createNameGroups.searchLocalIds, new Set(createNameGroupIds)),
        }
    }
}

function addGroupsToReplace(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            replaceLocalIds: new Set([...state.createNameGroups.replaceLocalIds, ...createNameGroupIds]),
        }
    }
}

function removeGroupsFromReplace(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            replaceLocalIds: difference(state.createNameGroups.replaceLocalIds, new Set(createNameGroupIds)),
        }
    }
}

function addGroupsToFilter(createNameGroupIds: string[], state: NamesState): NamesState {
    // also need to remove groups from selection if not in filter
    const newFilterIds = new Set([...state.createNameGroups.filterLocalIds, ...createNameGroupIds])
    const newSelectionIds = intersection(state.createNameGroups.selectionLocalIds, newFilterIds)
    const newSelectionState = addGroupsToSelection([...newSelectionIds], deselectAllGroups(state))

    return {
        ...newSelectionState,
        createNameGroups: {
            ...newSelectionState.createNameGroups,
            filterLocalIds: newFilterIds,
            filterEnabled: true,
        }
    }
}

function removeGroupsFromFilter(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            filterLocalIds: difference(state.createNameGroups.filterLocalIds, new Set(createNameGroupIds)),
        }
    }
}

function disableGroupFilter(state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            filterLocalIds: new Set(),
            filterEnabled: false,
        }
    }
}

function addCreateNameTokenValues(
    createNameTokenValues: CreateNameTokenValue[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newByCreateNameLocalId = { ...state.createNameTokenValues.byCreateNameLocalId }
    const newByTokenValueLocalId = { ...state.createNameTokenValues.byTokenValueLocalId }
    const newByRuleTokenLocalId = { ...state.createNameTokenValues.byRuleTokenLocalId }
    const newByLocalId = { ...state.createNameTokenValues.byLocalId }

    createNameTokenValues.forEach(cntv => {
        newByLocalId[cntv.localId] = cntv
        newByCreateNameLocalId[cntv.createNameLocalId] = [
            ...newByCreateNameLocalId[cntv.createNameLocalId] || [],
            cntv.localId
        ]
        newByTokenValueLocalId[cntv.tokenValueLocalId] = [
            ...newByTokenValueLocalId[cntv.tokenValueLocalId] || [],
            cntv.localId
        ]
        newByRuleTokenLocalId[cntv.ruleTokenLocalId] = [
            ...newByRuleTokenLocalId[cntv.ruleTokenLocalId] || [],
            cntv.localId
        ]
    })

    let newState = {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
            byRuleTokenLocalId: newByRuleTokenLocalId,
        }
    }

    const cngLocalIds = new Set<string>()

    for (const cntv of createNameTokenValues) {
        cngLocalIds.add(
            state.createNames.byLocalId[cntv.createNameLocalId].createNameGroupLocalId
        )
    }

    if (updateCompletionStatus === true) {
        newState = addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
            [...cngLocalIds],
            rulesStateSliceForCompleteCreateNames,
            newState,
        )
    }

    return newState
}

function deleteCreateNameTokenValues(
    createNameTokenValueLocalIds: string[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newByCreateNameLocalId = { ...state.createNameTokenValues.byCreateNameLocalId }
    const newByTokenValueLocalId = { ...state.createNameTokenValues.byTokenValueLocalId }
    const newByRuleTokenLocalId = { ...state.createNameTokenValues.byRuleTokenLocalId }
    const newByLocalId = { ...state.createNameTokenValues.byLocalId }

    createNameTokenValueLocalIds.forEach(cntvLocalId => {
        const cntv = state.createNameTokenValues.byLocalId[cntvLocalId]

        delete newByLocalId[cntv.localId]

        if (cntv.createNameLocalId in newByCreateNameLocalId) {
            newByCreateNameLocalId[cntv.createNameLocalId] = newByCreateNameLocalId[cntv.createNameLocalId]
                .filter(cntvLocalId => cntvLocalId != cntv.localId)

            if (newByCreateNameLocalId[cntv.createNameLocalId].length == 0) {
                delete newByCreateNameLocalId[cntv.createNameLocalId]
            }
        }
        if (cntv.tokenValueLocalId in newByTokenValueLocalId) {
            newByTokenValueLocalId[cntv.tokenValueLocalId] = newByTokenValueLocalId[cntv.tokenValueLocalId]
                .filter(cntvLocalId => cntvLocalId != cntv.localId)

            if (newByTokenValueLocalId[cntv.tokenValueLocalId].length == 0) {
                delete newByTokenValueLocalId[cntv.tokenValueLocalId]
            }
        }
        if (cntv.ruleTokenLocalId in newByRuleTokenLocalId) {
            newByRuleTokenLocalId[cntv.ruleTokenLocalId] = newByRuleTokenLocalId[cntv.ruleTokenLocalId]
                .filter(cntvLocalId => cntvLocalId != cntv.localId)

            if (newByRuleTokenLocalId[cntv.ruleTokenLocalId].length == 0) {
                delete newByRuleTokenLocalId[cntv.ruleTokenLocalId]
            }
        }
    })

    let newState = {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
            byRuleTokenLocalId: newByRuleTokenLocalId,
        }
    }

    const cngLocalIds = new Set<string>()

    for (const cntvLocalId of createNameTokenValueLocalIds) {
        const cnLocalId = state.createNameTokenValues.byLocalId[cntvLocalId].createNameLocalId

        cngLocalIds.add(
            state.createNames.byLocalId[cnLocalId].createNameGroupLocalId
        )
    }

    if (updateCompletionStatus === true) {
        newState = addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
            [...cngLocalIds],
            rulesStateSliceForCompleteCreateNames,
            newState,
        )
    }

    return newState
}

function setRuleCounterValueForCreateNames(
    ruleCounterValuesByCreateNameLocalId: Record<string, number | undefined>,
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newCreateNamesByLocalId = { ...state.createNames.byLocalId }

    for (const [cnLocalId, ruleCounterValue] of Object.entries(ruleCounterValuesByCreateNameLocalId)) {

        newCreateNamesByLocalId[cnLocalId] = {
            ...newCreateNamesByLocalId[cnLocalId],
            ruleCounterValue
        }
    }

    let newState = {
        ...state,
        createNames: {
            ...state.createNames,
            byLocalId: newCreateNamesByLocalId,
        },
    }

    const cngLocalIds = new Set<string>()

    for (const cnLocalId of Object.keys(ruleCounterValuesByCreateNameLocalId)) {
        cngLocalIds.add(
            state.createNames.byLocalId[cnLocalId].createNameGroupLocalId
        )
    }

    if (updateCompletionStatus === true) {
        newState = addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
            [...cngLocalIds],
            rulesStateSliceForCompleteCreateNames,
            newState,
        )
    }

    return newState
}

function addCompleteCreateNames(
    completeCreateNames: CompleteCreateName[],
    createNameCompleteCreateNames: CreateNameCompleteCreateName[],
    completeCreateNameParents: CompleteCreateNameParent[],
    state: NamesState,
): NamesState {
    const newCompleteCreateNames = {
        byLocalId: {
            ...state.completeCreateNames.byLocalId,
            ...listToIdObject(completeCreateNames, "localId"),
        }
    }

    const newCreateNameCompleteCreateNames = {
        byLocalId: { ...state.createNameCompleteCreateNames.byLocalId },
        byCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCreateNameLocalId },
        byCompleteCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCompleteCreateNameLocalId },
    }

    for (const cnccn of createNameCompleteCreateNames) {
        newCreateNameCompleteCreateNames.byLocalId[cnccn.localId] = cnccn
        newCreateNameCompleteCreateNames.byCreateNameLocalId[cnccn.createNameLocalId] = cnccn.localId

        if (!(cnccn.completeCreateNameLocalId in newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId)) {
            newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId] = []
        }
        newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId].push(cnccn.localId)
    }

    const newCompleteCreateNameParents = {
        byLocalId: {
            ...state.completeCreateNameParents.byLocalId,
            ...listToIdObject(completeCreateNameParents, "localId"),
        },
        byParentLocalId: { ...state.completeCreateNameParents.byParentLocalId },
        byChildLocalId: { ...state.completeCreateNameParents.byChildLocalId },
        byParentLocalIdByChildRenderedTokensByChildCounterValue: _.cloneDeep(
            state.completeCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
        ),
    }

    for (const ccnp of completeCreateNameParents) {
        if (!(ccnp.parentCompleteCreateNameLocalId in newCompleteCreateNameParents.byParentLocalId)) {
            newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId] = []
        }
        newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId].push(ccnp.localId)

        if (!(ccnp.completeCreateNameLocalId in newCompleteCreateNameParents.byChildLocalId)) {
            newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId] = []
        }
        newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId].push(ccnp.localId)


        const childCompleteCreateName = state.completeCreateNames.byLocalId[ccnp.completeCreateNameLocalId] = newCompleteCreateNames.byLocalId[ccnp.completeCreateNameLocalId]

        const hierarchyCounterValue = childCompleteCreateName.counterValue != undefined ? childCompleteCreateName.counterValue : UNSET_VALUE

        if (!(ccnp.parentCompleteCreateNameLocalId in newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue)) {
            newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue[ccnp.parentCompleteCreateNameLocalId] = {}
        }

        const byParentLocalId = newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue[ccnp.parentCompleteCreateNameLocalId]

        if (!(childCompleteCreateName.value in byParentLocalId)) {
            // @ts-ignore
            byParentLocalId[childCompleteCreateName.value] = {}
        }
        byParentLocalId[childCompleteCreateName.value][hierarchyCounterValue] = ccnp.localId
    }

    return {
        ...state,
        completeCreateNames: newCompleteCreateNames,
        createNameCompleteCreateNames: newCreateNameCompleteCreateNames,
        completeCreateNameParents: newCompleteCreateNameParents,
    }
}

function removeCompleteCreateNames(
    createNameCompleteCreateNames: CreateNameCompleteCreateName[],
    state: NamesState,
): NamesState {
    // find CompleteCreateNames and Parents from CreateNameCompleteCreateName
    const completeCreateNames: CompleteCreateName[] = []
    const completeCreateNameParents = new Set<CompleteCreateNameParent>()

    const newCreateNameCompleteCreateNames = {
        byLocalId: { ...state.createNameCompleteCreateNames.byLocalId },
        byCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCreateNameLocalId },
        byCompleteCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCompleteCreateNameLocalId },
    }

    for (const cnccn of createNameCompleteCreateNames) {
        delete newCreateNameCompleteCreateNames.byLocalId[cnccn.localId]
        delete newCreateNameCompleteCreateNames.byCreateNameLocalId[cnccn.createNameLocalId]

        if (cnccn.completeCreateNameLocalId in newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId) {

            newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId] =
                newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId]
                    .filter(candidateCnccnLocalId => candidateCnccnLocalId != cnccn.localId)

            if (newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId].length == 0) {
                delete newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId]

                // and delete CompleteCreateName and Parents since there aren't any CreateName associations left
                completeCreateNames.push(state.completeCreateNames.byLocalId[cnccn.completeCreateNameLocalId])

                const parentsToDelete = [
                    ...(state.completeCreateNameParents.byParentLocalId[cnccn.completeCreateNameLocalId] || []),
                    ...(state.completeCreateNameParents.byChildLocalId[cnccn.completeCreateNameLocalId] || []),
                ]
                parentsToDelete.forEach(ccnpLocalId => {
                    completeCreateNameParents.add(state.completeCreateNameParents.byLocalId[ccnpLocalId])
                })
            }
        }
    }

    const newCompleteCreateNames = {
        byLocalId: { ...state.completeCreateNames.byLocalId },
    }

    for (const ccn of completeCreateNames) {
        delete newCompleteCreateNames.byLocalId[ccn.localId]
    }

    const newCompleteCreateNameParents = {
        byLocalId: { ...state.completeCreateNameParents.byLocalId },
        byParentLocalId: { ...state.completeCreateNameParents.byParentLocalId },
        byChildLocalId: { ...state.completeCreateNameParents.byChildLocalId },
        byParentLocalIdByChildRenderedTokensByChildCounterValue: _.cloneDeep(
            state.completeCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
        ),
    }

    for (const ccnp of completeCreateNameParents) {
        delete newCompleteCreateNameParents.byLocalId[ccnp.localId]

        if (ccnp.parentCompleteCreateNameLocalId in newCompleteCreateNameParents.byParentLocalId) {

            newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId] =
                newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId]
                    .filter(candidateCcnpLocalId => candidateCcnpLocalId != ccnp.localId)

            if (newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId].length == 0) {
                delete newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId]
            }
        }

        if (ccnp.completeCreateNameLocalId in newCompleteCreateNameParents.byChildLocalId) {

            newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId] =
                newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId]
                    .filter(candidateCcnpLocalId => candidateCcnpLocalId != ccnp.localId)

            if (newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId].length == 0) {
                delete newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId]
            }
        }

        // remove from byParentLocalIdByChildRenderedTokensByChildCounterValue
        const childCompleteCreateName = state.completeCreateNames.byLocalId[ccnp.completeCreateNameLocalId]

        delete newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
        [ccnp.parentCompleteCreateNameLocalId][childCompleteCreateName.value]
        [childCompleteCreateName.counterValue != undefined ? childCompleteCreateName.counterValue : UNSET_VALUE]

        if (
            Object.keys(newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId][childCompleteCreateName.value]).length == 0
        ) {
            delete newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId][childCompleteCreateName.value]
        }

        if (
            Object.keys(newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId]).length == 0
        ) {
            delete newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId]
        }
    }

    return {
        ...state,
        completeCreateNames: newCompleteCreateNames,
        createNameCompleteCreateNames: newCreateNameCompleteCreateNames,
        completeCreateNameParents: newCompleteCreateNameParents,
    }
}

function isReplaceOrSearchGroup(createNameGroupLocalId: string, state: NamesState): boolean {
    return (
        state.createNameGroups.replaceLocalIds.has(createNameGroupLocalId)
        || state.createNameGroups.searchLocalIds.has(createNameGroupLocalId)
    )
}

// TODO: handle multiparent rules
// TODO: ID should be monotonically increasing
function addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
    createNameGroupsLocalIds: string[],
    rulesStateSlice: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    let newState = { ...state }

    // for each CreateNameGroup
    for (const cngLocalId of createNameGroupsLocalIds) {
        const completeCreateNamesToAdd: CompleteCreateName[] = []
        const createNameCompleteCreateNamesToAdd: CreateNameCompleteCreateName[] = []
        const createNameCompleteCreateNamesToRemove: CreateNameCompleteCreateName[] = []
        const completeCreateNameParentsToAdd: CompleteCreateNameParent[] = []

        if (isReplaceOrSearchGroup(cngLocalId, newState)) {
            continue
        }

        const createNameGroup = newState.createNameGroups.byLocalId[cngLocalId]
        const primaryRuleName = newState.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId].ruleName
        const ruleNamesHierarchicalList = rulesStateSlice.ruleNamesHierarchicalListByPrimaryRuleName[primaryRuleName]

        const createNames = newState.createNames.byCreateNameGroupLocalId[cngLocalId]
            .map(cnLocalId => newState.createNames.byLocalId[cnLocalId])

        const sortedCreateNames: CreateName[] = _.sortBy(createNames, [(cn: CreateName) => {
            return ruleNamesHierarchicalList.indexOf(cn.ruleName)
        }])

        // keep track of hierarchy completeness
        const completeCreateNames: (CompleteCreateName | undefined)[] = []
        let completeHierarchy = true

        // for each CreateName
        for (const [idx, createName] of sortedCreateNames.entries()) {
            const createNameTokenValues = (newState.createNameTokenValues.byCreateNameLocalId[createName.localId] || [])
                .map(cntvLocalId => newState.createNameTokenValues.byLocalId[cntvLocalId])

            const actualTokenCount = createNameTokenValues.length
            const targetTokenCount = rulesStateSlice.ruleTokenLocalIdsByRuleName[createName.ruleName].length
            const counterRequired = rulesStateSlice.counterRequiredByRuleName[createName.ruleName]
            const hasCounterValue = createName.ruleCounterValue != undefined

            let parentCompleteCreateName: CompleteCreateName | undefined = undefined

            if (idx > 0) {
                parentCompleteCreateName = completeCreateNames[idx - 1]

                if (parentCompleteCreateName == undefined) {
                    completeHierarchy = false
                }
            }

            // detect completeness
            if (
                targetTokenCount === actualTokenCount
                && (
                    !counterRequired
                    || hasCounterValue
                )
                && completeHierarchy
            ) {
                // is complete (at this level)

                let completeCreateNameParentLocalId: string | undefined = undefined
                const hierarchyRenderedTokenValue = renderTokens(
                    createNameTokenValues,
                    rulesStateSlice.ruleTokensByLocalId,
                    rulesStateSlice.tokenValuesByLocalId,
                )
                const hierarchyRenderedCounterValue = Number.parseInt(createName.ruleCounterValue != undefined ? renderCounter(createName) : UNSET_VALUE)

                // find completeCreateNameParentLocalId
                if (parentCompleteCreateName != undefined || idx == 0) {
                    const parentCompleteCreateNameLocalId = parentCompleteCreateName != undefined ? parentCompleteCreateName.localId : UNSET_VALUE

                    const renderedTokenValuesDict = newState.completeCreateNameParents
                        .byParentLocalIdByChildRenderedTokensByChildCounterValue[parentCompleteCreateNameLocalId]

                    if (renderedTokenValuesDict != undefined) {
                        const renderedCounterValuesDict = renderedTokenValuesDict[hierarchyRenderedTokenValue]

                        if (renderedCounterValuesDict != undefined) {
                            completeCreateNameParentLocalId = renderedCounterValuesDict[hierarchyRenderedCounterValue]
                        }
                    }
                }

                let completeCreateNameParent: CompleteCreateNameParent | undefined = undefined

                if (completeCreateNameParentLocalId) {
                    completeCreateNameParent = newState.completeCreateNameParents.byLocalId[completeCreateNameParentLocalId]
                }


                // create CompleteCreateName and CreateName association if doesn't exist
                let completeCreateName: CompleteCreateName | undefined = undefined

                if (completeCreateNameParent == undefined) {
                    completeCreateName = {
                        localId: createLocalId(),
                        value: renderTokens(
                            createNameTokenValues,
                            rulesStateSlice.ruleTokensByLocalId,
                            rulesStateSlice.tokenValuesByLocalId,
                        ),
                        counterValue: createName.ruleCounterValue,
                    }

                    completeCreateNamesToAdd.push(completeCreateName)

                    // remove existing association and its parents if exist
                    const existingAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId]
                    const existingAssociation = newState.createNameCompleteCreateNames.byLocalId[existingAssociationLocalId]

                    if (existingAssociation) {
                        createNameCompleteCreateNamesToRemove.push(existingAssociation)
                    }

                    // create new association
                    const completeCreateNameAssociation: CreateNameCompleteCreateName = {
                        localId: createLocalId(),
                        createNameLocalId: createName.localId,
                        completeCreateNameLocalId: completeCreateName.localId,
                    }
                    createNameCompleteCreateNamesToAdd.push(completeCreateNameAssociation)
                } else {
                    completeCreateName = newState.completeCreateNames.byLocalId[completeCreateNameParent.completeCreateNameLocalId]
                    const existingAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId]
                    const existingAssociation = newState.createNameCompleteCreateNames.byLocalId[existingAssociationLocalId]

                    if (existingAssociation?.completeCreateNameLocalId != completeCreateName.localId) {
                        // create association
                        const completeCreateNameAssociation: CreateNameCompleteCreateName = {
                            localId: createLocalId(),
                            createNameLocalId: createName.localId,
                            completeCreateNameLocalId: completeCreateName.localId,
                        }
                        createNameCompleteCreateNamesToAdd.push(completeCreateNameAssociation)

                        if (existingAssociation) {
                            createNameCompleteCreateNamesToRemove.push(existingAssociation)
                        }
                    }
                }

                // create CompleteCreateNameParent if doesn't exist
                if (completeCreateNameParent == undefined) {
                    const parent: CompleteCreateNameParent = {
                        localId: createLocalId(),
                        completeCreateNameLocalId: completeCreateName.localId,
                        parentCompleteCreateNameLocalId: parentCompleteCreateName?.localId || UNSET_VALUE,
                    }

                    completeCreateNameParentsToAdd.push(parent)
                }

                completeCreateNames.push(completeCreateName)

            } else {
                // is not complete (at this level)
                // remove all records if exist

                const completeCreateNameAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId]

                if (completeCreateNameAssociationLocalId != undefined) {
                    const completeCreateNameAssociation = newState.createNameCompleteCreateNames.byLocalId[completeCreateNameAssociationLocalId]
                    createNameCompleteCreateNamesToRemove.push(completeCreateNameAssociation)
                }

                completeCreateNames.push(undefined)
            }
        }

        newState = addCompleteCreateNames(
            completeCreateNamesToAdd,
            createNameCompleteCreateNamesToAdd,
            completeCreateNameParentsToAdd,
            removeCompleteCreateNames(createNameCompleteCreateNamesToRemove, newState),
        )
    }

    return newState
}

// TODO: handle multiparent rules
function removeCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
    createNameGroupsLocalIds: string[],
    state: NamesState,
): NamesState {
    let newState = { ...state }

    // for each CreateNameGroup
    for (const cngLocalId of createNameGroupsLocalIds) {
        const createNameCompleteCreateNamesToRemove: CreateNameCompleteCreateName[] = []

        if (isReplaceOrSearchGroup(cngLocalId, newState)) {
            continue
        }

        const createNames = newState.createNames.byCreateNameGroupLocalId[cngLocalId]
            .map(cnLocalId => newState.createNames.byLocalId[cnLocalId])

        for (const createName of createNames) {

            const completeCreateNameAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId]
            const completeCreateNameAssociation = newState.createNameCompleteCreateNames.byLocalId[completeCreateNameAssociationLocalId]

            if (completeCreateNameAssociation) {
                createNameCompleteCreateNamesToRemove.push(completeCreateNameAssociation)
            }
        }

        newState = removeCompleteCreateNames(createNameCompleteCreateNamesToRemove, newState)
    }

    return newState
}


function setComposeErrorForCreateNameGroup(createNameGroupLocalId: string, hasError: boolean, state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            composeErrorsByLocalId: {
                ...state.createNameGroups.composeErrorsByLocalId,
                [createNameGroupLocalId]: hasError,
            },
        },
    }
}