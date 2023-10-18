import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameGroup, CreateNameTokenValue, UNSET_VALUE } from '../models';
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
} from '../actions/names';
import { listToIdObject, listToIdGroupObject, defaultDict } from '../utils/object';
import { difference, intersection } from '../utils/set'
import { SearchReplaceCriteria, createSearchReplaceCriteriaFromGroups } from '../utils/names'



// TODO: change string[] to Set<string>?
interface CreateNamesState {
    byLocalId: Record<string, CreateName>
    byCreateNameGroupLocalId: Record<string, string[]>
}


/*
    Example byRenderedTokensByParentLocalIdByCounterValue
    - parent will be UNSET_VALUE if no parent
    - counterValue will be UNSET_VALUE if no counterValue on Rule
    { 
        [UNSET_VALUE]: {
            PROJECT: {
                [UNSET_VALUE]: localId-0,
            }
        },
        ccnp-localId-5: {
            HS: {
                1: localId-1,
                2: localId-2,
            },
        },
        ccnp-localId-6: {
            DNA: {
                0: localId-4,
            }
        },
    }
*/
export type CompleteCreateNamesByParentAndValues = Record<
    string | typeof UNSET_VALUE, Record<
        string, Record<
            number | typeof UNSET_VALUE, string
        >
    >
>


interface CompleteCreateNames {
    byLocalId: Record<string, CompleteCreateName>
    byLocalIdToCreateNameLocalIds: Record<string, Set<string>>
    byCreateNameLocalId: Record<string, string>
    byParentLocalIdbyRenderedTokensByCounterValue: CompleteCreateNamesByParentAndValues
}


interface CompleteCreateNameParents {
    byLocalId: Record<string, CompleteCreateNameParent>
    byParentLocalId: Record<string, Set<string>>
    byChildLocalId: Record<string, string>
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
}


export interface NamesState {
    createNames: CreateNamesState
    completeCreateNames: CompleteCreateNames
    completeCreateNameParents: CompleteCreateNameParents
    createNameTokenValues: CreateNameTokenValuesState
    createNameGroups: CreateNameGroupsState
}

const initialState: NamesState = {
    createNames: {
        byLocalId: {},
        byCreateNameGroupLocalId: {},
    },
    completeCreateNames: {
        byLocalId: {},
        byLocalIdToCreateNameLocalIds: {},
        byCreateNameLocalId: {},
        byParentLocalIdbyRenderedTokensByCounterValue: {},
    },
    completeCreateNameParents: {
        byLocalId: {},
        byParentLocalId: {},
        byChildLocalId: {},
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
    },
}


export function namesReducer(state: NamesState = initialState, action: ACTION_TYPE): NamesState {
    switch (action.type) {
        case ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES:
            return addNamesWithGroupsAndTokensValues(action.createNames, action.createNameGroups, action.createNameTokenValues, state)
        case DELETE_GROUPS_WITH_NAMES:
            return deleteGroupsWithNames(action.createNameGroupIds, state)
        case DELETE_SELECTED_GROUPS_WITH_NAMES:
            return deleteSelectedGroupsWithNames(state)
        case ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES:
            return addCreateNameTokenValues(
                action.newCreateNameTokenValues,
                deleteCreateNameTokenValues(action.oldCreateNameTokenValueLocalIds, state),
            )
        case DELETE_CREATE_NAME_TOKEN_VALUE:
            return deleteCreateNameTokenValues([action.createNameTokenValue.localId], state)
        case SET_COUNTER_VALUE_FOR_CREATE_NAMES:
            return setRuleCounterValueForCreateNames(action.ruleCounterValuesByCreateNameLocalId, state)
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
        default: {
            return state;
        }
    }
}


function addNamesWithGroupsAndTokensValues(
    createNames: CreateName[],
    createNameGroups: CreateNameGroup[],
    createNameTokenValues: CreateNameTokenValue[],
    state: NamesState
): NamesState {

    return {
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
        createNameTokenValues: addCreateNameTokenValues(createNameTokenValues, state).createNameTokenValues,
    }
}

function deleteGroupsWithNames(createNameGroupIds: string[], state: NamesState): NamesState {
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

    const { createNameTokenValues: newNameTokenValues } = deleteCreateNameTokenValues(cntvLocalIdsToDelete, state)

    return {
        ...state,
        createNames: {
            byLocalId: newNamesById,
            byCreateNameGroupLocalId: newNamesByGroupId,
        },
        createNameTokenValues: newNameTokenValues,
        createNameGroups: {
            ...state.createNameGroups,
            byLocalId: newGroupsById,
            searchLocalIds: removeGroupsFromSearch(createNameGroupIds, state).createNameGroups.searchLocalIds,
            filterLocalIds: removeGroupsFromFilter(createNameGroupIds, state).createNameGroups.filterLocalIds,
            replaceLocalIds: removeGroupsFromReplace(createNameGroupIds, state).createNameGroups.replaceLocalIds,
            selectionLocalIds: removeGroupsFromSelection(createNameGroupIds, state).createNameGroups.selectionLocalIds,
        },
    }
}

function deleteSelectedGroupsWithNames(state: NamesState): NamesState {
    return deleteGroupsWithNames([...state.createNameGroups.selectionLocalIds], state)
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

function addCreateNameTokenValues(createNameTokenValues: CreateNameTokenValue[], state: NamesState): NamesState {
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

    return {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
            byRuleTokenLocalId: newByRuleTokenLocalId,
        }
    }
}

function deleteCreateNameTokenValues(createNameTokenValueLocalIds: string[], state: NamesState): NamesState {
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

    return {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
            byRuleTokenLocalId: newByRuleTokenLocalId,
        }
    }
}

function setRuleCounterValueForCreateNames(ruleCounterValuesByCreateNameLocalId: Record<string, number | undefined>, state: NamesState): NamesState {
    const newCreateNamesByLocalId = { ...state.createNames.byLocalId }

    for (const [cnLocalId, ruleCounterValue] of Object.entries(ruleCounterValuesByCreateNameLocalId)) {

        newCreateNamesByLocalId[cnLocalId] = {
            ...newCreateNamesByLocalId[cnLocalId],
            ruleCounterValue
        }
    }

    return {
        ...state,
        createNames: {
            ...state.createNames,
            byLocalId: newCreateNamesByLocalId,
        },
    }
}