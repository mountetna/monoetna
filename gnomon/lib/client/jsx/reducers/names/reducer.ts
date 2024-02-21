import _ from 'lodash';

import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameCompleteCreateName, CreateNameGroup, CreateNameTokenValue, UNSET_VALUE } from '../../models';
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
    SET_MAGMA_NAMES_LIST_REQUEST,
    SET_SEARCH_VISIBILITY,
    SET_REPLACE_VISIBILITY,
    SET_MAGMA_CHECK_DUPLICATE_NAME_REQUEST,
    SET_MAGMA_INCREMENT_COUNTER_REQUEST,
} from '../../actions/names';
import { MagmaBulkGenerateResponse, MagmaListName, MagmaRequestState } from '../../utils/names';
import { Status } from '../../utils/models';
import { addGroupsToReplace, addGroupsToSearch, addGroupsToSelection, addNamesWithGroupsAndTokensValues, deleteGroupsWithNames, deleteSelectedGroupsWithNames, deselectAllGroups, disableGroupFilter, removeGroupsFromReplace, removeGroupsFromSearch, removeGroupsFromSelection, setGroupsFilterFromSearchCriteria, setGroupsSelectionFromSearchCriteria, setMagmaCheckDuplicateNameRequestForCreateNameGroup, setMagmaIncrementCounterRequestForCreateNameGroup, setRuleCounterValueForCreateNames } from './create-names-and-groups';
import { addCreateNameTokenValues, deleteCreateNameTokenValues } from './create-name-token-values';



interface CreateNamesState {
    byLocalId: Record<string, CreateName>
    byCreateNameGroupLocalId: Record<string, string[]>
    byRuleName: Record<string, string[]>
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


export interface MagmaCheckDuplicateNameRequestState {
    status: Status
    hasDuplicate?: boolean
}


export interface MagmaIncrementCounterRequestState {
    status: Status
}


interface CreateNameGroupsState {
    byLocalId: Record<string, CreateNameGroup>
    searchLocalIds: Set<string>
    searchVisible: boolean
    replaceLocalIds: Set<string>
    replaceVisible: boolean
    selectionLocalIds: Set<string>
    filterLocalIds: Set<string>
    filterEnabled: boolean
    magmaCheckDuplicateNameRequestsByLocalId: Record<string, MagmaCheckDuplicateNameRequestState>
    magmaIncrementCounterRequestsByLocalId: Record<string, MagmaIncrementCounterRequestState>
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
        byRuleName: {},
    },
    completeCreateNames: {
        byLocalId: {},
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
        searchVisible: false,
        replaceLocalIds: new Set(),
        replaceVisible: false,
        selectionLocalIds: new Set(),
        filterLocalIds: new Set(),
        filterEnabled: false,
        magmaCheckDuplicateNameRequestsByLocalId: {},
        magmaIncrementCounterRequestsByLocalId: {},
    },
    creationRequest: { status: 'idle' },
    magmaNamesListRequestsByRuleName: {},
};


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
            );
        case DELETE_GROUPS_WITH_NAMES:
            return deleteGroupsWithNames(
                action.createNameGroupIds,
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state,
            );
        case DELETE_SELECTED_GROUPS_WITH_NAMES:
            return deleteSelectedGroupsWithNames(
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state,
            );
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
            );
        case DELETE_CREATE_NAME_TOKEN_VALUE:
            return deleteCreateNameTokenValues(
                [action.createNameTokenValue.localId],
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state
            );
        case SET_COUNTER_VALUE_FOR_CREATE_NAMES:
            return setRuleCounterValueForCreateNames(
                action.ruleCounterValuesByCreateNameLocalId,
                action.updateCompletionStatus,
                action.rulesStateSliceForCompleteCreateNames,
                state
            );
        case ADD_CREATE_NAME_GROUPS_TO_SELECTION:
            return addGroupsToSelection(action.createNameGroupIds, state);
        case REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION:
            return removeGroupsFromSelection(action.createNameGroupIds, state);
        case CLEAR_CREATE_NAME_GROUPS_SELECTION:
            return deselectAllGroups(state);
        case CLEAR_CREATE_NAME_GROUPS_FILTER:
            return disableGroupFilter(state);
        case ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA:
            return addGroupsToSearch(action.createNameGroupLocalIds, state);
        case REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA:
            return removeGroupsFromSearch(action.createNameGroupLocalIds, state);
        case SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA:
            return setGroupsSelectionFromSearchCriteria(state);
        case SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA:
            return setGroupsFilterFromSearchCriteria(state);
        case ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA:
            return addGroupsToReplace(action.createNameGroupLocalIds, state);
        case REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA:
            return removeGroupsFromReplace(action.createNameGroupLocalIds, state);
        case SET_MAGMA_NAMES_CREATION_REQUEST:
            return {
                ...state,
                creationRequest: {
                    ...state.creationRequest,
                    ..._.omit(action, ['type'])
                }
            };
        case SET_MAGMA_CHECK_DUPLICATE_NAME_REQUEST:
            return setMagmaCheckDuplicateNameRequestForCreateNameGroup(
                action.createNameGroupLocalId,
                action.status,
                state,
                action.hasDuplicate,
            );
        case SET_MAGMA_NAMES_LIST_REQUEST:
            return {
                ...state,
                magmaNamesListRequestsByRuleName: {
                    ...state.magmaNamesListRequestsByRuleName,
                    [action.ruleName]: {
                        ...state.magmaNamesListRequestsByRuleName[action.ruleName] || {},
                        ..._.omit(action, ['type', 'ruleName'])
                    },
                }
            };
        case SET_MAGMA_INCREMENT_COUNTER_REQUEST:
            return setMagmaIncrementCounterRequestForCreateNameGroup(
                action.createNameGroupLocalId,
                action.status,
                state,
            );
        case SET_SEARCH_VISIBILITY:
            return {
                ...state,
                createNameGroups: {
                    ...state.createNameGroups,
                    searchVisible: action.visible,
                },
            };
        case SET_REPLACE_VISIBILITY:
            return {
                ...state,
                createNameGroups: {
                    ...state.createNameGroups,
                    replaceVisible: action.visible,
                },
            };
        default: {
            return state;
        }
    }
}