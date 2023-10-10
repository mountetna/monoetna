import { CreateName, CreateNameGroup, CreateNameTokenValue, TokenValue } from '../models';
import {
    ACTION_TYPE,
    ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES,
    ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE,
    DELETE_CREATE_NAME_TOKEN_VALUE,
    SET_COUNTER_VALUE_FOR_CREATE_NAME,
    SET_CREATE_NAME_GROUPS_SELECTED,
    DELETE_GROUPS_WITH_NAMES,
    DELETE_SELECTED_GROUPS_WITH_NAMES,
} from '../actions/names';
import { listToIdObject, listToIdGroupObject } from '../utils/object';
import { selectSelectedCreateNameGroupIds } from '../selectors/names';


// TODO: change string[] to Set<string>?
interface CreateNamesState {
    byLocalId: Record<string, CreateName>
    byCreateNameGroupLocalId: Record<string, string[]>
}


interface CreateNameTokenValuesState {
    byLocalId: Record<string, CreateNameTokenValue>
    byCreateNameLocalId: Record<string, string[]>
    byTokenValueLocalId: Record<string, string[]>
}


export interface NamesState {
    createNames: CreateNamesState
    createNameTokenValues: CreateNameTokenValuesState
    createNameGroups: Record<string, CreateNameGroup>
}

const initialState: NamesState = {
    createNames: { byLocalId: {}, byCreateNameGroupLocalId: {} },
    createNameTokenValues: { byLocalId: {}, byCreateNameLocalId: {}, byTokenValueLocalId: {} },
    createNameGroups: {},
}


export function namesReducer(state: NamesState = initialState, action: ACTION_TYPE): NamesState {
    switch (action.type) {
        case ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES:
            return addNamesWithGroupsAndTokensValues(action.createNames, action.createNameGroups, action.createNameTokenValues, state)
        case DELETE_GROUPS_WITH_NAMES:
            return deleteGroupsWithNames(action.createNameGroupIds, state)
        case DELETE_SELECTED_GROUPS_WITH_NAMES:
            return deleteSelectedGroupsWithNames(state)
        case ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE:
            if (action.oldCreateNameTokenValue) {
                state = deleteCreateNameTokenValues([action.oldCreateNameTokenValue.localId], state)
            }
            return addCreateNameTokenValues([action.newCreateNameTokenValue], state)
        case DELETE_CREATE_NAME_TOKEN_VALUE:
            return deleteCreateNameTokenValues([action.createNameTokenValue.localId], state)
        case SET_COUNTER_VALUE_FOR_CREATE_NAME:
            return {
                ...state,
                createNames: {
                    ...state.createNames,
                    byLocalId: {
                        ...state.createNames.byLocalId,
                        [action.createNameLocalId]: {
                            ...state.createNames.byLocalId[action.createNameLocalId],
                            ruleCounterValue: action.ruleCounterValue,
                        }
                    }
                },
            }
        case SET_CREATE_NAME_GROUPS_SELECTED:
            return setGroupsSelected(action.createNameGroupIds, action.selected, state)
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
            ...listToIdObject(createNameGroups, "localId"),
        },
        createNameTokenValues: addCreateNameTokenValues(createNameTokenValues, state).createNameTokenValues,
    }
}


function deleteGroupsWithNames(createNameGroupIds: string[], state: NamesState): NamesState {
    const newGroups = { ...state.createNameGroups }
    const newNames = { ...state.createNames }

    const cntvLocalIdsToDelete: string[] = []

    createNameGroupIds.forEach((cngId) => {
        delete newGroups[cngId]

        state.createNames.byCreateNameGroupLocalId[cngId].forEach((cnId) => {
            delete newNames.byLocalId[cnId]

            // keep track of CreateNameTokenValue.localIds to cleanup
            cntvLocalIdsToDelete.push(...(state.createNameTokenValues.byCreateNameLocalId[cnId] || []))
        })
        delete newNames.byCreateNameGroupLocalId[cngId]
    })

    const { createNameTokenValues: newNameTokenValues } = deleteCreateNameTokenValues(cntvLocalIdsToDelete, state)

    return {
        ...state,
        createNames: newNames,
        createNameTokenValues: newNameTokenValues,
        createNameGroups: newGroups,
    }
}

function deleteSelectedGroupsWithNames(state: NamesState): NamesState {
    const selectedCreateNameGroupIds: string[] = selectSelectedCreateNameGroupIds({ names: state })

    return deleteGroupsWithNames(selectedCreateNameGroupIds, state)
}

function setGroupsSelected(createNameGroupIds: string[], selected: boolean, state: NamesState): NamesState {
    const updatedGroups = { ...state.createNameGroups }

    createNameGroupIds.forEach((cngId) => {
        updatedGroups[cngId].selected = selected
    })

    return {
        ...state,
        createNameGroups: updatedGroups
    }
}

function addCreateNameTokenValues(createNameTokenValues: CreateNameTokenValue[], state: NamesState): NamesState {
    const newByCreateNameLocalId = { ...state.createNameTokenValues.byCreateNameLocalId }
    const newByTokenValueLocalId = { ...state.createNameTokenValues.byTokenValueLocalId }
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
    })

    return {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
        }
    }
}

function deleteCreateNameTokenValues(createNameTokenValueLocalIds: string[], state: NamesState): NamesState {
    const newByCreateNameLocalId = { ...state.createNameTokenValues.byCreateNameLocalId }
    const newByTokenValueLocalId = { ...state.createNameTokenValues.byTokenValueLocalId }
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
    })

    return {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
        }
    }
}