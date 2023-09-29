import { CreateName, CreateNameGroup, CreateNameGroupItem, CreateNameTokenValue } from '../models';
import {
    ACTION_TYPE,
    ADD_NAMES_WITH_GROUP,
    SET_TOKEN_VALUE_FOR_CREATE_NAME,
    SET_COUNTER_VALUE_FOR_CREATE_NAME,
    SET_CREATE_NAME_GROUPS_SELECTED,
    DELETE_GROUPS_WITH_NAMES,
    DELETE_SELECTED_GROUPS_WITH_NAMES,
} from '../actions/names';
import { listToIdObject } from './utils';


// all these child states have duplicate values
// but they're all refs
interface CreateNameTokenValues {
    "byCreateNameLocalId": Record<string, CreateNameTokenValue[]>
    "byTokenValueName": Record<string, CreateNameTokenValue[]>
}


interface CreateNameGroupNamesState {
    "byCreateNameLocalId": Record<string, CreateNameGroupItem>
    "byCreateNameGroupLocalId": Record<string, CreateNameGroupItem>
}


export interface NamesState {
    createNames: Record<string, CreateName>
    createNameTokenValues: CreateNameTokenValues
    createNameGroups: Record<string, CreateNameGroup>
    createNameGroupItems: CreateNameGroupNamesState
}

const initialState: NamesState = {
    createNames: {},
    createNameTokenValues: { byCreateNameLocalId: {}, byTokenValueName: {} },
    createNameGroups: {},
    createNameGroupItems: { byCreateNameLocalId: {}, byCreateNameGroupLocalId: {} }
}


export function namesReducer(state: NamesState = initialState, action: ACTION_TYPE): NamesState {
    switch (action.type) {
        case ADD_NAMES_WITH_GROUP:
            return {
                ...state,
                createNames: {
                    ...state.createNames,
                    ...listToIdObject(action.createNames, "localId"),
                },
                createNameGroups: {
                    ...state.createNameGroups,
                    [action.createNameGroup.localId]: action.createNameGroup,
                },
                createNameGroupItems: {
                    ...state.createNameGroupItems,
                    byCreateNameLocalId: {
                        ...state.createNameGroupItems.byCreateNameLocalId,
                        ...listToIdObject(action.createNameGroupItems, "createNameLocalId"),
                    },
                    byCreateNameGroupLocalId: {
                        ...state.createNameGroupItems.byCreateNameGroupLocalId,
                        ...listToIdObject(action.createNameGroupItems, "createNameLocalId"),
                    },
                }
            }
        case DELETE_GROUPS_WITH_NAMES:
            return deleteGroupsWithNames(action.createNameGroupIds, state)
        case DELETE_SELECTED_GROUPS_WITH_NAMES:
            return deleteSelectedGroupsWithNames(state)
        case SET_TOKEN_VALUE_FOR_CREATE_NAME:
            // TODOTODOTODOTODOTODOTODOTODOTODOTODOTODOTODOTODO
            const tokenValues = [...state.createNames[action.createNameLocalId].tokenValues]
            tokenValues[action.tokenIdx] = action.tokenValue

            return {
                ...state,
                createNameTokenValues: {
                    ...state.createNameTokenValues,
                    
                }
            }
        case SET_COUNTER_VALUE_FOR_CREATE_NAME:
            return {
                ...state,
                createNames: {
                    ...state.createNames,
                    [action.createNameLocalId]: {
                        ...state.createNames[action.createNameLocalId],
                        counterValue: action.counterValue,
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


function deleteGroupsWithNames(createNameGroupIds: string[], state: NamesState): NamesState {
    const newGroups = { ...state.createNameGroups }
    const newNames = { ...state.createNames }

    createNameGroupIds.forEach((cngId) => {
        delete newGroups[cngId]

        state.createNameGroups[cngId].createNameIds.forEach((cnId) => {
            delete newNames[cnId]
        })
    })

    return {
        ...state,
        createNames: newNames,
        createNameGroups: newGroups,
    }
}

function deleteSelectedGroupsWithNames(state: NamesState): NamesState {
    const newGroups = { ...state.createNameGroups }
    const newNames = { ...state.createNames }

    Object.values(state.createNameGroups).forEach((group) => {
        if (group.selected) {
            delete newGroups[group.localId]

            group.createNameIds.forEach((cnId) => {
                delete newNames[cnId]
            })
        }
    })

    return {
        ...state,
        createNames: newNames,
        createNameGroups: newGroups,
    }
}

function setGroupsSelected(createNameGroupIds: string[], selected: boolean, state: NamesState) {
    const updatedGroups = { ...state.createNameGroups }

    createNameGroupIds.forEach((cngId) => {
        updatedGroups[cngId].selected = selected
    })

    return {
        ...state,
        createNameGroups: updatedGroups
    }
}