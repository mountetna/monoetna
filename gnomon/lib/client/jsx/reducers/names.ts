import { CreateName, CreateNameGroup } from '../models';
import {
    ACTION_TYPE,
    ADD_NAMES_WITH_GROUP,
    SET_TOKEN_VALUE_FOR_CREATE_NAME,
    SET_COUNTER_VALUE_FOR_CREATE_NAME,
    SET_CREATE_NAME_GROUPS_SELECTED,
    DELETE_GROUPS_WITH_NAMES,
} from '../actions/names';
import { listToIdObject } from './utils';


export interface NamesState {
    createNames: Record<string, CreateName>
    createNameGroups: Record<string, CreateNameGroup>
}

const initialState: NamesState = {
    createNames: {},
    createNameGroups: {},
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
                }
            }
        case DELETE_GROUPS_WITH_NAMES:
            return deleteGroupsWithNames(action.createNameGroupIds, state)
        case SET_TOKEN_VALUE_FOR_CREATE_NAME:
            const tokenValues = [...state.createNames[action.createNameLocalId].tokenValues]
            tokenValues[action.tokenIdx] = action.tokenValue

            return {
                ...state,
                createNames: {
                    ...state.createNames,
                    [action.createNameLocalId]: {
                        ...state.createNames[action.createNameLocalId],
                        tokenValues,
                    }
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