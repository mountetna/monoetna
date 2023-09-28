import { CreateName, CreateNameGroup } from '../models';
import {
    ACTION_TYPE,
    ADD_NAMES_WITH_GROUP,
    SET_TOKEN_VALUE_FOR_CREATE_NAME,
    SET_COUNTER_VALUE_FOR_CREATE_NAME,
    SET_CREATE_NAME_GROUP_SELECTED,
    DELETE_GROUP_WITH_NAMES,
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
        case DELETE_GROUP_WITH_NAMES:
            const newGroups = { ...state.createNameGroups }
            delete newGroups[action.createNameGroupId]

            const newNames = { ...state.createNames }
            state.createNameGroups[action.createNameGroupId].createNameIds.forEach((cnId) => {
                delete newNames[cnId]
            })

            return {
                ...state,
                createNames: newNames,
                createNameGroups: newGroups,
            }
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
        case SET_CREATE_NAME_GROUP_SELECTED:
            return {
                ...state,
                createNameGroups: {
                    ...state.createNameGroups,
                    [action.createNameGroupId]: {
                        ...state.createNameGroups[action.createNameGroupId],
                        selected: action.selected,
                    }
                }
            }
        default: {
            return state;
        }
    }
}