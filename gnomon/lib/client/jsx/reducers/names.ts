import { CreateName, CreateNameGroup } from '../models';
import { ACTION_TYPE, ADD_NAMES_WITH_GROUP, SET_TOKEN_VALUE_FOR_CREATE_NAME } from '../actions/names';
import { listToIdObject } from './utils';



export interface NamesState {
    createNames: Record<string, CreateName>
    createNameGroups: Record<string, CreateNameGroup>
}

const initialState: NamesState = {
    createNames: {},
    createNameGroups: {}
}


export function namesReducer(state: NamesState = initialState, action: ACTION_TYPE): NamesState {
    switch (action.type) {
        case ADD_NAMES_WITH_GROUP:
            return {
                ...state,
                createNames: {
                    ...state.createNames,
                    ...listToIdObject(action.createNames, "localId")
                },
                createNameGroups: {
                    ...state.createNameGroups,
                    [action.createNameGroup.localId]: action.createNameGroup
                }
            }
        case SET_TOKEN_VALUE_FOR_CREATE_NAME:
            const createName = state.createNames[action.createNameLocalId]
            const tokenValues = [...createName.tokenValues]

            tokenValues[action.tokenIdx] = action.tokenValue

            return {
                ...state,
                createNames: {
                    ...state.createNames,
                    [action.createNameLocalId]: {
                        ...createName,
                        tokenValues
                    }
                }
            }
        default: {
            return state;
        }
    }
}