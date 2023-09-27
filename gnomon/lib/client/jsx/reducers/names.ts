import { CreateName, CreateNameGroup } from '../models';
import { ACTION_TYPE, ADD_NAMES_WITH_GROUP } from '../actions/names';
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
        default: {
            return state;
        }
    }
}