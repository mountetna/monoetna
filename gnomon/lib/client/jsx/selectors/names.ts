import { createSelector } from 'reselect'
import * as _ from "lodash"

import { CreateName, CreateNameGroup } from "../models";
import { NamesState } from "../reducers/names"
import { defaultDict } from "../utils/object"


interface State {
    names: NamesState
}


export const selectCreateNameGroups = (state: State): Record<string, CreateNameGroup> => state.names.createNameGroups

export const selectCreateNameGroupById = (localId: string) => (state: State): CreateNameGroup => {
    return state.names.createNameGroups[localId]
}

export const selectCreateNameGroupByIds = (localIds: string[]) => (state: State): CreateNameGroup[] => {
    return localIds.map(localId => state.names.createNameGroups[localId])
}

export const selectCreateNameGroupIdsByPrimaryRule = createSelector(
    [(state: State) => state.names],
    (names: NamesState) => {
        const result = defaultDict<string, string[]>(_ => [])

        Object.values(names.createNameGroups).forEach(createNameGroup => {
            const createName: CreateName = names.createNames[createNameGroup.primaryCreateNameId]
            result[createName.ruleName].push(createNameGroup.localId)
        })

        return result
    }
)

export const selectCreateNames = (state: State): Record<string, CreateName> => state.names.createNames

export const selectCreateNameById = (localId: string) => (state: State): CreateName => {
    return state.names.createNames[localId]
}

export const selectCreateNamesByIds = (localIds: string[]) => (state: State): CreateName[] => {
    return localIds.map(localId => state.names.createNames[localId])
}

export const selectCounterValuesByRuleName = createSelector(
    [(state: State) => state.names],
    (names: NamesState): Record<string, Record<string, number>> => {
        const counterValuesByRuleName = defaultDict<string, Record<string, number>>(
            () => defaultDict<string, number>(() => 0)
        )

        Object.values(names.createNames).forEach((createName) => {
            if (createName.counterValue != undefined) {
                ++counterValuesByRuleName[createName.ruleName][createName.counterValue]
            }
        })

        return counterValuesByRuleName
    }
)

export const selectSelectedCreateNameGroupIds: ((state: State) => string[]) = createSelector(
    [(state: State) => state.names],
    (names: NamesState): string[] => {
        return _.filter(names.createNameGroups, { selected: true })
    }
)