import { createSelector } from 'reselect'
import * as _ from "lodash"

import { CreateName, CreateNameGroup, CreateNameTokenValue } from "../models";
import { NamesState } from "../reducers/names"
import { defaultDict } from "../utils/object"


interface State {
    names: NamesState
}


export const selectCreateNameGroups = (state: State): Record<string, CreateNameGroup> => state.names.createNameGroups

export const selectCreateNameGroupWithLocalId = (localId: string) => (state: State): CreateNameGroup => {
    return state.names.createNameGroups[localId]
}

export const selectCreateNameGroupsWithLocalIds = (localIds: string[]) => (state: State): CreateNameGroup[] => {
    return localIds.map(localId => state.names.createNameGroups[localId])
}

export const selectCreateNameGroupIdsWithPrimaryRule = createSelector(
    [(state: State) => state.names],
    (names: NamesState) => {
        const result = defaultDict<string, string[]>(_ => [])

        Object.values(names.createNameGroups).forEach(createNameGroup => {
            const createName: CreateName = names.createNames.byLocalId[createNameGroup.primaryCreateNameId]
            result[createName.ruleName].push(createNameGroup.localId)
        })

        return result
    }
)

export const selectCreateNamesByLocalId = (state: State): Record<string, CreateName> => state.names.createNames.byLocalId

export const selectCreateNameWithLocalId = (localId: string) => (state: State): CreateName => {
    return state.names.createNames.byLocalId[localId]
}

export const selectCreateNamesWithLocalIds = (localIds: string[]) => (state: State): CreateName[] => {
    return localIds.map(localId => state.names.createNames.byLocalId[localId])
}

export const selectCreateNameLocalIdsWithGroupId = (groupId: string) => (state: State): string[] => state.names.createNames.byCreateNameGroupLocalId[groupId]

export const selectCounterValuesWithRuleName = createSelector(
    [(state: State) => state.names],
    (names: NamesState): Record<string, Record<string, number>> => {
        const counterValuesByRuleName = defaultDict<string, Record<string, number>>(
            () => defaultDict<string, number>(() => 0)
        )

        Object.values(names.createNames.byLocalId).forEach((createName) => {
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

export const selectCreateNameTokenValuesByLocalId = (state: State): Record<string, CreateNameTokenValue> => state.names.createNameTokenValues.byLocalId

export const selectCreateNameTokenValueLocalIdsWithCreateNameLocalId = (createNameLocalId: string) => (state: State): string[] => {
    return state.names.createNameTokenValues.byCreateNameLocalId[createNameLocalId]
}