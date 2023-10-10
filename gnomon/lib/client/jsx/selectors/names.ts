import { createSelector } from 'reselect'
import * as _ from "lodash"

import { CreateName, CreateNameGroup, CreateNameTokenValue } from "../models";
import { NamesState } from "../reducers/names"
import { defaultDict } from "../utils/object"


interface State {
    names: NamesState
}


export const selectCreateNameGroupsByLocalId = (state: State): Record<string, CreateNameGroup> => state.names.createNameGroups

export const selectCreateNameGroupWithLocalId = (localId: string) => (state: State): CreateNameGroup => {
    return state.names.createNameGroups[localId]
}

export const selectCreateNameGroupsWithLocalIds = (localIds: string[]) => (state: State): CreateNameGroup[] => {
    return localIds.map(localId => state.names.createNameGroups[localId])
}

export const selectCreateNameGroupIdsWithPrimaryRule: (state: State, localOnlyValue?: boolean) => Record<string, string[]> = createSelector(
    [
        (state: State) => state.names,
        (state, localOnlyValue) => localOnlyValue,
    ],
    (names: NamesState, localOnlyValue: boolean): Record<string, string[]> => {
        const result = defaultDict<string, string[]>(_ => [])

        Object.values(names.createNameGroups).forEach(createNameGroup => {
            if (localOnlyValue != undefined && createNameGroup.localOnly != localOnlyValue) {
                return
            }
            const createName: CreateName = names.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId]
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

export const selectCreateNameLocalIdsByGroupId = (state: State): Record<string, string[]> => state.names.createNames.byCreateNameGroupLocalId

export const selectCreateNameLocalIdsWithGroupId = (groupId: string) => (state: State): string[] => state.names.createNames.byCreateNameGroupLocalId[groupId]

export const selectRuleCounterValuesbyRuleName: (state: State) => Record<string, Record<string, number>> = createSelector(
    [(state: State) => state.names],
    (names: NamesState): Record<string, Record<string, number>> => {
        const counterValuesByRuleName = defaultDict<string, Record<string, number>>(
            () => defaultDict<string, number>(() => 0)
        )

        Object.values(names.createNames.byLocalId).forEach((createName) => {
            if (createName.ruleCounterValue != undefined) {
                ++counterValuesByRuleName[createName.ruleName][createName.ruleCounterValue]
            }
        })

        return counterValuesByRuleName
    }
)

export const selectSelectedCreateNameGroupIds: (state: State) => string[] = createSelector(
    [(state: State) => state.names],
    (names: NamesState): string[] => {
        return _.filter(names.createNameGroups, { selected: true }).map(cng => cng.localId)
    }
)

export const selectCreateNameTokenValuesByLocalId = (state: State): Record<string, CreateNameTokenValue> => state.names.createNameTokenValues.byLocalId

export const selectCreateNameTokenValueLocalIdsByCreateNameLocalId = (state: State): Record<string, string[]> => {
    return state.names.createNameTokenValues.byCreateNameLocalId
}

export const selectCreateNameTokenValueLocalIdsWithCreateNameLocalId = (createNameLocalId: string) => (state: State): string[] => {
    return state.names.createNameTokenValues.byCreateNameLocalId[createNameLocalId]
}