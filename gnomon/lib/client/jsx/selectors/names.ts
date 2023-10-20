import { createSelector } from 'reselect'
import * as _ from "lodash"

import { CreateName, CreateNameCompleteCreateName, CreateNameGroup, CreateNameTokenValue } from "../models";
import { CompleteCreateNameParentsByRenderedValues, NamesState } from "../reducers/names"
import { defaultDict } from "../utils/object"
import { State } from '../store';



export const selectNamesState = (state: State): NamesState => state.names

export const selectCreateNameGroupsByLocalId = (state: State): Record<string, CreateNameGroup> => state.names.createNameGroups.byLocalId

export const selectCreateNameGroupWithLocalId = (localId: string) => (state: State): CreateNameGroup => {
    return state.names.createNameGroups.byLocalId[localId]
}

export const selectCreateNameGroupsWithLocalIds = (localIds: string[]) => (state: State): CreateNameGroup[] => {
    return localIds.map(localId => state.names.createNameGroups.byLocalId[localId])
}

export const selectCreateNameGroupIdsByPrimaryRule: (state: State, omitSearchAndReplaceIds?: boolean, respectFilter?: boolean) => Record<string, string[]> = createSelector(
    [
        (state: State) => state.names,
        (state, omitSearchAndReplaceIds = true) => omitSearchAndReplaceIds,
        (state, omitSearchAndReplaceIds, respectFilter = true) => respectFilter,
    ],
    (names: NamesState, omitSearchAndReplaceIds: boolean, respectFilter: boolean): Record<string, string[]> => {

        const result = defaultDict<string, string[]>(_ => [])

        Object.values(names.createNameGroups.byLocalId).forEach(createNameGroup => {
            if (
                omitSearchAndReplaceIds
                && (
                    names.createNameGroups.searchLocalIds.has(createNameGroup.localId)
                    || names.createNameGroups.replaceLocalIds.has(createNameGroup.localId)
                )
            ) {
                return
            }
            // only ignore non-filter IDs if there are any in filter
            if (
                respectFilter
                && names.createNameGroups.filterEnabled === true
                && !names.createNameGroups.filterLocalIds.has(createNameGroup.localId)) {
                return
            }
            const createName: CreateName = names.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId]
            result[createName.ruleName].push(createNameGroup.localId)
        })

        return result
    }
)

export const selectFilterStatus = (state: State): boolean => state.names.createNameGroups.filterEnabled

export const selectCreateNamesByLocalId = (state: State): Record<string, CreateName> => state.names.createNames.byLocalId

export const selectCreateNameWithLocalId = (localId: string) => (state: State): CreateName => {
    return state.names.createNames.byLocalId[localId]
}

export const selectCreateNamesWithLocalIds = (localIds: string[]) => (state: State): CreateName[] => {
    return localIds.map(localId => state.names.createNames.byLocalId[localId])
}

export const selectCreateNameLocalIdsByGroupId = (state: State): Record<string, string[]> => state.names.createNames.byCreateNameGroupLocalId

export const selectCreateNameLocalIdsWithGroupId = (groupId: string) => (state: State): string[] => state.names.createNames.byCreateNameGroupLocalId[groupId]

export const selectSelectedCreateNameGroupIds = (state: State): Set<string> => state.names.createNameGroups.selectionLocalIds

export const selectSearchCreateNameGroupIds = (state: State): Set<string> => state.names.createNameGroups.searchLocalIds

export const selectCreateNameTokenValuesByLocalId = (state: State): Record<string, CreateNameTokenValue> => state.names.createNameTokenValues.byLocalId

export const selectCreateNameTokenValueLocalIdsByCreateNameLocalId = (state: State): Record<string, string[]> => {
    return state.names.createNameTokenValues.byCreateNameLocalId
}

export const selectCreateNameTokenValueLocalIdsWithCreateNameLocalId = (createNameLocalId: string) => (state: State): string[] => {
    return state.names.createNameTokenValues.byCreateNameLocalId[createNameLocalId]
}

export const selectCompleteCreateNamesByParentAndValues = (state: State): CompleteCreateNameParentsByRenderedValues => {
    return state.names.completeCreateNames.byParentLocalIdbyRenderedTokensByCounterValue
}

export const selectRenderedCompleteCreateNameWithLocalId: (state: State, completeCreateNameLocalId: string) => string = createSelector(
    [
        (state: State) => state.names,
        (state, completeCreateNameLocalId) => completeCreateNameLocalId,
    ],
    (names: NamesState, completeCreateNameLocalId: string): string => {

        let renderedValues: string[] = []
        const localIdsToProcess: string[] = [completeCreateNameLocalId]

        while (localIdsToProcess.length) {
            const completeCnLocalId = localIdsToProcess.pop()
            if (!completeCnLocalId) {
                break
            }

            const completeCn = names.completeCreateNames.byLocalId[completeCnLocalId]

            renderedValues.unshift(completeCn.value + completeCn.counterValue || "")

            // grab parents to process
            // maybe recursive call to allow caching?
            if (completeCn.localId in names.completeCreateNameParents.byChildLocalId) {

                const parentCompleteCnLocalIds = [...names.completeCreateNameParents.byChildLocalId[completeCn.localId]]
                    .map(parentLocalId => {
                        return names.completeCreateNameParents.byLocalId[parentLocalId].parentCompleteCreateNameLocalId
                    })

                // TODO: properly handle multiple parents
                // (will be ok with one for now)
                localIdsToProcess.unshift(...parentCompleteCnLocalIds)
            }
        }

        return renderedValues.reduce((prev, curr) => prev + curr, "")
    }
)

export const selectRenderedCompleteCreateNamesByLocalId: (state: State) => Record<string, string> = createSelector(
    [
        (state: State) => state.names,
    ],
    (names: NamesState): Record<string, string> => {

        const renderedNamesByLocalId: Record<string, string> = {}

        // TODO:
        // go through all CCNs
        // use selectRenderedCompleteCreateNameWithLocalId as cache

        return renderedNamesByLocalId
    }
)

export const selectCreateNameCompleteCreateNameLocalIdsByCreateNameLocalId = (state: State): Record<string, string> => {
    return state.names.createNameCompleteCreateNames.byCreateNameLocalId
}

export const selectCreateNameCompleteCreateNamesByLocalId = (state: State): Record<string, CreateNameCompleteCreateName> => {
    return state.names.createNameCompleteCreateNames.byLocalId
}