import { createSelector } from 'reselect'
import * as _ from "lodash"

import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameCompleteCreateName, CreateNameGroup, CreateNameTokenValue, UNSET_VALUE } from "../models";
import { CompleteCreateNameParentsByRenderedValues, NamesState } from "../reducers/names"
import { defaultDict } from "../utils/object"
import { State } from '../store';
import { selectRuleNamesHierarchicalListByPrimaryRuleName } from './rules';



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

export const selectCompleteCreateNameParentLocalIdsByRenderedValues = (state: State): CompleteCreateNameParentsByRenderedValues => {
    return state.names.completeCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
}

export const selectRenderedCompleteCreateNamesByLocalId: (state: State) => Record<string, string> = createSelector(
    [
        (state: State) => state,
    ],
    (state: State): Record<string, string> => {

        const renderedNamesByLocalId: Record<string, string> = {}
        const namesState = state.names

        for (const ccn of Object.values(namesState.completeCreateNames.byLocalId)) {

            let renderedValues: string[] = []
            const localIdsToProcess: string[] = [ccn.localId]

            while (localIdsToProcess.length) {
                const completeCnLocalId = localIdsToProcess.pop()
                if (!completeCnLocalId) {
                    break
                }

                const completeCn = namesState.completeCreateNames.byLocalId[completeCnLocalId]

                renderedValues.unshift(
                    completeCn.value
                    + (completeCn.counterValue != undefined ? String(completeCn.counterValue) : "")
                )

                // grab parents to process
                if (completeCn.localId in namesState.completeCreateNameParents.byChildLocalId) {

                    for (const ccnpLocalId of [...namesState.completeCreateNameParents.byChildLocalId[completeCn.localId]]) {
                        const parentCcnLocalId = namesState.completeCreateNameParents.byLocalId[ccnpLocalId].parentCompleteCreateNameLocalId

                        if (parentCcnLocalId == UNSET_VALUE) {
                            continue
                            // grab precomputed value if exists
                        } else if (parentCcnLocalId in renderedNamesByLocalId) {
                            renderedValues.unshift(renderedNamesByLocalId[parentCcnLocalId])
                            // queue to compute
                        } else {
                            localIdsToProcess.unshift(parentCcnLocalId)
                        }
                    }
                }
            }

            const renderedValue = renderedValues.reduce((prev, curr) => prev + curr, "")
            renderedNamesByLocalId[ccn.localId] = renderedValue
        }

        return renderedNamesByLocalId
    }
)

export const selectRenderedCompleteCreateNamesByCreateNameGroupLocalId: (state: State) => Record<string, string> = createSelector(
    [
        (state: State) => state
    ],
    (state: State): Record<string, string> => {

        const renderedNamesByGroupLocalId: Record<string, string> = {}
        const renderedCompleteCreateNamesByLocalId = selectRenderedCompleteCreateNamesByLocalId(state)

        for (const createNameGroup of Object.values(state.names.createNameGroups.byLocalId)) {

            const targetCompleteCreateNameLocalId: string | undefined = state.names.createNameCompleteCreateNames
                .byCreateNameLocalId[createNameGroup.primaryCreateNameLocalId]

            if (targetCompleteCreateNameLocalId != undefined) {
                renderedNamesByGroupLocalId[createNameGroup.localId] = renderedCompleteCreateNamesByLocalId[targetCompleteCreateNameLocalId]
            }
        }

        return renderedNamesByGroupLocalId
    }
)

// TODO: handle multiple parents
export const selectSortedCompleteCreateNamesWithCreateNameGroupLocalId = (createNameGroupLocalId: string) => (
    state: State
): (CompleteCreateName | undefined)[] => {

    const ruleNamesHierarchicalListByPrimaryRuleName = selectRuleNamesHierarchicalListByPrimaryRuleName(state)
    const createNameGroup = state.names.createNameGroups.byLocalId[createNameGroupLocalId]
    const createNames = state.names.createNames.byCreateNameGroupLocalId[createNameGroupLocalId]
        .map(createNameLocalId => state.names.createNames.byLocalId[createNameLocalId])
    const primaryCreateName = state.names.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId]

    return ruleNamesHierarchicalListByPrimaryRuleName[primaryCreateName.ruleName].map(ruleName => {

        const createName = createNames.find(cn => cn.ruleName == ruleName)
        if (createName == undefined) {
            throw new Error(`can't find CreateName.ruleName equal to ${ruleName}`)
        }

        const completeCreateNameAssociationLocalId = state.names.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId]
        const completeCreateNameAssociation: CreateNameCompleteCreateName | undefined = state.names.createNameCompleteCreateNames
            .byLocalId[completeCreateNameAssociationLocalId]

        if (completeCreateNameAssociation == undefined) {
            return undefined
        }
        return state.names.completeCreateNames.byLocalId[completeCreateNameAssociation.completeCreateNameLocalId]
    })
}

export const selectCreateNameCompleteCreateNameLocalIdsByCreateNameLocalId = (state: State): Record<string, string> => {
    return state.names.createNameCompleteCreateNames.byCreateNameLocalId
}

export const selectCreateNameCompleteCreateNamesByLocalId = (state: State): Record<string, CreateNameCompleteCreateName> => {
    return state.names.createNameCompleteCreateNames.byLocalId
}

export const selectCompleteCreateNameParentLocalIdsByChildLocalId = (state: State): Record<string, string[]> => {
    return state.names.completeCreateNameParents.byChildLocalId
}

export const selectCompleteCreateNameParentsByLocalId = (state: State): Record<string, CompleteCreateNameParent> => {
    return state.names.completeCreateNameParents.byLocalId
}