import * as _ from "lodash"

import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameCompleteCreateName, CreateNameGroup, CreateNameTokenValue, RuleParent, RuleToken, TokenValue, UNSET_VALUE } from "../models"
import { makeActionObject } from "./utils"
import { createLocalId } from "../utils/models"
import { defaultDict } from "../utils/object"
import { NamesState } from "../reducers/names"
import { createSearchReplaceCriteriaFromGroups, renderCounter, renderTokens } from "../utils/names"
import { State } from "../store"
import { RulesStateSliceForCompleteCreateNames } from "../selectors/global"



export const ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES = "ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES"
export const ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES = "ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES"
export const DELETE_CREATE_NAME_TOKEN_VALUE = "DELETE_CREATE_NAME_TOKEN_VALUE"
export const SET_COUNTER_VALUE_FOR_CREATE_NAMES = "SET_COUNTER_VALUE_FOR_CREATE_NAMES"
export const ADD_CREATE_NAME_GROUPS_TO_SELECTION = "ADD_CREATE_NAME_GROUPS_TO_SELECTION"
export const REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION = "REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION"
export const CLEAR_CREATE_NAME_GROUPS_SELECTION = "CLEAR_CREATE_NAME_GROUPS_SELECTION"
export const CLEAR_CREATE_NAME_GROUPS_FILTER = "CLEAR_CREATE_NAME_GROUPS_FILTER"
export const DELETE_GROUPS_WITH_NAMES = "DELETE_GROUPS_WITH_NAMES"
export const DELETE_SELECTED_GROUPS_WITH_NAMES = "DELETE_SELECTED_GROUPS_WITH_NAMES"
export const ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA = "ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA"
export const REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA = "REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA"
export const SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA = "SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA"
export const SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA = "SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA"
export const ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA = "ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA"
export const REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA = "REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA"
export const ADD_OR_REPLACE_COMPLETE_CREATE_NAMES = "ADD_OR_REPLACE_COMPLETE_CREATE_NAMES"
export const REMOVE_COMPLETE_CREATE_NAMES = "REMOVE_COMPLETE_CREATE_NAMES"


export function addNamesWithGroupsWithTokenValues(createNames: CreateName[], createNameGroups: CreateNameGroup[], createNameTokenValues: CreateNameTokenValue[]) {
    return makeActionObject(ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES, { createNames, createNameGroups, createNameTokenValues })
}

export function createNamesWithGroupForRule(
    primaryRuleName: string,
    state: State,
) {

    const rulesState = state.rules
    const createNames: CreateName[] = []
    const ruleNamesToScan: string[] = [primaryRuleName]

    const createNameGroup: CreateNameGroup = {
        localId: createLocalId(),
        primaryCreateNameLocalId: "TODO",
    }

    const createNameTokenValues: CreateNameTokenValue[] = []

    // make all CreateNames
    while (ruleNamesToScan.length) {
        const ruleName = ruleNamesToScan.pop()
        // typescript needs this check even though
        // we're already checking ruleNamesToScan.length
        if (!ruleName) {
            break
        }
        const createName: CreateName = {
            localId: createLocalId(),
            ruleName,
            createNameGroupLocalId: createNameGroup.localId
        }

        createNames.push(createName)
        if (ruleName == primaryRuleName) {
            createNameGroup.primaryCreateNameLocalId = createName.localId
        }

        // make CreateNameTokenValue if only one value for a given token
        const ruleTokens = rulesState.ruleTokens.byRuleName[ruleName].map(rtLocalId => rulesState.ruleTokens.byLocalId[rtLocalId])
        ruleTokens.forEach(ruleToken => {
            const tokenValueLocalIds = rulesState.tokenValues.byTokenName[ruleToken.tokenName]

            if (tokenValueLocalIds.length == 1) {
                createNameTokenValues.push({
                    localId: createLocalId(),
                    tokenValueLocalId: tokenValueLocalIds[0],
                    createNameLocalId: createName.localId,
                    ruleTokenLocalId: ruleToken.localId,
                })
            }
        })

        if (rulesState.ruleParents.byRuleName[ruleName]) {
            const parentRuleNames = rulesState.ruleParents.byRuleName[ruleName].map(rpLocalId => {
                return rulesState.ruleParents.byLocalId[rpLocalId].parentRuleName
            })
            // TODO: properly handle multiple parents
            // (will be ok with one for now)
            ruleNamesToScan.push(...parentRuleNames)
        }
    }

    return addNamesWithGroupsWithTokenValues(createNames, [createNameGroup], createNameTokenValues)
}


interface DuplicateCreateNameGroupReturn {
    createNames: CreateName[]
    createNameGroup: CreateNameGroup
    createNameTokenValues: CreateNameTokenValue[]
}


function _duplicateCreateNameGroup(
    createNameGroup: CreateNameGroup,
    namesState: NamesState,
    ruleCounterValuesByRuleName: Record<string, number> = {},
): DuplicateCreateNameGroupReturn {

    const newCng: CreateNameGroup = {
        localId: createLocalId(),
        primaryCreateNameLocalId: "TODO",
    }
    const newCntvs: CreateNameTokenValue[] = []

    const newCns: CreateName[] = namesState.createNames.byCreateNameGroupLocalId[createNameGroup.localId].map(oldCnLocalId => {
        const oldCn = namesState.createNames.byLocalId[oldCnLocalId]

        const newCn: CreateName = {
            ...oldCn,
            localId: createLocalId(),
            createNameGroupLocalId: newCng.localId,
        }

        if (newCn.ruleName in ruleCounterValuesByRuleName) {
            newCn.ruleCounterValue = ruleCounterValuesByRuleName[newCn.ruleName]
        }

        if (createNameGroup.primaryCreateNameLocalId == oldCn.localId) {
            newCng.primaryCreateNameLocalId = newCn.localId
        }

        const _newCnTvs: CreateNameTokenValue[] = (namesState.createNameTokenValues.byCreateNameLocalId[oldCn.localId] || []).map(oldCntvLocalId => {
            const oldCntv = namesState.createNameTokenValues.byLocalId[oldCntvLocalId]

            return {
                ...oldCntv,
                localId: createLocalId(),
                createNameLocalId: newCn.localId,
            }
        })
        newCntvs.push(..._newCnTvs)

        return newCn
    })

    return { createNames: newCns, createNameGroup: newCng, createNameTokenValues: newCntvs }
}

// maybe just pass entire state?
export function duplicateCreateNameGroups(
    createNameGroups: CreateNameGroup[],
    state: State,
    quantity: number = 1,
    overrideRuleCounterValuesByRuleName: Record<string, number> = {},
) {

    const newCreateNames: CreateName[] = []
    const newCreateNameGroups: CreateNameGroup[] = []
    const newCreateNameTokenValues: CreateNameTokenValue[] = []

    createNameGroups.forEach(cng => {

        for (let index = 1; index <= quantity; index++) {

            const {
                createNames: _newCreateNames,
                createNameGroup: _newCreateNameGroup,
                createNameTokenValues: _newCreateNameTokenValues
            } = _duplicateCreateNameGroup(
                cng,
                state.names,
                overrideRuleCounterValuesByRuleName,
            )

            newCreateNames.push(..._newCreateNames)
            newCreateNameGroups.push(_newCreateNameGroup)
            newCreateNameTokenValues.push(..._newCreateNameTokenValues)
        }
    })

    return addNamesWithGroupsWithTokenValues(newCreateNames, newCreateNameGroups, newCreateNameTokenValues)
}

export function iterateOnCreateNameGroupsByRule(
    createNameGroups: CreateNameGroup[],
    ruleName: string,
    startIdx: number,
    endIdx: number,
    state: State,
) {

    const createNameGroupsByPrimaryRuleName = defaultDict<string, CreateNameGroup[]>(_ => [])

    createNameGroups.forEach(cng => {
        const primaryRuleName = state.names.createNames.byLocalId[cng.primaryCreateNameLocalId].ruleName

        createNameGroupsByPrimaryRuleName[primaryRuleName].push(cng)
    })

    const newCreateNames: CreateName[] = []
    const newCreateNameGroups: CreateNameGroup[] = []
    const newCreateNameTokenValues: CreateNameTokenValue[] = []

    Object.entries(createNameGroupsByPrimaryRuleName).forEach(([_, createNameGroups]) => {

        for (let index = startIdx; index <= endIdx; index++) {

            createNameGroups.forEach(cng => {

                const overrideRuleCounterValuesByRuleName = { [ruleName]: index }

                const {
                    createNames: _newCreateNames,
                    createNameGroup: _newCreateNameGroup,
                    createNameTokenValues: _newCreateNameTokenValues
                } = _duplicateCreateNameGroup(
                    cng,
                    state.names,
                    overrideRuleCounterValuesByRuleName,
                )

                newCreateNames.push(..._newCreateNames)
                newCreateNameGroups.push(_newCreateNameGroup)
                newCreateNameTokenValues.push(..._newCreateNameTokenValues)
            })
        }
    })

    return addNamesWithGroupsWithTokenValues(newCreateNames, newCreateNameGroups, newCreateNameTokenValues)
}

export function addOrReplaceCreateNameTokenValues(newCreateNameTokenValues: CreateNameTokenValue[], oldCreateNameTokenValueLocalIds: string[]) {
    return makeActionObject(ADD_AND_REMOVE_CREATE_NAME_TOKEN_VALUES, { newCreateNameTokenValues, oldCreateNameTokenValueLocalIds })
}

export function deleteCreateNameTokenValue(createNameTokenValue: CreateNameTokenValue) {
    return makeActionObject(DELETE_CREATE_NAME_TOKEN_VALUE, { createNameTokenValue })
}

export function setCreateNameRuleCounterValues(ruleCounterValuesByCreateNameLocalId: Record<string, number | undefined>) {
    return makeActionObject(SET_COUNTER_VALUE_FOR_CREATE_NAMES, { ruleCounterValuesByCreateNameLocalId })
}

export function addCreateNameGroupsToSelection(createNameGroupIds: string[]) {
    return makeActionObject(ADD_CREATE_NAME_GROUPS_TO_SELECTION, { createNameGroupIds })
}

export function removeCreateNameGroupsFromSelection(createNameGroupIds: string[]) {
    return makeActionObject(REMOVE_CREATE_NAME_GROUPS_FROM_SELECTION, { createNameGroupIds })
}

export function clearCreateNameGroupsSelection() {
    return makeActionObject(CLEAR_CREATE_NAME_GROUPS_SELECTION, {})
}

export function setCreateNameGroupsSelectionFromSearchCriteria() {
    return makeActionObject(SET_CREATE_NAME_GROUPS_SELECTION_FROM_SEARCH_CRITERIA, {})
}

export function setCreateNameGroupsFilterFromSearchCriteria() {
    return makeActionObject(SET_CREATE_NAME_GROUPS_FILTER_FROM_SEARCH_CRITERIA, {})
}

export function clearCreateNameGroupsFilter() {
    return makeActionObject(CLEAR_CREATE_NAME_GROUPS_FILTER, {})
}

export function deleteGroupsWithNames(createNameGroupIds: string[]) {
    return makeActionObject(DELETE_GROUPS_WITH_NAMES, { createNameGroupIds })
}

export function deleteSelectedGroupsWithNames() {
    return makeActionObject(DELETE_SELECTED_GROUPS_WITH_NAMES, {})
}

export function addCreateNameGroupsToSearchCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(ADD_CREATE_NAME_GROUPS_TO_SEARCH_CRITERIA, { createNameGroupLocalIds })
}

export function removeCreateNameGroupsFromSearchCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(REMOVE_CREATE_NAME_GROUPS_FROM_SEARCH_CRITERIA, { createNameGroupLocalIds })
}

export function addCreateNameGroupsToReplaceCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(ADD_CREATE_NAME_GROUPS_TO_REPLACE_CRITERIA, { createNameGroupLocalIds })
}

export function removeCreateNameGroupsFromReplaceCriteria(createNameGroupLocalIds: string[]) {
    return makeActionObject(REMOVE_CREATE_NAME_GROUPS_FROM_REPLACE_CRITERIA, { createNameGroupLocalIds })
}

export function addOrReplaceCreateNameTokenValuesAndRuleCounterValuesFromReplaceCriteria(state: State) {
    const namesState = state.names

    // TODO: support multiple CreateNameGroups for criteria?
    const replaceCriteria = createSearchReplaceCriteriaFromGroups(
        namesState, namesState.createNameGroups.replaceLocalIds
    )[0].byRuleName

    const newCntvs: CreateNameTokenValue[] = []
    const oldCntvLocalIds: string[] = []
    const newCounterValuesByCnLocalId: Record<string, number | undefined> = {}
    const modifiedCreateNameGroupIds: Set<string> = new Set()

    for (const selectionGroupId of namesState.createNameGroups.selectionLocalIds) {
        const selectionCnLocalIds = namesState.createNames.byCreateNameGroupLocalId[selectionGroupId]

        for (const selectionCnLocalId of selectionCnLocalIds) {
            const selectionCn = namesState.createNames.byLocalId[selectionCnLocalId]

            const selectionCntvs = (namesState.createNameTokenValues.byCreateNameLocalId[selectionCnLocalId] || [])
                .map(cntvLocalId => namesState.createNameTokenValues.byLocalId[cntvLocalId])
            const replaceCntvs = (replaceCriteria[selectionCn.ruleName]?.createNameTokenValueLocalIds || [])
                .map(cntvLocalId => namesState.createNameTokenValues.byLocalId[cntvLocalId])

            for (const replaceCntv of replaceCntvs) {

                // create or replace TokenValues
                const selectionCntv: CreateNameTokenValue | undefined = selectionCntvs.find(cntv => {
                    return cntv.ruleTokenLocalId == replaceCntv.ruleTokenLocalId
                })

                if (selectionCntv?.tokenValueLocalId == replaceCntv.tokenValueLocalId) {
                    continue
                }

                newCntvs.push({
                    localId: createLocalId(),
                    tokenValueLocalId: replaceCntv.tokenValueLocalId,
                    createNameLocalId: selectionCn.localId,
                    ruleTokenLocalId: replaceCntv.ruleTokenLocalId,
                })

                if (selectionCntv) {
                    oldCntvLocalIds.push(selectionCntv.localId)
                }

                modifiedCreateNameGroupIds.add(selectionGroupId)
            }

            // update counterValue
            const newCounterValue = replaceCriteria[selectionCn.ruleName]?.ruleCounterValue

            if (newCounterValue != undefined && selectionCn.ruleCounterValue != newCounterValue) {
                newCounterValuesByCnLocalId[selectionCn.localId] = newCounterValue

                modifiedCreateNameGroupIds.add(selectionGroupId)
            }
        }
    }

    return [
        addOrReplaceCreateNameTokenValues(newCntvs, oldCntvLocalIds),
        setCreateNameRuleCounterValues(newCounterValuesByCnLocalId),
    ]
}

// TODO: handle multiparent rules
export function addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
    createNameGroupsLocalIds: string[],
    rulesState: RulesStateSliceForCompleteCreateNames,
    namesState: NamesState,
) {
    const completeCreateNamesToAdd: CompleteCreateName[] = []
    const completeCreateNamesToRemove: CompleteCreateName[] = []
    const createNameCompleteCreateNamesToAdd: CreateNameCompleteCreateName[] = []
    const createNameCompleteCreateNamesToRemove: CreateNameCompleteCreateName[] = []
    const completeCreateNameParentsToAdd: CompleteCreateNameParent[] = []
    const completeCreateNameParentsToRemove: CompleteCreateNameParent[] = []

    // for each CreateNameGroup
    for (const cngLocalId of createNameGroupsLocalIds) {

        const createNameGroup = namesState.createNameGroups.byLocalId[cngLocalId]
        const primaryRuleName = namesState.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId].ruleName
        const ruleNamesHierarchicalList = rulesState.ruleNamesHierarchicalListByPrimaryRuleName[primaryRuleName]

        const createNames = namesState.createNames.byCreateNameGroupLocalId[cngLocalId]
            .map(cnLocalId => namesState.createNames.byLocalId[cnLocalId])

        const sortedCreateNames: CreateName[] = _.sortBy(createNames, [(cn: CreateName) => {
            return ruleNamesHierarchicalList.indexOf(cn.ruleName)
        }])

        // keep track of hierarchy completeness
        const completeCreateNames: (CompleteCreateName | undefined)[] = []
        let completeHierarchy = true

        // for each CreateName
        for (const [idx, createName] of sortedCreateNames.entries()) {
            const createNameTokenValues = (namesState.createNameTokenValues.byCreateNameLocalId[createName.localId] || [])
                .map(cntvLocalId => namesState.createNameTokenValues.byLocalId[cntvLocalId])

            const actualTokenCount = createNameTokenValues.length
            const targetTokenCount = rulesState.ruleTokenLocalIdsByRuleName[createName.ruleName].length
            const counterRequired = rulesState.counterRequiredByRuleName[createName.ruleName]
            const hasCounterValue = createName.ruleCounterValue != undefined

            let parentCompleteCreateName: CompleteCreateName | undefined = undefined

            if (idx > 0) {
                parentCompleteCreateName = completeCreateNames[idx - 1]

                if (parentCompleteCreateName == undefined) {
                    completeHierarchy = false
                }
            }

            // detect completeness
            if (
                targetTokenCount === actualTokenCount
                && (
                    !counterRequired
                    || hasCounterValue
                )
                && completeHierarchy
            ) {
                // is complete

                let completeCreateNameParentLocalId: string | undefined = undefined
                const hierarchyRenderedTokenValue = renderTokens(
                    createNameTokenValues,
                    rulesState.ruleTokensByLocalId,
                    rulesState.tokenValuesByLocalId,
                )
                const hierarchyRenderedCounterValue = createName.ruleCounterValue != undefined ? renderCounter(createName) : UNSET_VALUE

                // find completeCreateNameParentLocalId
                if (parentCompleteCreateName != undefined) {
                    const renderedTokenValuesDict = namesState.completeCreateNameParents
                        .byParentLocalIdByChildRenderedTokensByChildCounterValue[parentCompleteCreateName.localId]

                    if (renderedTokenValuesDict != undefined) {
                        const renderedCounterValuesDict = renderedTokenValuesDict[hierarchyRenderedTokenValue]

                        if (renderedCounterValuesDict != undefined) {
                            completeCreateNameParentLocalId = renderedCounterValuesDict[hierarchyRenderedCounterValue]
                        }
                    }
                }

                let completeCreateNameParent: CompleteCreateNameParent | undefined = undefined

                if (completeCreateNameParentLocalId) {
                    completeCreateNameParent = namesState.completeCreateNameParents.byLocalId[completeCreateNameParentLocalId]
                }


                // create CompleteCreateName and CreateName association if doesn't exist
                let completeCreateName: CompleteCreateName | undefined = undefined

                if (completeCreateNameParent == undefined) {
                    completeCreateName = {
                        localId: createLocalId(),
                        value: renderTokens(
                            createNameTokenValues,
                            rulesState.ruleTokensByLocalId,
                            rulesState.tokenValuesByLocalId,
                        ),
                        counterValue: createName.ruleCounterValue,
                    }

                    completeCreateNamesToAdd.push(completeCreateName)

                    const completeCreateNameAssociation: CreateNameCompleteCreateName = {
                        localId: createLocalId(),
                        createNameLocalId: createName.localId,
                        completeCreateNameLocalId: completeCreateName.localId,
                    }

                    createNameCompleteCreateNamesToAdd.push(completeCreateNameAssociation)
                } else {
                    completeCreateName = namesState.completeCreateNames.byLocalId[completeCreateNameParent.completeCreateNameLocalId]
                }


                // create CompleteCreateNameParent if doesn't exist
                if (
                    completeCreateNameParent == undefined
                    && parentCompleteCreateName != undefined
                ) {
                    const parent: CompleteCreateNameParent = {
                        localId: createLocalId(),
                        completeCreateNameLocalId: completeCreateName.localId,
                        parentCompleteCreateNameLocalId: parentCompleteCreateName.localId,
                    }

                    completeCreateNameParentsToAdd.push(parent)
                }

                completeCreateNames.push(completeCreateName)

            } else {
                // is not complete
                // remove all records if exist

                const completeCreateNameAssociationLocalId = namesState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId]

                if (completeCreateNameAssociationLocalId != undefined) {
                    const completeCreateNameAssociation = namesState.createNameCompleteCreateNames.byLocalId[completeCreateNameAssociationLocalId]
                    const completeCreateName = namesState.completeCreateNames.byLocalId[completeCreateNameAssociation.completeCreateNameLocalId]
                    const completeCreateNameParents = [
                        ...namesState.completeCreateNameParents.byParentLocalId[completeCreateName.localId],
                        ...namesState.completeCreateNameParents.byChildLocalId[completeCreateName.localId],
                    ]
                        .map(ccnpLocalId => namesState.completeCreateNameParents.byLocalId[ccnpLocalId])

                    createNameCompleteCreateNamesToRemove.push(completeCreateNameAssociation)
                    completeCreateNamesToRemove.push(completeCreateName)
                    completeCreateNameParentsToRemove.push(...completeCreateNameParents)
                }

                completeCreateNames.push(undefined)
            }
        }
    }

    return makeActionObject(
        ADD_OR_REPLACE_COMPLETE_CREATE_NAMES,
        {
            completeCreateNamesToAdd,
            completeCreateNameParentsToAdd,
            completeCreateNameLocalIdsToRemove: completeCreateNamesToRemove.map(x => x.localId),
            completeCreateNameParentLocalIdsToRemove: completeCreateNameParentsToRemove.map(x => x.localId),
        }
    )
}

export function removeCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
    createNameGroupsLocalIds: string[],
    namesState: NamesState,
) {
    const completeCreateNamesToRemove: CompleteCreateName[] = []
    const createNameCompleteCreateNamesToRemove: CreateNameCompleteCreateName[] = []
    const completeCreateNameParentsToRemove: CompleteCreateNameParent[] = []

    // for each CreateNameGroup
    for (const cngLocalId of createNameGroupsLocalIds) {

        const createNames = namesState.createNames.byCreateNameGroupLocalId[cngLocalId]
            .map(cnLocalId => namesState.createNames.byLocalId[cnLocalId])

        for (const createName of createNames) {

            const completeCreateNameAssociationLocalId = namesState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId]

            if (completeCreateNameAssociationLocalId != undefined) {
                const completeCreateNameAssociation = namesState.createNameCompleteCreateNames.byLocalId[completeCreateNameAssociationLocalId]
                const completeCreateName = namesState.completeCreateNames.byLocalId[completeCreateNameAssociation.completeCreateNameLocalId]
                const completeCreateNameParents = [
                    ...namesState.completeCreateNameParents.byParentLocalId[completeCreateName.localId],
                    ...namesState.completeCreateNameParents.byChildLocalId[completeCreateName.localId],
                ]
                    .map(ccnpLocalId => namesState.completeCreateNameParents.byLocalId[ccnpLocalId])

                createNameCompleteCreateNamesToRemove.push(completeCreateNameAssociation)
                completeCreateNamesToRemove.push(completeCreateName)
                completeCreateNameParentsToRemove.push(...completeCreateNameParents)
            }
        }
    }


    return makeActionObject(
        REMOVE_COMPLETE_CREATE_NAMES,
        {
            completeCreateNameLocalIdsToRemove: completeCreateNamesToRemove.map(x => x.localId),
            createNameCompleteCreateNameLocalIdsToRemove: completeCreateNamesToRemove.map(x => x.localId),
            completeCreateNameParentLocalIdsToRemove: completeCreateNameParentsToRemove.map(x => x.localId),
        }
    )
}

export type ACTION_TYPE =
    | ReturnType<typeof addNamesWithGroupsWithTokenValues>
    | ReturnType<typeof duplicateCreateNameGroups>
    | ReturnType<typeof addOrReplaceCreateNameTokenValues>
    | ReturnType<typeof deleteCreateNameTokenValue>
    | ReturnType<typeof setCreateNameRuleCounterValues>
    | ReturnType<typeof addCreateNameGroupsToSelection>
    | ReturnType<typeof removeCreateNameGroupsFromSelection>
    | ReturnType<typeof clearCreateNameGroupsSelection>
    | ReturnType<typeof setCreateNameGroupsSelectionFromSearchCriteria>
    | ReturnType<typeof setCreateNameGroupsFilterFromSearchCriteria>
    | ReturnType<typeof clearCreateNameGroupsFilter>
    | ReturnType<typeof deleteGroupsWithNames>
    | ReturnType<typeof deleteSelectedGroupsWithNames>
    | ReturnType<typeof addCreateNameGroupsToSearchCriteria>
    | ReturnType<typeof removeCreateNameGroupsFromSearchCriteria>
    | ReturnType<typeof addCreateNameGroupsToReplaceCriteria>
    | ReturnType<typeof removeCreateNameGroupsFromReplaceCriteria>
    | ReturnType<typeof addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds>
    | ReturnType<typeof removeCompleteCreateNamesAndParentsForCreateNameGroupLocalIds>
