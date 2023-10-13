import { CreateName, CreateNameGroup, CreateNameTokenValue, RuleParent, RuleToken } from "../models"
import { makeActionObject } from "./utils"
import { createLocalId } from "../utils/models"
import { defaultDict } from "../utils/object"



export const ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES = "ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES"
export const ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE = "ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE"
export const DELETE_CREATE_NAME_TOKEN_VALUE = "DELETE_CREATE_NAME_TOKEN_VALUE"
export const SET_COUNTER_VALUE_FOR_CREATE_NAME = "SET_COUNTER_VALUE_FOR_CREATE_NAME"
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
export const REPLACE_VALUES_FROM_REPLACE_CRITERIA = "REPLACE_VALUES_FROM_REPLACE_CRITERIA"


// must add CreateNameTokenValues per RuleToken if only one TokenValue per Token
export function addNamesWithGroupsWithTokenValues(createNames: CreateName[], createNameGroups: CreateNameGroup[], createNameTokenValues: CreateNameTokenValue[]) {
    return makeActionObject(ADD_CREATE_NAMES_WITH_GROUPS_WITH_TOKEN_VALUES, { createNames, createNameGroups, createNameTokenValues })
}

export function createNamesWithGroupForRule(
    primaryRuleName: string,
    ruleParentLocalIdsByRuleName: Record<string, string[]>,
    ruleParentsByLocalId: Record<string, RuleParent>,
    ruleTokenLocalIdsByRuleName: Record<string, string[]>,
    ruleTokensByLocalId: Record<string, RuleToken>,
    tokenValueLocalIdsByTokenName: Record<string, string[]>,
): ReturnType<typeof addNamesWithGroupsWithTokenValues> {

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
        const ruleTokens = ruleTokenLocalIdsByRuleName[ruleName].map(rtLocalId => ruleTokensByLocalId[rtLocalId])
        ruleTokens.forEach(ruleToken => {
            const tokenValueLocalIds = tokenValueLocalIdsByTokenName[ruleToken.tokenName]

            if (tokenValueLocalIds.length == 1) {
                createNameTokenValues.push({
                    localId: createLocalId(),
                    tokenValueLocalId: tokenValueLocalIds[0],
                    createNameLocalId: createName.localId,
                    ruleTokenLocalId: ruleToken.localId,
                })
            }
        })

        if (ruleParentLocalIdsByRuleName[ruleName]) {
            const parentRuleNames = ruleParentLocalIdsByRuleName[ruleName].map(rpLocalId => {
                return ruleParentsByLocalId[rpLocalId].parentRuleName
            })
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
    createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]>,
    createNamesByLocalId: Record<string, CreateName>,
    createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]>,
    createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue>,
    ruleCounterValuesByRuleName: Record<string, number> = {},
): DuplicateCreateNameGroupReturn {

    const newCng: CreateNameGroup = {
        localId: createLocalId(),
        primaryCreateNameLocalId: "TODO",
    }
    const newCntvs: CreateNameTokenValue[] = []

    const newCns: CreateName[] = createNameLocalIdsByCreateNameGroupLocalId[createNameGroup.localId].map(oldCnLocalId => {
        const oldCn = createNamesByLocalId[oldCnLocalId]

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

        const _newCnTvs: CreateNameTokenValue[] = (createNameTokenValueLocalIdsByCreateNameLocalId[oldCn.localId] || []).map(oldCntvLocalId => {
            const oldCntv = createNameTokenValuesByLocalId[oldCntvLocalId]

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
    createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]>,
    createNamesByLocalId: Record<string, CreateName>,
    createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]>,
    createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue>,
    quantity: number = 1,
    overrideRuleCounterValuesByRuleName: Record<string, number> = {},
): ReturnType<typeof addNamesWithGroupsWithTokenValues> {

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
                createNameLocalIdsByCreateNameGroupLocalId,
                createNamesByLocalId,
                createNameTokenValueLocalIdsByCreateNameLocalId,
                createNameTokenValuesByLocalId,
                overrideRuleCounterValuesByRuleName
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
    createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]>,
    createNamesByLocalId: Record<string, CreateName>,
    createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]>,
    createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue>,
): ReturnType<typeof addNamesWithGroupsWithTokenValues> {

    const createNameGroupsByPrimaryRuleName = defaultDict<string, CreateNameGroup[]>(_ => [])

    createNameGroups.forEach(cng => {
        const primaryRuleName = createNamesByLocalId[cng.primaryCreateNameLocalId].ruleName

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
                    createNameLocalIdsByCreateNameGroupLocalId,
                    createNamesByLocalId,
                    createNameTokenValueLocalIdsByCreateNameLocalId,
                    createNameTokenValuesByLocalId,
                    overrideRuleCounterValuesByRuleName
                )

                newCreateNames.push(..._newCreateNames)
                newCreateNameGroups.push(_newCreateNameGroup)
                newCreateNameTokenValues.push(..._newCreateNameTokenValues)
            })
        }
    })

    return addNamesWithGroupsWithTokenValues(newCreateNames, newCreateNameGroups, newCreateNameTokenValues)
}

export function addOrReplaceCreateNameTokenValue(newCreateNameTokenValue: CreateNameTokenValue, oldCreateNameTokenValue?: CreateNameTokenValue) {
    return makeActionObject(ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE, { newCreateNameTokenValue, oldCreateNameTokenValue })
}

export function deleteCreateNameTokenValue(createNameTokenValue: CreateNameTokenValue) {
    return makeActionObject(DELETE_CREATE_NAME_TOKEN_VALUE, { createNameTokenValue })
}

export function setCreateNameRuleCounterValue(createNameLocalId: string, ruleCounterValue: number | undefined) {
    return makeActionObject(SET_COUNTER_VALUE_FOR_CREATE_NAME, { createNameLocalId, ruleCounterValue })
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

export function replaceValuesFromReplaceCriteria() {
    return makeActionObject(REPLACE_VALUES_FROM_REPLACE_CRITERIA, {})
}


export type ACTION_TYPE =
    | ReturnType<typeof addNamesWithGroupsWithTokenValues>
    | ReturnType<typeof duplicateCreateNameGroups>
    | ReturnType<typeof addOrReplaceCreateNameTokenValue>
    | ReturnType<typeof deleteCreateNameTokenValue>
    | ReturnType<typeof setCreateNameRuleCounterValue>
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
    | ReturnType<typeof replaceValuesFromReplaceCriteria>
