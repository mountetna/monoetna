import { CreateName, CreateNameGroup, CreateNameTokenValue, RuleParent, RuleToken } from "../models"
import { makeActionObject } from "./utils"
import { createLocalId } from "../utils/models"



export const ADD_CREATE_NAMES_WITH_GROUP = "ADD_CREAT_NAMES_WITH_GROUP"
export const DUPLICATE_CREATE_NAME_GROUP = "DUPLICATE_CREATE_NAME_GROUP"
export const ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE = "ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE"
export const DELETE_CREATE_NAME_TOKEN_VALUE = "DELETE_CREATE_NAME_TOKEN_VALUE"
export const SET_COUNTER_VALUE_FOR_CREATE_NAME = "SET_COUNTER_VALUE_FOR_CREATE_NAME"
export const SET_CREATE_NAME_GROUPS_SELECTED = "SET_CREATE_NAME_GROUPS_SELECTED"
export const DELETE_GROUPS_WITH_NAMES = "DELETE_GROUPS_WITH_NAMES"
export const DELETE_SELECTED_GROUPS_WITH_NAMES = "DELETE_SELECTED_GROUPS_WITH_NAMES"


interface AddNamePayload {
    createNames: CreateName[]
    createNameGroup: CreateNameGroup
}


export function createNamesWithGroupForRule(
    primaryRuleName: string,
    ruleParentLocalIdsByRuleName: Record<string, string[]>,
    ruleParentsByLocalId: Record<string, RuleParent>,
    ruleTokenLocalIdsByRuleName: Record<string, string[]>,
    ruleTokensByLocalId: Record<string, RuleToken>,
    tokenValueNamesByTokenName: Record<string, string[]>,
): AddNamePayload {

    const createNames: CreateName[] = []
    const ruleNamesToScan: string[] = [primaryRuleName]

    const createNameGroup: CreateNameGroup = {
        localId: createLocalId(),
        primaryCreateNameId: "TODO",
        selected: false
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
            createNameGroup.primaryCreateNameId = createName.localId
        }

        // make CreateNameTokenValue if only one value for a given token
        const ruleTokens = ruleTokenLocalIdsByRuleName[ruleName].map(rtLocalId => ruleTokensByLocalId[rtLocalId])
        ruleTokens.forEach(ruleToken => {
            const tokenValueNames = tokenValueNamesByTokenName[ruleToken.tokenName]

            if (tokenValueNames.length == 1) {
                createNameTokenValues.push({
                    localId: createLocalId(),
                    tokenValueName: tokenValueNames[0],
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

    return addNamesWithGroup(createNames, createNameGroup, createNameTokenValues)
}

// must add CreateNameTokenValues per RuleToken if only one TokenValue per Token
export function addNamesWithGroup(createNames: CreateName[], createNameGroup: CreateNameGroup, createNameTokenValues: CreateNameTokenValue[]) {
    return makeActionObject(ADD_CREATE_NAMES_WITH_GROUP, { createNames, createNameGroup, createNameTokenValues })
}

export function duplicateCreateNameGroup(
    createNameGroup: CreateNameGroup,
    createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]>,
    createNamesByLocalId: Record<string, CreateName>,
    createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]>,
    createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue>,
) {
    const newCng: CreateNameGroup = {
        localId: createLocalId(),
        primaryCreateNameId: "TODO",
        selected: false
    }
    const newCntvs: CreateNameTokenValue[] = []

    const newCns: CreateName[] = createNameLocalIdsByCreateNameGroupLocalId[createNameGroup.localId].map(oldCnLocalId => {
        const oldCn = createNamesByLocalId[oldCnLocalId]

        const newCn: CreateName = {
            ...oldCn,
            localId: createLocalId(),
            createNameGroupLocalId: newCng.localId,
        }

        if (createNameGroup.primaryCreateNameId == oldCn.localId) {
            newCng.primaryCreateNameId = newCn.localId
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

    return addNamesWithGroup(newCns, newCng, newCntvs)
}

export function addOrReplaceCreateNameTokenValue(newCreateNameTokenValue: CreateNameTokenValue, oldCreateNameTokenValue?: CreateNameTokenValue) {
    return makeActionObject(ADD_OR_REPLACE_CREATE_NAME_TOKEN_VALUE, { newCreateNameTokenValue, oldCreateNameTokenValue })
}

export function deleteCreateNameTokenValue(createNameTokenValue: CreateNameTokenValue) {
    return makeActionObject(DELETE_CREATE_NAME_TOKEN_VALUE, { createNameTokenValue })
}

export function setCreateNameCounterValue(createNameLocalId: string, counterValue: number | undefined) {
    return makeActionObject(SET_COUNTER_VALUE_FOR_CREATE_NAME, { createNameLocalId, counterValue })
}

export function setCreateNameGroupsSelected(createNameGroupIds: string[], selected: boolean) {
    return makeActionObject(SET_CREATE_NAME_GROUPS_SELECTED, { createNameGroupIds, selected })
}

export function deleteGroupsWithNames(createNameGroupIds: string[]) {
    return makeActionObject(DELETE_GROUPS_WITH_NAMES, { createNameGroupIds })
}

export function deleteSelectedGroupsWithNames() {
    return makeActionObject(DELETE_SELECTED_GROUPS_WITH_NAMES, {})
}


export type ACTION_TYPE =
    | ReturnType<typeof addNamesWithGroup>
    | ReturnType<typeof duplicateCreateNameGroup>
    | ReturnType<typeof addOrReplaceCreateNameTokenValue>
    | ReturnType<typeof deleteCreateNameTokenValue>
    | ReturnType<typeof setCreateNameCounterValue>
    | ReturnType<typeof setCreateNameGroupsSelected>
    | ReturnType<typeof deleteGroupsWithNames>
    | ReturnType<typeof deleteSelectedGroupsWithNames>
