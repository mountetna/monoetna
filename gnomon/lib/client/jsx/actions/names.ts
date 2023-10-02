import { v4 as uuidv4 } from 'uuid';

import { CreateName, CreateNameGroup, CreateNameTokenValue, RuleParent } from "../models"
import { makeActionObject } from "./utils"



export const ADD_NAMES_WITH_GROUP = "ADD_NAMES_WITH_GROUP"
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
): AddNamePayload {

    const createNames: CreateName[] = []
    const ruleNamesToScan: string[] = [primaryRuleName]

    const createNameGroup: CreateNameGroup = {
        localId: uuidv4(),
        primaryCreateNameId: "TODO",
        selected: false
    }

    // make all CreateNames
    while (ruleNamesToScan.length) {
        const ruleName = ruleNamesToScan.pop()
        // typescript needs this check even though
        // we're already checking ruleNamesToScan.length
        if (!ruleName) {
            break
        }
        const createName: CreateName = {
            localId: uuidv4(),
            ruleName,
            createNameGroupLocalId: createNameGroup.localId
        }

        createNames.push(createName)
        if (ruleName == primaryRuleName) {
            createNameGroup.primaryCreateNameId = createName.localId
        }

        if (ruleParentLocalIdsByRuleName[ruleName]) {
            const parentRuleNames = ruleParentLocalIdsByRuleName[ruleName].map(rpLocalId => {
                return ruleParentsByLocalId[rpLocalId].parentRuleName
            })
            ruleNamesToScan.push(...parentRuleNames)
        }
    }

    return addNamesWithGroup(createNames, createNameGroup)
}

export function addNamesWithGroup(createNames: CreateName[], createNameGroup: CreateNameGroup) {
    return makeActionObject(ADD_NAMES_WITH_GROUP, { createNames, createNameGroup })
}

export function duplicateCreateNameGroup(createNameGroupId: string) {
    return makeActionObject(DUPLICATE_CREATE_NAME_GROUP, { createNameGroupId })
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
