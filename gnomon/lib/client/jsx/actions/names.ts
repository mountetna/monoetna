import { v4 as uuidv4 } from 'uuid';

import { CreateName, CreateNameGroup, Rule, TOKEN_VALUE_PLACEHOLDER } from "../models"
import { makeActionObject } from "./utils"



export const ADD_NAMES_WITH_GROUP = "ADD_NAMES_WITH_GROUP"


interface AddNamePayload {
    createNames: CreateName[]
    createNameGroup: CreateNameGroup
}


export function createNamesWithGroupForRule(primaryRuleName: string, allRules: Record<string, Rule>): AddNamePayload {
    const createNames: CreateName[] = []
    const ruleNamesToScan: string[] = [primaryRuleName]
    let primaryCreateNameId: string | undefined = undefined

    // make all CreateNames
    while (ruleNamesToScan.length) {
        const ruleName = ruleNamesToScan.pop()
        // typescript needs this check even though
        // we're already checking ruleNamesToScan.length
        if (!ruleName) {
            break
        }
        const rule = allRules[ruleName];
        const createName = {
            localId: uuidv4(),
            tokenValues: rule.tokenNames.map(_ => TOKEN_VALUE_PLACEHOLDER),
            ruleName,

        } as CreateName;

        createNames.push(createName)
        if (ruleName == primaryRuleName) {
            primaryCreateNameId = createName.localId
        }

        if (rule.parentRuleNames.length) {
            ruleNamesToScan.push(...rule.parentRuleNames)
        }
    }

    if (!primaryCreateNameId) {
        throw new Error("Error creating names with group. Missing primary name.")
    }

    return {
        createNames: createNames,
        createNameGroup: {
            localId: uuidv4(),
            createNameIds: createNames.map(cn => cn.localId),
            primaryCreateNameId: primaryCreateNameId,
            selected: false,
        },
    }
}


export function addNamesWithGroup(createNames: CreateName[], createNameGroup: CreateNameGroup) {
    return makeActionObject(ADD_NAMES_WITH_GROUP, { createNames, createNameGroup })
}


export type ACTION_TYPE =
    | ReturnType<typeof addNamesWithGroup>
