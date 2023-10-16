export interface TokenValue {
    localId: string
    name: string
    label: string
    tokenName: string
}

export const UNSET_VALUE = "UNSET"

export const UNSET_TOKEN_VALUE = {
    localId: UNSET_VALUE,
    name: UNSET_VALUE,
    label: UNSET_VALUE,
    tokenName: UNSET_VALUE,
}

export interface Token {
    name: string
    label: string
}

export interface Rule {
    name: string
    hasCounter: boolean
}

export interface RuleToken {
    localId: string
    ruleName: string
    tokenName: string
    ord: number
}

export interface RuleParent {
    localId: string
    ruleName: string
    parentRuleName: string
    ord: number
}

export interface Synonym {
    value: string
    tokenName: string
}

export interface CreateName {
    localId: string
    ruleCounterValue?: number
    ruleName: string
    createNameGroupLocalId: string
}

export interface CreateNameTokenValue {
    localId: string
    tokenValueLocalId: string
    createNameLocalId: string
    ruleTokenLocalId: string
}

export interface CreateNameGroup {
    localId: string
    primaryCreateNameLocalId: string
}