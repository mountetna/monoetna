// TODO: make this a module instead?
export interface TokenValue {
    name: string
    label: string
    tokenName: string
}

export const TOKEN_VALUE_PLACEHOLDER: TokenValue = {
    name: "$$",
    label: "PLACEHOLDER",
    tokenName: "PLACEHOLDER",
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
    // this is only necessary if tokens
    // can appear on rules multiple times
    localId: string
    ruleName: string
    tokenName: string
    ord: number
}

export interface RuleParent {
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
    counterValue?: number
    ruleName: string
}

export interface CreateNameTokenValue {
    value: string
    tokenValueName: string
    createNameLocalId: string
    ruleTokenLocalId: string  // is this necessary?
}

export interface CreateNameGroup {
    localId: string
    primaryCreateNameId: string
    selected: boolean
}

export interface CreateNameGroupItem {
    createNameLocalId: string
    createNameGroupLocalId: string
}