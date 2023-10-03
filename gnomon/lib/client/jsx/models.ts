// TODO: make this a module instead?
export interface TokenValue {
    name: string
    label: string
    tokenName: string
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
    counterValue?: number
    ruleName: string
    createNameGroupLocalId: string
}

export interface CreateNameTokenValue {
    localId: string
    tokenValueName: string
    createNameLocalId: string
    ruleTokenLocalId: string
}

export interface CreateNameGroup {
    localId: string
    primaryCreateNameId: string
    selected: boolean
}