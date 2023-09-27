// TODO: make this a module instead?
export interface TokenValue {
    name: string
    label: string
}

export interface Token {
    name: string
    label: string
    values: TokenValue[]
}

export type Synonym = Array<string>

export interface Rule {
    name: string
    tokenNames: string[]
    hasCounter: boolean
    parentRuleNames: string[]
}

export const TOKEN_VALUE_PLACEHOLDER = "$$"

export interface CreateName {
    localId: string
    // string=tokenValue number=counterValue
    tokenValues: string[]
    counterValue?: number
    ruleName: string
}

export interface CreateNameGroup {
    localId: string
    // CreateName.localId[]
    // Cached here to prevent needing to search.
    // (Value will never change.)
    createNameIds: string[]
    primaryCreateNameId: string
    selected: boolean
}