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

export interface TokenRef {
    tokenName: string
}

export type Synonym = Array<string>

export interface RuleRef {
    ruleName: string
}

export interface Rule {
    name: string
    tokens: TokenRef[]
    hasCounter: boolean
    parents: RuleRef[]
}

export const ELEMENT_PLACEHOLDER = "$$"

export interface CreateName {
    localId: string
    // string=tokenValue number=counterValue
    tokenValues: string[]
    counterValue?: number
    ruleName: string
    parents: string[]  // CreateName.localId[]
}

export interface CreateNameGroup {
    localId: string
    // CreateName.localId[]
    // Cached here to prevent needing to search.
    // (Value will never change.)
    createNames: string[]
    selected: boolean
}