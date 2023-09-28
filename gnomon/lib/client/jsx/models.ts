// TODO: make this a module instead?
export interface TokenValue {
    name: string
    label: string
}

export const TOKEN_VALUE_PLACEHOLDER = { name: "$$", label: "PLACEHOLDER" } as TokenValue

export interface Token {
    name: string
    label: string
    // TODO make this TokenValue.name[]
    values: TokenValue[]
}

export type Synonym = Array<string>

export interface Rule {
    name: string
    tokenNames: string[]
    hasCounter: boolean
    parentRuleNames: string[]
}

export interface CreateName {
    localId: string
    // TODO make this TokenValue.name[]
    tokenValues: TokenValue[]
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