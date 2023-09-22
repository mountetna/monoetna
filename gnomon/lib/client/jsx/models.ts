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

export interface Synonym {
    [index: number]: string
}

export interface Counter {}

export interface RuleRef {
    value: string
}

export interface Rule {
    name: string
    elements: (Rule | Token | Counter)[]
}

export interface CreateName {
    // string=token_value number=index
    elements: (string | number)[]
    rule_name: string
    selected: boolean
}