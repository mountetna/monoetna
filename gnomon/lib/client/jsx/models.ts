// TODO: make this a module instead?
export interface Separator {
    text: string
}

export interface TokenBase {
    name: string
    index?: number
}

export interface RuleToken extends TokenBase {
    values: string[]
}

export interface SelectableToken extends TokenBase {
    value?: string
}

export interface Rule {
    elements: (RuleToken | Separator)[]
}

export interface CreateName {
    elements: (SelectableToken | Separator)[]
    rule_name: string
    selected: boolean
}