import { flattenObject, unflattenObject } from "./object"

const flattenDelimiter = '.'

export function parseSearchParams(searchParams: URLSearchParams): Record<string, any> {
    const parsed: Record<string, any> = {}

    for (const [k, v] of searchParams) {
        parsed[k] = JSON.parse(v)
    }

    return unflattenObject(parsed, flattenDelimiter)
}

export function toSearchParamsString(queryParams: Record<string, any>): string {
    const stringified: Record<string, string> = {}

    for (const [k, v] of Object.entries(flattenObject(queryParams, flattenDelimiter))) {
        stringified[k] = JSON.stringify(v)
    }

    return (new URLSearchParams(stringified)).toString()
}