export function parseSearchParams(searchParams: URLSearchParams): Record<string, any> {
    const parsed: Record<string, any> = {}

    searchParams.forEach((v, k) => {
        parsed[k] = JSON.parse(v)
    })

    return parsed
}

export function toSearchParamsString(queryParams: Record<string, any>): string {
    const stringified: Record<string, string> = {}

    Object.entries(queryParams).forEach(([k, v]) => {
        stringified[k] = JSON.stringify(v)
    })
    
    return (new URLSearchParams(stringified)).toString()
}