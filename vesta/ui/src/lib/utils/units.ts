export interface SIValue {
    rawValue: number
    value: number
    siUnitPrefix: string
}

const SI_PREFIXES = ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

export function roundValueToNearestSIPrefix(value: number, precision: number = 3): SIValue {
    const initValue = value
    const base = 1000
    let index = 0

    while (Math.abs(value) >= base && index < SI_PREFIXES.length - 1) {
        value /= base
        index++
    }

    return {
        rawValue: initValue,
        value: Number(Number(value.toFixed(0)).toPrecision(precision)),
        siUnitPrefix: SI_PREFIXES[index]
    }
}
