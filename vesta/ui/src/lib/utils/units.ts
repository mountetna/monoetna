import { unit, format } from 'mathjs'


export interface SIValue {
    value: number
    SIUnitPrefix: string
}

export function roundValueToNearestSIUnit(value: number, precision: number = 3): SIValue {
    const simpleSIUnit = 'Hz'
    const formatted = format(unit(value, simpleSIUnit), precision).split(' ')

    return {
        value: Number(formatted[0]),
        SIUnitPrefix: formatted[1].replace(simpleSIUnit, '')
    }
}