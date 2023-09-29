import { defaultDict } from "../utils/object";



// TODO: type to enforce T[idPropName]: string
export function listToIdObject<T>(list: T[], idPropName: string): Record<string, T> {
    const idObject = {};
    list.forEach(item => {
        idObject[item[idPropName]] = item
    })

    return idObject;
}

export function listToIdListObject<T>(list: T[], idPropName: string): Record<string, T[]> {
    const idObject = defaultDict<string, T[]>(_ => []);

    list.forEach(item => idObject[item[idPropName]].push(item))

    return { ...idObject } // remove default value getter
}