import { defaultDict } from "../utils/object";



// TODO: type to enforce T[idPropName]: string
export function listToIdObject<T>(list: T[], idPropName: string): Record<string, T> {
    const idObject = {};
    list.forEach(item => {
        idObject[item[idPropName]] = item
    })

    return idObject;
}

export function listToIdGroupObject(list: object[], groupPropName: string, idPropName: string): Record<string, string[]> {
    const idObject = defaultDict<string, string[]>(_ => []);

    list.forEach(item => idObject[item[groupPropName]].push(item[idPropName]))

    return { ...idObject } // remove default value getter
}