export function listToIdObject<T>(list: T[], idPropName: string): Record<string, T> {
    const idObject = {};
    list.forEach(item => idObject[item[idPropName]] = item)
    
    return idObject;
}