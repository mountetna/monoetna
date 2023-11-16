export function defaultDict<K extends keyof any, V>(createValue: (property: any) => V): Record<K, V> {
    return new Proxy(Object.create(null), {
        get(storage: Record<any, any>, property: any) {
            if (!(property in storage)) {
                storage[property] = createValue(property);
            }
            return storage[property];
        }
    });
}

// TODO: type to enforce idPropName in T
export function listToIdObject<T extends Record<any, any>>(list: T[], idPropName: string): Record<any, T> {
    const idObject: Record<any, T> = {};
    list.forEach(item => {
        idObject[item[idPropName]] = item
    })

    return idObject;
}

// TODO: type to enforce groupPropName in T and idPropName in T
export function listToIdGroupObject<T extends Record<any, any>>(list: T[], groupPropName: string, idPropName: string): Record<any, string[]> {
    const idObject = defaultDict<string, string[]>(_ => []);

    list.forEach(item => idObject[item[groupPropName]].push(item[idPropName]))

    return { ...idObject } // remove default value getter
}