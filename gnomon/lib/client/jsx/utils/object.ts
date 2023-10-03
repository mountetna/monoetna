export function defaultDict<K extends string, V>(createValue: (property: string) => any): Record<K, V> {
    return new Proxy(Object.create(null), {
        get(storage: object, property: string) {
            if (!(property in storage))
                storage[property] = createValue(property);
            return storage[property];
        }
    });
}

// TODO: type to enforce idPropName in T
export function listToIdObject<T>(list: T[], idPropName: string): Record<string, T> {
    const idObject = {};
    list.forEach(item => {
        idObject[item[idPropName]] = item
    })

    return idObject;
}

// TODO: type to enforce groupPropName in T and idPropName in T
export function listToIdGroupObject(list: object[], groupPropName: string, idPropName: string): Record<string, string[]> {
    const idObject = defaultDict<string, string[]>(_ => []);

    list.forEach(item => idObject[item[groupPropName]].push(item[idPropName]))

    return { ...idObject } // remove default value getter
}