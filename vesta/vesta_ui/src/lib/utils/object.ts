// From https://www.30secondsofcode.org/js/s/flatten-unflatten-object/
// Flattening and unflattening objects might result in data loss--especially if the keys contain the delimiter.
// Always be careful about handling and transforming data.
export function flattenObject(obj: Record<string, any>, delimiter = '.', prefix = '', flattenArrays: boolean = false): Record<string, any> {
    return Object.keys(obj).reduce((acc: Record<string, any>, k) => {
        const pre = prefix.length ? `${prefix}${delimiter}` : '';
        if (
            typeof obj[k] === 'object' &&
            (!Array.isArray(obj[k]) || flattenArrays) &&
            obj[k] !== null &&
            Object.keys(obj[k]).length > 0
        )
            Object.assign(acc, flattenObject(obj[k], delimiter, pre + k));
        else acc[pre + k] = obj[k];
        return acc;
    }, {})
}

export function unflattenObject(obj: Record<string, any>, delimiter = '.'): Record<string, any> {
    return Object.keys(obj).reduce((res, k) => {
        k.split(delimiter).reduce(
            (acc: Record<string, any>, e, i, keys) =>
                acc[e] ||
                (acc[e] = isNaN(Number(keys[i + 1]))
                    ? keys.length - 1 === i
                        ? obj[k]
                        : {}
                    : []),
            res
        );
        return res;
    }, {})
}

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

export function listToIdObject<T extends Record<any, any>>(list: T[], idPropName: keyof T): Record<any, T> {
    const idObject: Record<any, T> = {};
    list.forEach(item => {
        idObject[item[idPropName]] = item;
    });

    return idObject;
}

export function listToGroupObject<T extends Record<any, any>>(list: T[], groupPropName: keyof T): Record<string, T[]> {
    const idObject = defaultDict<string, T[]>(_ => []);

    list.forEach(item => idObject[item[groupPropName]].push(item));

    return { ...idObject }; // removes default value getter
}