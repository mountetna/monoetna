export function defaultDict<K extends string, V>(createValue: (property: string) => any): Record<K, V> {
    return new Proxy(Object.create(null), {
        get(storage: object, property: string) {
            if (!(property in storage))
                storage[property] = createValue(property);
            return storage[property];
        }
    });
}