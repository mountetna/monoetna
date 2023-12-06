export function makeActionObject<T extends string, P>(type: T, payload: P): { type: T } & P {
    return { type, ...payload };
}