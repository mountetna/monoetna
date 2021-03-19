// TODO: In the future, we'd want to convert etna to typescript, or atleast type its .d.ts files there, rather
// than each project having its own hand rolled .d.ts.

declare module 'etna-js/utils/fetch' {
    export function checkStatus(response: Response): any;
    export function handleFetchSuccess(response: Response): Promise<any>;
    export function handleFetchError(err: any): Promise<[string]>;
    export type HeaderType = 'json' | 'csrf' | 'auth';
    export function headers(...types: HeaderType[]): {[k: string]: string};
    export function isJSON(response: Response): boolean;
}