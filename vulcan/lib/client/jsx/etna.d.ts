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

declare module 'etna-js/hooks/useActionInvoker' {
    export function useActionInvoker<T extends { type: string }, R = any>(): (a: T) => R;
}

declare module 'etna-js/actions/message_actions' {
    export function showMessages(messages: string[]): { type: 'SHOW_MESSAGES', messages: string[] };
}

declare module 'etna-js/utils/cancellable' {
    export class Cancellable {
        race<T>(p: Promise<T>): Promise<{ result?: T, cancelled?: true }>;
        cancel(): void;
    }
}

declare module 'javascript-color-gradient' {
    export = class Gradient {
        setMidpoint(point: number): void
        setGradient(s: string, e: string): void;
        getColor(p: number): string;
    }
}

declare module 'etna-js/plots/models/vector' {
    export = class Vector<T = any> {
        constructor(items: {value: T, label: string | null}[]) {
        }
    }
}

declare module 'etna-js/plots/components/xy_plot/xy_plot' {
    export = any;
}

declare module 'etna-js/utils/colors' {
    export function autoColors(n: number): string[];
}

declare module 'etna-js/utils/markdown' {
    export = function(md: string): string {};
}

declare module 'etna-js/components/icon' {
    export = function Icon(params: {className?: string, icon: string, disabled?: boolean, overlay?: string, title?: string, onClick?: () => void}): any {};
}

declare module 'etna-js/components/inputs/numeric_input' {
    export function FloatInput(...args: any[]): any;
    export function IntegerInput(...args: any[]): any;
}

declare module 'etna-js/components/inputs/list_input' {
    export default function ListInput(...args: any[]): any;
}

declare module 'etna-js/components/inputs/dropdown_input' {
    export default function DropdownInput(...args: any[]): any;
}

declare module 'etna-js/components/inputs/dropdown_autocomplete' {
    export default function DropdownAutocomplete(...args: any[]): any;
}

declare module 'etna-js/components/inputs/slow_text_input' {
    export default function SlowTextInput(...args: any[]): any;
}
