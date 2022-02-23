// TODO: In the future, we'd want to convert etna to typescript, or atleast type its .d.ts files there, rather
// than each project having its own hand rolled .d.ts.

declare module 'etna-js/utils/fetch' {
  export function checkStatus(response: Response): any;

  export function handleFetchSuccess(response: Response): Promise<any>;

  export function handleFetchError(err: any): Promise<[string]>;

  export type HeaderType = 'json' | 'csrf' | 'auth';

  export function headers(...types: HeaderType[]): {[k: string]: string};

  export function isJSON(response: Response): boolean;

  export function json_get(path: string, params?: any): Promise<any>;

  export function json_post(path: string, params?: any): Promise<any>;

  export function json_delete(path: string): Promise<any>;
}

declare module 'etna-js/hooks/useActionInvoker' {
  export function useActionInvoker<
    T extends {type: string; [others: string]: any},
    R = any
  >(): (a: T) => R;
}

declare module 'etna-js/actions/message_actions' {
  export function showMessages(messages: string[]): {
    type: 'SHOW_MESSAGES';
    messages: string[];
  };
  export function dismissMessages(): {type: 'DISMISS_MESSAGES'};
}

declare module 'etna-js/utils/cancellable' {
  export type CancelOrResult<T> = {result?: T; cancelled?: true};
  export class Cancellable {
    race<T>(p: Promise<T>): Promise<{result?: T; cancelled?: true}>;
    cancel(): void;
    async run<Result, P>(
      gen: AsyncGenerator<Result, P>
    ): Promise<{result?: Result; cancelled?: true}>;
  }
  export function cancelledAsMaybe<T>(cancelled: {
    result?: T;
    cancelled?: true;
  }): [T] | null;
  export type AsyncGenerator<Result, P> = Generator<Promise<P>, Result, P>;
}

declare module 'etna-js/utils/cancellable_helpers' {
  import {Cancellable, AsyncGenerator} from 'etna-js/utils/cancellable';

  export function useWithContext(
    fn: (context: Cancellable) => void,
    deps: any[] = []
  );
  export function useAsync<Result>(
    fn: () => AsyncGenerator<Result, any>,
    deps: any[] = [],
    cleanup: () => void = () => null
  ): [[Result] | null, [any] | null];
  export function useAsyncCallback<Result, Args extends any[]>(
    fn: (...args: Args) => AsyncGenerator<Result, any>,
    deps: any[] = [],
    cleanup: () => void = () => null
  ): [
    (...args: Args) => Promise<{result?: Result; cancelled?: true}>,
    () => void
  ];

  export function* runAsync<T>(fn: () => Promise<T>): AsyncGenerator<T, T>;
  export function* runPromise<T>(v: Promise<T>): AsyncGenerator<T, T>;
  export function* runGen<T>(
    v: Generator<unknown, T, unknown>
  ): AsyncGenerator<T, T>;
}

declare module 'etna-js/utils/retryable' {
  import {AsyncGenerator} from 'etna-js/utils/cancellable';

  export async function withRetries<T>(
    b: () => Promise<T>,
    isRetryable?: (e: any) => number,
    maxRetries?: number
  ): Promise<T>;
  export function runAttempts<T>(
    b: () => Promise<T>,
    isRetryable?: (e: any) => number,
    maxRetries?: number
  ): AsyncGenerator<T, any>;
}

declare module 'etna-js/utils/semaphore' {
  export class Trigger<T = void> {
    public resolve(v: T): void;
    public reject(e: any): void;
    public promise: Promise<T>;
  }

  export class UnbufferedChannel<T = void> {
    closed: boolean;
    public send(v: T): void;
    public reject(e: any): void;
    receive(): Promise<T>;
    close(v: T): void;
  }

  export class BufferedChannel<T = void> extends UnbufferedChannel<T> {
    constructor(bufferSize = 0);
    drainPending(): Iterable<Promise<T>>;
    buffer: V[];
  }

  export class Semaphore {
    async ready<T>(work: () => Promise<T>): Promise<T>;
  }
}

declare module 'etna-js/utils/blob' {
  export function readTextFile(accept: string): Promise<string>;
  export function downloadBlob(params: {
    data: string;
    filename?: string;
    contentType?: string;
  }): void;
}

declare module 'javascript-color-gradient' {
  export = class Gradient {
    setMidpoint(point: number): void;

    setGradient(s: string, e: string): void;

    getColor(p: number): string;
  };
}

declare module 'etna-js/plots/models/vector' {
  export = class Vector<T = any> {
    constructor(items: {value: T; label: string | null}[]) {}
  };
}

declare module 'etna-js/plots/components/xy_plot/xy_plot' {
  export = any;
}

declare module 'etna-js/utils/colors' {
  export function autoColors(n: number): string[];
}

declare module 'etna-js/utils/markdown' {
  export = function (md: string): string {};
}

declare module 'etna-js/components/flat-button' {
  export = function FlatButton(params: {
    className?: string;
    icon?: string;
    disabled?: boolean;
    title?: string;
    label?: string;
    onClick?: () => void;
  }): any {};
}

declare module 'etna-js/components/Notifications' {
  import {ReactElement} from 'react';
  export function Notifications(props: any): ReactElement;
}

declare module 'etna-js/components/ModalDialogContainer' {
  import {ReactElement} from 'react';
  export function ModalDialogContainer(props: any): ReactElement;
}

declare module 'etna-js/components/messages' {
  import {ReactElement} from 'react';
  export = function Messages(props: any): ReactElement {};
}

declare module 'etna-js/components/icon' {
  export = function Icon(params: {
    className?: string;
    icon: string;
    disabled?: boolean;
    overlay?: string;
    title?: string;
    onClick?: () => void;
  }): any {};
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

declare module 'etna-js/components/inputs/dropdown_autocomplete_wrapper' {
  export default function DropdownAutocompleteInput(...args: any[]): any;
}

declare module 'etna-js/components/inputs/slow_text_input' {
  export default function SlowTextInput(...args: any[]): any;
}

declare module 'etna-js/components/inputs/text_input' {
  export default function TextInput(...args: any[]): any;
}

declare module 'etna-js/components/inputs/select_input' {
  type Value = string | {text: string; value: any};
  export default function SelectInput(props: {
    showNone?: 'disabled' | boolean;
    defaultValue?: any;
    values: Value[];
    onChange: (v: any) => void;
  }): any;
}

declare module 'etna-js/components/link' {
  export default function Link(p: {link: string; children?: any}): any;
}

declare module 'etna-js/hooks/useReduxState' {
  export function useReduxState(): any;
}

declare module 'etna-js/spec/helpers' {
  import {Store} from 'redux';

  export function delay(time: number): Promise<void>;

  export function mockStore(state: any): Store;

  export function joinedDeferredPromises(
    ...promiseChains: (Promise<() => Promise<any>> | (() => Promise<any>))[][]
  ): Promise<Promise[]>;

  export function stubUrl(params: {
    verb: 'get' | 'post' | 'delete' | 'put';
    path?: string;
    url?: string;
    host?: string;
    headers?: {[k: string]: string};
    status?: number;
    response:
      | Object
      | ((
          uri: string,
          body: any,
          cb: (err: any | null, result: [number, string | Buffer]) => void
        ) => void);
    request?: Object | string;
  }): Promise<void>;
}

declare module 'etna-js/selectors/magma' {
  export function selectModelNames(state: any): string[];
  export function selectModels(state: any): any;
  export function selectTemplate(state: any, model_name: string): any;
}

declare module 'etna-js/actions/magma_actions' {
  export function requestModels(): any;
  export function requestAnswer(question: any): any;
}

declare module 'etna-js/utils/debounce' {
  export default function debounce(
    func: any,
    wait: number,
    immediate?: boolean
  ): any;
}

declare module 'etna-js/utils/debouncer' {
  export class Debouncer {
    constructor(params: any);
    reset();
    ready(func: any);
    timeout();
  }
}

declare module 'etna-js/api/magma_api' {
  export function getAnswer(question: any, exchange: any): Promise<T>;
  export function getDocuments(doc_args: any, fetch: Function): Promise<T>;
}

declare module 'etna-js/utils/copy' {
  export function copyText(text: string): void;
}

declare module 'etna-js/actions/exchange_actions' {
  export class Exchange {
    dispatch: any;
    exchange_name: string;

    constructor(dispatch: any, exchange_name: string);

    fetch: Function;
  }
}

declare module 'etna-js/utils/tsv' {
  export type MatrixDatum = {[key: string]: any};
  export type MatrixData = MatrixDatum[];

  export function downloadTSV(
    data: MatrixData,
    fields: string[],
    fileName: string
  ): void;
}

declare module 'etna-js/utils/fetch';

declare module 'etna-js/contexts/magma-context';

declare module 'etna-js/components/ModalDialogContainer' {
  export function useModal(): {openModal: any; dismissModal: any};
}

declare module 'etna-js/utils/file' {
  export function filePath(parent: string, file_name: string): string;
  export function includesFolders(path: string | null): boolean;
}

declare module 'etna-js/actions/location_actions' {
  export function pushLocation(link: string): {
    type: 'UPDATE_LOCATION';
    link: string;
  };
}
