import {act} from "react-test-renderer";

export function createFakeStorage(): Storage {
  const storage: { [k: string]: string } = {};

  return {
    setItem(key: string, value: string) {
      storage[key] = value + "";
    },
    get length() {
      return Object.keys(storage).length;
    },
    clear() {
      Object.keys(storage).forEach(k => delete storage[k]);
    },
    getItem(key: string): string | null {
      if (key in storage) return storage[key];
      return null;
    },
    removeItem(key: string) {
      delete storage[key];
    },
    key(index: number): string | null {
      return Object.keys(storage)[index];
    }
  }
}

export type PromiseArgs<R> = [(r: R) => Promise<void>, (e: any) => Promise<void>];
export class AsyncMock<R> {
  constructor(public jestMock: jest.Mock, private promiseArgs: PromiseArgs<R>[]) {
  }

  hasPendingRequest() {
    return this.promiseArgs.length > 0;
  }

  reset() {
    this.promiseArgs.length = 0;
    this.jestMock.mockClear();
  }

  async awaitCall(reset: boolean): Promise<[PromiseArgs<R>, any[]]> {
    if (reset) {
      this.reset();
    }

    while (!this.hasPendingRequest()) {
      await new Promise((resolve) => setTimeout(resolve, 1000));
    }

    const next = this.promiseArgs.shift();
    const call = this.jestMock.mock.calls.shift();
    if (next && call) {
      return [next, call];
    }

    throw new Error('Should never have been reached');
  }
}

export function asyncFn<R>(): AsyncMock<R> {
  const mock = jest.fn();
  const promiseArgs: PromiseArgs<R>[] = [];

  mock.mockImplementation((...args: any[]) => {
    return new Promise<R>((resolve, reject) => {
      promiseArgs.push([(v) => act(async () => resolve(v)), (e) => act(async () => reject(e))]);
    })
  });

  return new AsyncMock<R>(mock, promiseArgs);
}