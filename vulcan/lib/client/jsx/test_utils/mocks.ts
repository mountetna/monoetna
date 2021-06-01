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

export type PromiseArgs<R> = [(r: R) => void, (e: any) => void];
export type AsyncMock<R> = [jest.Mock, PromiseArgs<R>[]];
export function asyncFn<R>(): AsyncMock<R> {
  const mock = jest.fn();
  const promiseArgs: PromiseArgs<R>[] = [];

  mock.mockImplementation((...args: any[]) => {
    return new Promise<R>((resolve, reject) => {
      promiseArgs.push([resolve, reject]);
    })
  });

  return [mock, promiseArgs];
}