import {BufferedChannel, Trigger} from "etna-js/utils/semaphore";

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

export function makeBlockingAsyncMock<A extends any[], T>(f: (...a: A) => Promise<T>): {
  channel: BufferedChannel<[A, Trigger<T>]>
  mock: (...a: A) => Promise<T>,
} {

  const channel = new BufferedChannel<[A, Trigger<T>]>();

  return {
    channel,
    async mock(...a: A) {
      const trigger = new Trigger<T>();
      console.log('awaiting to send through the mock')
      await channel.send([a, trigger]);
      console.log('awaiting the response')
      return await trigger.promise;
    }
  }
}

export function countIter(iter: Iterable<any>): number {
  let i = 0;
  for (let elem of iter) {
    i++;
  }

  return i;
}