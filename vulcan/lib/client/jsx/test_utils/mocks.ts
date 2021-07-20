import {BufferedChannel, Trigger, UnbufferedChannel} from "etna-js/utils/semaphore";

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

class AsyncMock<A extends any[], T, C extends UnbufferedChannel<[A, Trigger<T>]>> {
  constructor(f: (...a: A) => Promise<T>,
    public channel: C,
    protected scheduler: (f: () => Promise<void>) => Promise<void> = async f => await f()) {
  }

  mock = async (...a: A) => {
    const trigger = new Trigger<T>();
    await this.channel.send([a, trigger]);
    return await trigger.promise;
  }

  async respond(f: (...a: A) => Promise<T> | T, q: Iterable<Promise<[A, Trigger<T>]>> = [this.channel.receive()]): Promise<void> {
    await this.scheduler(async () => {
      for (let request of q) {
        const [a, trigger] = await request;
        await Promise.resolve(f(...a)).then(trigger.resolve);
      }
    });
  }

  async reject(f: (...a: A) => never, q: Iterable<Promise<[A, Trigger<T>]>> = [this.channel.receive()]): Promise<void> {
    await this.scheduler(async () => {
      for (let request of q) {
        const [a, trigger] = await request;
        try {
          await Promise.resolve(f(...a)).catch(trigger.reject);
        } catch (e) {
          trigger.reject(e);
        }
      }
    });
  }
}

export class BlockingAsyncMock<A extends any[], T> extends AsyncMock<A, T, BufferedChannel<[A, Trigger<T>]>> {
  constructor(f: (...a: A) => Promise<T>, scheduler: (f: () => Promise<void>) => Promise<void> = async f => await f()) {
    super(f, new BufferedChannel(Infinity), scheduler);
  }

  pendingCount() {
    return this.channel.buffer.length;
  }

  async respond(f: (...a: A) => Promise<T> | T): Promise<void> {
    const {channel} = this;
    await super.respond(f, function* () {
      let found = false;
      for (let request of channel.drainPending()) {
        yield request;
        found = true;
      }

      if (!found) yield channel.receive();
    }());
  }

  async reject(f: (...a: A) => never): Promise<void> {
    const {channel} = this;
    await super.reject(f, function* () {
      let found = false;
      for (let request of channel.drainPending()) {
        yield request;
        found = true;
      }

      if (!found) yield channel.receive();
    }());
  }
}

export class UnbufferedAsyncMock<A extends any[], T> extends AsyncMock<A, T, UnbufferedChannel<[A, Trigger<T>]>> {
  constructor(f: (...a: A) => Promise<T>, scheduler: (f: () => Promise<void>) => Promise<void> = async f => await f()) {
    super(f, new UnbufferedChannel(), scheduler);
  }
}

export function countIter(iter: Iterable<any>): number {
  let i = 0;
  for (let elem of iter) {
    i++;
  }

  return i;
}