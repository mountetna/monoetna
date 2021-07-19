export class Trigger {
  resolve = () => null;
  reject = () => null;
  promise = new Promise((resolve, reject) => {
    this.resolve = resolve;
    this.reject = reject;
  });

  constructor() {
  }
}

export class UnbufferedChannel {
  trigger = new Trigger()
  closed = false;

  send(v) {
    if (this.closed) return;
    this.trigger.resolve(v);
    this.trigger = new Trigger();
  }

  reject(e) {
    if (this.closed) return;
    this.trigger.reject(e);
    this.trigger = new Trigger();
  }

  async receive() {
    return await this.trigger.promise;
  }

  close() {
    this.trigger.reject(new Error('Channel was closed'));
    this.closed = true;
  }
}

export class BufferedChannel extends UnbufferedChannel {
  sendSemaphore = new Semaphore();
  receiveSemaphore = new Semaphore();

  close() {
    this.sendSemaphore.close();
    this.receiveSemaphore.close();
    return super.close();
  }

  async send(v) {
    const { release } = await this.sendSemaphore.acquire();
    try {
      await this.receiveSemaphore.closed();
      return await super.send(v);
    } finally {
      release();
    }
  }

  async receive() {
    const { release } = await this.receiveSemaphore.acquire();

    try {
      return await super.receive();
    } finally {
      release();
    }
  }

  pending() {
    const self = this;

    function* genPending() {
      while (!self.sendSemaphore.open && self.receiveSemaphore.open) {
        yield self.receive();
      }
    }

    return genPending();
  }
}

export async function allAcquired(...semaphores) {
  while (true) {
    const awaiting = semaphores.find(({ open }) => !open);
    if (!awaiting) return;
    await awaiting.closed();
  }
}

export class Semaphore {
  channel = new UnbufferedChannel();
  open = true;

  close() {
    this.channel.close();
  }

  async ready(work) {
    const { release } = this.acquire();
    try {
      return await work();
    } finally {
      release();
    }
  }

  async closed() {
    while (!this.channel.closed) {
      if (!this.open) return;
      await this.channel.receive();
    }
  }

  // Probably could be more efficient.
  async acquire() {
    while (!this.channel.closed) {
      const acquired = this.maybeAcquire();
      if (acquired) return acquired[0];
      await this.channel.receive();
    }
  }

  maybeAcquire() {
    if (this.open) {
      this.open = false;
      this.channel.send();

      return [{
        release: () => {
          this.open = true;
          this.channel.send();
        }
      }];
    }

    return null;
  }
}