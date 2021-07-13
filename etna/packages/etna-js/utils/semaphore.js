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

    close(v) {
        this.trigger.resolve(v);
        this.closed = true;
    }
}

export async function allAcquired(...semaphores) {
    while (true) {
        const awaiting = semaphores.find(({open}) => !open);
        if (!awaiting) return;
        await awaiting.closed();
    }
}

export class Semaphore {
    channel = new UnbufferedChannel();
    open = true;

    async ready(work) {
          const {release} = this.acquire();
          try {
              return await work();
          } finally {
              release();
          }
    }

    async closed() {
        while (true) {
            if (!this.open) return;
            await this.channel.receive();
        }
    }

    // Probably could be more efficient.
    async acquire() {
        while (true) {
            const acquired = this.maybeAcquire();
            if (acquired) return acquired[0];
            await this.channel.receive();
        }
    }

    maybeAcquire() {
        if (this.open) {
            this.open = false;
            return [{ release: () => { this.open = true; this.channel.send(); } }];
        }

        return null;
    }
}