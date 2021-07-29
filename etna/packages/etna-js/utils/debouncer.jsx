export class Debouncer {
  constructor({eager= false, windowMs = 100, maxGating = null} = {}) {
    this.windowMs = windowMs;
    this.lastTimeout = null;
    this.eager = eager;
    this.maxGating = maxGating;
    this.reset();
  }

  reset() {
    this.promise = new Promise((resolve) => {
      this.resolve = resolve;
    });
    if (this.lastTimeout) clearTimeout(this.lastTimeout);
    this.lastTimeout = null;
    this.lastFire = Date.now();
    this.f = null;
  }

  ready(f = null) {
    if (f) this.f = f;
    clearTimeout(this.lastTimeout);

    const { promise } = this;

    const fire = () => {
      const { resolve, f } = this;
      this.reset();

      if (f) f();
      resolve();
    };

    if (this.timeout() === 0) {
      fire();
    }

    this.lastTimeout = setTimeout(fire, this.timeout());

    return promise;
  }

  timeout() {
    if (this.lastTimeout == null && this.eager) {
      return 0;
    }

    const remainingGatingTime = this.maxGating != null ?
      Math.max(this.maxGating - Date.now() + this.lastFire, 0)
      : this.windowMs;
    return Math.min(remainingGatingTime, this.windowMs);
  }
}