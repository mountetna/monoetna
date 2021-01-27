// A subscription is just a collection of callbacks that can be invoked
// together with 'end'
// It's main purpose to collect clean up tasks for asynchronous work so that tests can
// clean up between runs and prevent test spoilage.
export class Subscription {
  constructor() {
    this.callbacks = [];
  }

  // Returns a wrapped callback that will behave similarly to the original cb, adds cleanup to this subscription
  // so that when end is called, it will no longer pass invocations on to the original cb.
  addSubscribedCallback(cb) {
    let cancelled = false;
    function wrapped() {
      if (cancelled) return;
      return cb.apply(this, arguments);
    }
    this.addCleanup(() => cancelled = true);
    return wrapped;
  }

  // Adds an arbitrary cleanup function.
  addCleanup(f) {
    if (f == null) return;
    if (f instanceof Subscription) {
      this.callbacks.push(() => f.end());
    } else {
      this.callbacks.push(f);
    }
  }

  end() {
    const {callbacks} = this;
    this.callbacks = [];
    callbacks.forEach(cb => {
      try {
        cb();
      } catch(e) {
        console.error(e);
      }
    });
  }
}

// A top level subscription for application level work to be cancelled
// Used primarily by tests.
export const appSubscription = new Subscription();