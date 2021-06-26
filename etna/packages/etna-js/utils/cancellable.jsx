export function useCancelledOnDismount() {
  const [cancellable, _] = useState(() => new Cancellable());
  useEffect(() => cancellable.cancel());
  return cancellable;
}

export class Cancellable {
  constructor() {
    this.cancelledPromise = new Promise((resolve, reject) => {
      this.cancel = () => resolve({cancelled: true});
    });
  }

  race(promise) {
    return Promise.race([
      promise.then(result => ({result})),
      this.cancelledPromise,
    ]);
  }

  async run(gen) {
    if (gen instanceof Function) gen = gen();

    let {result, cancelled} = await this.race(Promise.resolve());
    let done = false;
    let value = null;

    while (!cancelled) {
      ({done, value} = gen.next(result));

      if (done) {
        return {result: value};
      }

      ({result, cancelled} = await this.race(value));
    }

    if (cancelled) {
      return {cancelled};
    }

    return {result};
  }
}