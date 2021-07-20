export function cancelledAsMaybe(v) {
  if (!v) return null;
  if ('cancelled' in v && v.cancelled) return null;
  if (!('result' in v)) return null;
  return [v.result];
}

export class Cancellable {
  constructor() {
    this.cancelledPromise = new Promise((resolve, reject) => {
      this.cancel = () => resolve({cancelled: true});
    });
  }

  race(v) {
    return Promise.race([
      Promise.resolve(v).then(result => ({result})),
      this.cancelledPromise,
    ]);
  }

  async run(gen) {
    if (gen instanceof Function) gen = gen();

    let {result, cancelled} = await this.race(Promise.resolve());
    let done = false;
    let value = null;
    let error = null;

    while (!cancelled) {
      if (error) {
        const innerError = error[0];
        error = null;
        ({done, value} = gen.throw(innerError));
      } else {
        ({done, value} = gen.next(result));
      }

      if (done) {
        return {result: value};
      }

      try {
        ({ result, cancelled } = await this.race(value));
      } catch (e) {
        error = [e];
      }
    }

    if (cancelled) {
      return {cancelled};
    }

    return {result};
  }
}
