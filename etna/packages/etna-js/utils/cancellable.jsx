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

    console.log('pausing here...')
    let {result, cancelled} = await this.race(Promise.resolve());
    console.log('going in')
    let done = false;
    let value = null;
    let error = null;

    while (!cancelled) {
      if (error) {
        const innerError = error[0];
        error = null;
        console.log('gen.error');
        ({done, value} = gen.error(innerError));
      } else {
        console.log('gen.next', result);
        ({done, value} = gen.next(result));
      }

      if (done) {
        console.log('done')
        return {result: value};
      }

      try {
        console.log('racing...', value);
        ({ result, cancelled } = await this.race(value));
        console.log('got', result);
      } catch (e) {
        error = [e];
        console.log('errrrr')
      }
    }

    if (cancelled) {
      console.log('cancelled')
      return {cancelled};
    }

    return {result};
  }
}
