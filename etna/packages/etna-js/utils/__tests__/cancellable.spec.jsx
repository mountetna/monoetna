import {Cancellable} from "../cancellable";

describe('Cancellable', () => {
  function* yieldThese(numbers) {
    const result = [];
    for (let i = 0; i < numbers.length; ++i) {
      result.push(yield Promise.resolve(numbers[i]));
    }

    return result;
  }

  it('can be cancelled before making any iterations', async () => {
    const cancellable = new Cancellable();

    async function runTest() {
      let fn = jest.fn();

      function* gen() {
        fn();
      }

      await cancellable.run(gen());
      return fn;
    }


    let genInner = await runTest();
    expect(genInner).toHaveBeenCalled();

    cancellable.cancel();
    genInner = await runTest();
    expect(genInner).not.toHaveBeenCalled();
  });

  it('does not continue after being cancelled', async () => {
    async function runTest(stopAfter) {
      let fn = jest.fn();
      const cancellable = new Cancellable();

      function* gen() {
        for (let v of [1, 2, 3, 4, 5, 6, 7]) {
          fn(yield Promise.resolve(v).then(v => {
            if (--stopAfter < 0) cancellable.cancel();
            return v;
          }));
        }

        return 'done';
      }

      fn(await cancellable.run(gen()));
      return fn;
    }


    let genInner = await runTest(3);
    expect(genInner).toHaveBeenCalledWith(1);
    expect(genInner).toHaveBeenCalledWith(2);
    expect(genInner).toHaveBeenCalledWith(3);
    expect(genInner).not.toHaveBeenCalledWith(4);
    expect(genInner).toHaveBeenCalledWith({cancelled: true});
    expect(genInner).not.toHaveBeenCalledWith({result: 'done'});

    genInner = await runTest(100);
    expect(genInner).toHaveBeenCalledWith(1);
    expect(genInner).toHaveBeenCalledWith(2);
    expect(genInner).toHaveBeenCalledWith(3);
    expect(genInner).toHaveBeenCalledWith(4);
    expect(genInner).toHaveBeenCalledWith(5);
    expect(genInner).toHaveBeenCalledWith(6);
    expect(genInner).toHaveBeenCalledWith(7);
    expect(genInner).not.toHaveBeenCalledWith({cancelled: true});
    expect(genInner).toHaveBeenCalledWith({result: 'done'});
  });

  it('can await all registered work to finish', async () => {
    async function runTest(stopAfter) {
      let fn = jest.fn();
      const cancellable = new Cancellable();

      function* gen() {
        for (let v of [1, 2, 3, 4, 5, 6, 7]) {
          fn(yield Promise.resolve(v).then(v => {
            if (--stopAfter < 0) cancellable.cancel();
            return v;
          }));
        }

        return 'done';
      }

      fn(await cancellable.run(gen()));
      return fn;
    }


    let genInner = await runTest(3);
    expect(genInner).toHaveBeenCalledWith(1);
    expect(genInner).toHaveBeenCalledWith(2);
    expect(genInner).toHaveBeenCalledWith(3);
    expect(genInner).not.toHaveBeenCalledWith(4);
    expect(genInner).toHaveBeenCalledWith({cancelled: true});
    expect(genInner).not.toHaveBeenCalledWith({result: 'done'});

    genInner = await runTest(100);
    expect(genInner).toHaveBeenCalledWith(1);
    expect(genInner).toHaveBeenCalledWith(2);
    expect(genInner).toHaveBeenCalledWith(3);
    expect(genInner).toHaveBeenCalledWith(4);
    expect(genInner).toHaveBeenCalledWith(5);
    expect(genInner).toHaveBeenCalledWith(6);
    expect(genInner).toHaveBeenCalledWith(7);
    expect(genInner).not.toHaveBeenCalledWith({cancelled: true});
    expect(genInner).toHaveBeenCalledWith({result: 'done'});
  });
});