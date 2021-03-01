import {Debouncer} from "../debouncer";
import {delay} from "../../spec/helpers";

describe('Debouncer', () => {
  describe('in eager mode with 0 maxGating', () => {
    let debouncer = new Debouncer({ eager: true, maxGating: 0 });
    let cb = jest.fn();

    // Every ready call is synchronous
    debouncer.ready(cb);
    expect(cb).toHaveBeenCalledTimes(1);
    debouncer.ready(cb);
    expect(cb).toHaveBeenCalledTimes(2);
    debouncer.ready(cb);
    expect(cb).toHaveBeenCalledTimes(3);
  })

  describe('in eager mode', () => {
    it('fires immediately on any invocation having spent > window time since the last', async () => {
      let debouncer = new Debouncer({ eager: true, windowMs: 1000000 });
      let cb = jest.fn();

      let promise = debouncer.ready(cb);
      // Synchronous execution of cb
      expect(cb).toHaveBeenCalledTimes(1);
      // The promise should resolve immediately thanks to the eager fire.  Otherwise this test timesout and the
      // implementation is wrong.
      await promise;

      debouncer.windowMs = 500;


      promise = debouncer.ready(cb);
      // Bounced
      expect(cb).toHaveBeenCalledTimes(1);

      await delay(200);
      // Bounced
      debouncer.ready(cb);
      expect(cb).toHaveBeenCalledTimes(1);

      await delay(500);
      expect(cb).toHaveBeenCalledTimes(2);
      debouncer.ready(cb);

      // Eager synchronous callback again
      expect(cb).toHaveBeenCalledTimes(3);
    });
  })
})