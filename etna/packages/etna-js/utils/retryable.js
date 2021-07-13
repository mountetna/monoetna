import {runPromise} from "./cancellable_helpers";

export function delay(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

export async function withRetries(asyncFn, isRetryable, maxRetries = 5) {
  for (let i = 0; i < maxRetries; ++i) {
    try {
      return await asyncFn();
    } catch(e) {
      const retryable = isRetryable ? isRetryable(e) : true;
      if (retryable === 0) {
        throw e;
      }

      if (retryable < 0) {
        i -= 1;
      }

      await delay(1000 * (i + 1) ** 2);
    }
  }

  return await asyncFn();
}

export function* runAttempts(asyncFn, isRetryable, maxRetries = 5) {
  for (let i = 0; i < maxRetries; ++i) {
    try {
      return yield* runPromise(b());
    } catch(e) {
      const retryable = isRetryable ? isRetryable(e) : null;
      if (retryable === 0) {
        throw e;
      }

      if (retryable < 0) {
        i -= 1;
      }

      yield delay(1000 * (i + 1) ** 2);
    }
  }


  return yield* runPromise(b());
}
