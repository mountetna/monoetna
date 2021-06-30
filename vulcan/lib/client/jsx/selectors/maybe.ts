export type Maybe<T> = null | [T];

export function withDefault<T>(v: Maybe<T>, d: T) {
  if (v == null) return d;
  return v[0];
}

export function isSome<T>(a: Maybe<T>): a is [T] {
  return a != null;
}

export function some<T>(v: T): Maybe<T> {
  return [v];
}

export function alternate<T>(a: Maybe<T>, ...alts: Maybe<T>[]): Maybe<T> {
  if (isSome(a)) return a;

  for(let alt of alts) {
    if (isSome(alt)) return alt;
  }

  return null;
}

export function maybeOfNullable<T>(v: T | null | undefined): Maybe<T> {
  if (v == null) return null;
  return [v];
}

export function mapSome<A, B>(a: Maybe<A>, f: (a: A) => B): Maybe<B> {
  if (isSome(a)) {
    return [f(a[0])];
  }

  return null;
}

export function bindSome<A, B>(a: Maybe<A>, f: (a: A) => Maybe<B>): Maybe<B> {
  if (isSome(a)) {
    return f(a[0]);
  }

  return null;
}

export function applySome<A, B>(f: Maybe<(a: A) => B>, a: Maybe<A>): Maybe<B> {
  if (isSome(f)) {
    if (isSome(a)) {
      return [f[0](a[0])];
    }
  }

  return null;
}