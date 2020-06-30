export function assertIsSome(thing) {
  const keys = Object.keys(thing);

  if (keys.length === 1) {
    const key = keys[0];
    const value = thing[key];
    if (!value) {
      throw new Error(`Expected ${key} to be present, value was ${value}`);
    }
    return value;
  }

  if (!thing) {
    throw new Error(`Expected value, got ${thing}`);
  }
  return thing;
}