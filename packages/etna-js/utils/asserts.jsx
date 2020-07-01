export function assertIsSome(thing) {
  const keys = Object.keys(thing);

  if (keys.length === 1) {
    const value = thing[keys[0]];
    if (value == null) {
      throw new Error(`Expected ${key} to be present, value was ${value}`);
    }
    return value;
  } else if (keys.length > 1) {
    keys.forEach(key => {
      const value = thing[key];
      if (value == null) {
        throw new Error(`Expected ${key} to be present, value was ${value}`);
      }
    })
  }

  if (thing == null) {
    throw new Error(`Expected value, got ${thing}`);
  }
  return thing;
}