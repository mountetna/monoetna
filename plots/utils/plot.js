export const empty = (i) => i == null || i == undefined || i == '';

export const validDomain = (min, max, vector) =>
  !empty(min) && !empty(max)
    ? [parseFloat(min), parseFloat(max)]
    : vector.values;
