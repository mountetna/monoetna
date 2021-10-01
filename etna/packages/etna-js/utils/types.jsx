import Vector from '../plots/models/vector';
import Matrix from '../plots/models/matrix';

export const isPrimitiveType = (value) =>
  typeof value === 'string' ||
  typeof value === 'number' ||
  typeof value === 'boolean' ||
  typeof value === 'undefined' ||
  typeof value === null;

export const isVector = (value) => value instanceof Vector;

export const isMatrix = (value) => value instanceof Matrix;

export const flatten = (array) => [...new Set([].concat(...array))];

export const mapObject = (o, f) =>
  Object.assign({}, ...Object.keys(o).map((k) => ({[k]: f(k, o[k])})));

export const getType = v => Object.prototype.toString.call(v).match(/^\[object (.*)\]$/)[1];
