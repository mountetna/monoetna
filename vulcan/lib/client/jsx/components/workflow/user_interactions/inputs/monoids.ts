import {DataEnvelope} from "./input_types";

export type Monoid<DataElement, C extends DataElement> = {
  unit: C,
  concat(a: DataElement, b: DataElement): C,
};

export function flattener<DataElement, C extends DataElement = DataElement>({unit, concat}: Monoid<DataElement, C>) {
  return function flattenData(data: DataEnvelope<DataElement> | undefined | null): C {
    if (!data) return unit;
    return Object.values(data).reduce<C>(concat, unit);
  }
}

const collator = new Intl.Collator(undefined, {
  numeric: true,
  sensitivity: 'base'
});
export type StringOptions = string | string[];
export const flattenStringOptions = flattener<StringOptions, string[]>({
  unit: [],
  concat(a, b) {
    return [
      ...(typeof a === "string" ? [a] : a),
      ...(typeof b === "string" ? [b] : b),
    ].sort(collator.compare)
  }
})

export function defaultSelector<T>(unit: T) {
  return function selectDefault(data: DataEnvelope<T> | null | undefined): T {
    if (data) {
      for (let k in data) {
        return data[k];
      }
    }

    return unit;
  }
}

export const selectDefaultNumber = defaultSelector(0);
export const selectDefaultString = defaultSelector("");
export const selectDefaultBoolean = defaultSelector(false);

export function flattenNesting<T>(data: DataEnvelope<DataEnvelope<T>> | undefined | null): DataEnvelope<T> {
  if (!data) return {};
  return Object.values(data).reduce<DataEnvelope<T>>((a, b) => ({...a, ...b}), {} as DataEnvelope<T>);
}
