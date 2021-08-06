import {Maybe, some, withDefault} from "../../../../selectors/maybe";
import {useEffect} from "react";

export function useSetsDefault<T>(_default: T, value: Maybe<T>, onChange: (v: Maybe<T>) => void) {
  useEffect(() => {
    if (!value) onChange(some(_default));
  }, [_default, onChange, value])

  return withDefault(value, _default);
}