import {Maybe, some, withDefault} from '../../../selectors/maybe';
import {useEffect} from 'react';

export function useSetsDefault<T>(_default: T, value: Maybe<T>, onChange: (v: Maybe<T>) => void, key: string) {
  useEffect(() => {
    // console.log(`useSetsDefault onChange triggered for ${key}`)
    if (!value) {
      onChange({[key]: some(_default)});
    } else if (!(key in value) || (value[key] == null && _default != null)) {
      const newVals = {...value};
      newVals[key] = some(_default);
      onChange(newVals);
    }
  }, [_default, onChange, value, key]);

  return key in value ? withDefault(value[key], _default) : _default;
}