import {InputSpecification} from './input_types';
import {Dispatch, SetStateAction, useEffect, useRef, useState} from 'react';

export function useBufferedInputState<T>(
  input: InputSpecification,
  d: T
): [T, Dispatch<SetStateAction<T>>] {
  const {value: inputValue} = input;
  const [value, setValue] = useState(inputValue === null ? d : inputValue);

  // Using a ref here safely skips the eslint check for the useEffect, allowing us to ONLY update when the input.default
  // value changes, and not get stuck in a loop.
  const ref = useRef(value);
  ref.current = inputValue;

  useEffect(() => {
    if (null != inputValue && inputValue !== ref.current) {
      setValue(inputValue);
    }
  }, [inputValue, setValue]);

  return [value, setValue];
}
