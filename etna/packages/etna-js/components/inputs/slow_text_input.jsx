import React, {Component, useCallback, useEffect, useRef, useState} from 'react';
import debounce from '../../utils/debounce';
import {Debouncer} from "../../utils/debouncer";

export default function SlowTextInput(props) {
  const {onChange, waitTime, eager, defaultValue, followDefault, ...inputProps} = props;
  const inputRef = useRef();
  const [value, setValue] = useState(defaultValue == null ? '' : defaultValue);
  const [prevDefault, setPrevDefault] = useState(() => defaultValue);
  const [debouncer, setDebouncer] = useState(() => new Debouncer({windowMs: waitTime, eager}));
  // Clear the existing debouncer and accept any new changes to the settings
  useEffect(() => {
    const debouncer = new Debouncer({windowMs: waitTime, eager})
    setDebouncer(debouncer);
    return () => debouncer.reset();
  }, [waitTime, eager]);

  const onChangeWithDebounce = useCallback(() => {
    const v = inputRef.current.value;
    debouncer.ready(() => onChange(v));
    setValue(v);
  }, [onChange, debouncer]);

  // When the default value changes, follow it
  useEffect(() => {
    if (followDefault && defaultValue !== prevDefault) {
      debouncer.reset();
      setValue(defaultValue);
      setPrevDefault(defaultValue);
    }
  }, [defaultValue, followDefault, debouncer, setPrevDefault, prevDefault]);


  return <input
    type='text'
    ref={inputRef}
    onChange={onChangeWithDebounce}
    value={value}
    {...inputProps}
    />;
}
