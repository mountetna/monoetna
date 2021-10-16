import React, { useState, useEffect, useCallback, useContext } from 'react';

export const RedcapContext = React.createContext();

export const RedcapProvider = (props) => {
  const [ schema, setSchema ] = useState({});
  return (<RedcapContext.Provider value={{
    schema, setSchema
  }}>
    {props.children}
  </RedcapContext.Provider>);
}
