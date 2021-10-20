import React, { useState, useEffect, useCallback, useContext } from 'react';

type RedcapContextType = {
  schema: any,
  setSchema: Function
}

export const RedcapContext = React.createContext<RedcapContextType>({ schema: {}, setSchema: () => {} });

export const RedcapProvider = ({children}:{children:React.ReactNode}) => {
  const [ schema, setSchema ] = useState({});
  return (<RedcapContext.Provider value={{
    schema, setSchema
  }}>
    {children}
  </RedcapContext.Provider>);
}
