import React, { useState, useEffect, useCallback, useContext } from 'react';

type SchemaContextType = {
  schema: any,
  setSchema: Function
}

export const SchemaContext = React.createContext<SchemaContextType>({ schema: {}, setSchema: () => {} });

export const SchemaProvider = ({children}:{children:React.ReactNode}) => {
  const [ schema, setSchema ] = useState({});
  return (<SchemaContext.Provider value={{
    schema, setSchema
  }}>
    {children}
  </SchemaContext.Provider>);
};
