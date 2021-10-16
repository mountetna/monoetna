
import React, { useContext, useState } from 'react';


export const MagmaContext = React.createContext();

export const MagmaProvider = (props) => {
  const [ models, setModels ] = useState<MagmaModels | {}>({});

  return (<MagmaContext.Provider value={{
    models,
    setModels
  }}>
    {props.children}
  </MagmaContext.Provider>);
}
