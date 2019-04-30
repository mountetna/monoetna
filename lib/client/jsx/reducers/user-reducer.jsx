
const parseToken = (token) => {
  let [ header, params, signature ] = token.split(/\./);

  return JSON.parse(atob(params));
}

const user = (state, action) => {
  if (!state) state = { };

  switch(action.type){
    case 'ADD_USER':
      return {
        ...state,
        ...parseToken(action.token)
      }

    default:
      return state;
  }
};

export default user;
