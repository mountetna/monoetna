const dialog = (state, action) => {
  if (!state) state = { };

  switch(action.type) {
    case 'SHOW_DIALOG':
      return {
        ...action.dialog
      };
    case 'DISMISS_DIALOG':
      return {};
    default:
      return state;
  }
};

export default dialog;

