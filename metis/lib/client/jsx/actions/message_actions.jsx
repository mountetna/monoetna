export const message = (message_type, title, message) => ({
  type: 'SHOW_DIALOG',
  dialog: { type: 'message', title, message, message_type }
});

export const errorMessage = (dispatch, message_type, title, message_handler) =>
  response => {
    if (response instanceof Promise) {
response.then(
        ({error}) => dispatch(message(message_type, title, message_handler(error)))
      );
}
    else {dispatch(message(message_type, title, message_handler(response)));}
  };
