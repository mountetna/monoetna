import {
  UPDATE_UPLOADS,
} from "../workers/uploader";

// Uploader manages its own state.  Largely we are just reflecting that state in the host reducer.
const uploads = (old_uploads, action) => {
  switch(action.type) {
    case UPDATE_UPLOADS:
      return action.uploads;
  }
};

export default uploads;
