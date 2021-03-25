import * as Reselect from 'reselect';
import { isEditor } from '../utils/janus';

export const selectUser = ({user})=>( user || {} );

export const selectUserPermissions = Reselect.createSelector(
  selectUser,
  (user)=>{
    return ('permissions' in user) ? user.permissions : [];
  }
);

export const selectUserProjectRole = Reselect.createSelector(
  selectUserPermissions,
  (permissions) => {
    let perm = permissions[CONFIG.project_name];
    return perm ? perm.role : null
  }
);

export const selectIsEditor = Reselect.createSelector(
  selectUser,
  user => user.permissions && isEditor(user, CONFIG.project_name)
);
