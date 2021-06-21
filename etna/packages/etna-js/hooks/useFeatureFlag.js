import {useState} from 'react'
import {useReduxState} from "./useReduxState";

export function useFeatureFlag(flagName) {
  const user = useReduxState(({ user }) => user);

  return useState(() => {
    if (user && user.flags && user.flags.includes(flagName)) return true;
    return location.search.includes(`${flagName}=1`);
  })[0];
}