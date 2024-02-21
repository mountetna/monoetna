import { TypedUseSelectorHook, useSelector } from 'react-redux';
import type { State } from './store';



// Use throughout your app instead of plain `useSelector`
export const useAppSelector: TypedUseSelectorHook<State> = useSelector;