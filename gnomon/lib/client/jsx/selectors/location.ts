import { State } from "../store"



export const selectPathParts = (state: State) => state.location.path.split("/").filter(el => el.length)