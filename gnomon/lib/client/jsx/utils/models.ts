import { v4 as uuidv4 } from "uuid";


export const createLocalId = (): string => uuidv4()


export type Status = "idle" | "inProgress" | "success" | "error"