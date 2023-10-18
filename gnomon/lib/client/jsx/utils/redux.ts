import { useDispatch as useDispatchBase, batch } from 'react-redux'



export const useDispatch = (...args: any) => {
    const dispatch = useDispatchBase(...args)

    const _dispatch = (dispatchee: any) => {
        if (Array.isArray(dispatchee)) {
            batch(() => {
                for (const payload of dispatchee) {
                    dispatch(payload)
                }
            })
            return

        }

        dispatch(dispatchee)
    }

    return _dispatch
}