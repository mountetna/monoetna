import { useDispatch as useDispatchBase, batch } from 'react-redux';



export const useDispatch = (...args: Parameters<typeof useDispatchBase>): (
    (dispatchee: Record<any, any> | Record<any, any>[]) => void
) => {
    const dispatch = useDispatchBase(...args);

    const _dispatch = (dispatchee: Record<any, any> | Record<any, any>[]) => {
        if (Array.isArray(dispatchee)) {
            batch(() => {
                for (const payload of dispatchee) {
                    dispatch(payload);
                }
            });
            return;
        }

        dispatch(dispatchee);
    };

    return _dispatch;
};