import {VulcanAction} from "../actions/vulcan_actions";
import VulcanReducer, {defaultVulcanState, VulcanState} from "../reducers/vulcan_reducer";
import {renderHook, act} from "@testing-library/react-hooks";
import {Dispatch} from "react";


export function stateFromActions<T = void>(actions: VulcanAction[],
                                           hook?: (dispatch: Dispatch<VulcanAction>,
                                                   state: VulcanState) => T): { state: VulcanState, dispatch: (action: VulcanAction) => VulcanState } {
    let state = defaultVulcanState;

    let rerender: (() => void) | null = null;
    const dispatch = (action: VulcanAction) => act(() => {
        state = VulcanReducer(state, action);
        if (rerender) rerender();
    });

    if (hook) {
        ({rerender} = renderHook(() => hook(dispatch, state)));
    }

    for (let action of actions) {
        state = VulcanReducer(state, action);
        if (rerender) rerender();
    }

    return {
        state, dispatch(action) {
            state = VulcanReducer(state, action);
            if (rerender) rerender();
            return state;
        }
    };
}