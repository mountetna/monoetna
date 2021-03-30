import {VulcanAction} from "../actions/vulcan";
import VulcanReducer, {defaultVulcanState, VulcanState} from "../reducers/vulcan_reducer";
import {renderHook} from "@testing-library/react-hooks";


export function stateFromActions<T = void>(actions: VulcanAction[], hook?: (state: VulcanState) => T): VulcanState {
    let state = defaultVulcanState;
    const { rerender } = hook ? renderHook(() => hook(state)) : { rerender : null };

    for (let action of actions) {
        state = VulcanReducer(state, action);
        if (rerender) rerender();
    }

    return state;
}