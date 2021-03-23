import {MutableRefObject, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan";
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from "etna-js/actions/message_actions";

export const defaultApiHelpers = {
    isLoading: false,
    scheduleWork<T>(work: Promise<T>): Promise<T> {
        work.catch((e) => {
            console.error(e);
        })

        return work;
    }
}

export function useApi(
    state: MutableRefObject<VulcanState>,
    props: {} = {},
    ): typeof defaultApiHelpers {

    const invoke = useActionInvoker();

    const [workCount, setWorkCount] = useState(0);
    const workRef = useRef(workCount);

    return {
        isLoading: workCount > 0,
        scheduleWork<T>(work: Promise<T>): Promise<T> {
            setWorkCount(++workRef.current);

            work.catch(e => {
                if (!(e instanceof Array)) {
                    e = [`${e}`];
                }

                console.error(e);
                invoke(showMessages(e));
            }).finally(() => {
                setWorkCount(--workRef.current);
            });

            return work;
        }
    }
}