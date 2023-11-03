export async function asyncSetTimeout(waitIntervalMs: number) {
    await new Promise(resolve => setTimeout(resolve, waitIntervalMs))
}


export function createFnConcurrencyWrapper<TFunc extends (...args: any[]) => any>(
    func: TFunc,
    concurrencyLimit: number,
    waitIntervalMs: number = 100
): (...args: Parameters<TFunc>) => Promise<Awaited<ReturnType<typeof func>>> {

    let concurrentCalls = 0

    return async function (...args: Parameters<TFunc>) {
        while (concurrentCalls >= concurrencyLimit) {
            await new Promise(resolve => setTimeout(resolve, waitIntervalMs))
        }

        try {
            concurrentCalls += 1
            const res = func(...args)

            if (res instanceof Promise) {
                return await res
            }
            return res

        } finally {
            concurrentCalls -= 1
        }
    };
}