// taken from https://github.com/vercel/next.js/blob/canary/examples/api-routes-rate-limit/pages/api/user.ts
import { LRUCache } from "lru-cache"

type Options = {
    uniqueTokenPerInterval: number;
    interval: number;
};

export default function rateLimit(options: Options) {
    const tokenCache = new LRUCache({
        max: options.uniqueTokenPerInterval,
        ttl: options.interval,
    });

    return {
        check: (limit: number, token: string) =>
            new Promise<void>((resolve, reject) => {
                const tokenCount = (tokenCache.get(token) as number[]) || [0];
                if (tokenCount[0] === 0) {
                    tokenCache.set(token, tokenCount);
                }
                tokenCount[0] += 1;

                const currentUsage = tokenCount[0];
                const isRateLimited = currentUsage >= limit;
                return isRateLimited ? reject() : resolve();
            }),
    };
}